import shelve
import copy as cp
import pylab
from pprint import pprint as p
from numpy import *





class CWTPeakPicker:
    def __init__( self,
                  data,
                  waveletCoefficients,
                  parDict = {}):
        
        self.__data = array( data)
        
        if parDict != {}:
            wScales = arange( parDict['scale_to_scan_min'], 
                              parDict['scale_to_scan_max'], 
                              parDict['scale_to_scan_step'])[::-1]
            localNoiseWindow = parDict['rangeForNoise']
            noiseQuantile    = parDict['noiseQuantile']
            snrMin           = parDict['snr_threshold']
            scaleMin         = parDict['scale_min_threshold']
            scaleMax         = parDict['scale_max_threshold']
        else:
            #defaults        
            self.__wScales           = arange(1,50,1)[::-1]
            self.__localNoiseWindow  = 0.2
            self.__noiseQuantile     = 0.95
            self.__snrMin            = 5.0
            self.__scaleMin          = 3
            self.__scaleMax          = 50

        if waveletCoefficients is None:
            self.__wCoefs = self.__multipleScalesCWT()
        else:
            self.__wCoefs = waveletCoefficients
        self.__wLocMax     = self.__findLocalMaximaMultipleScales()
        self.__wRidges     = self.__findRidges()
        self.__wAllPeaks   = self.__extractPeaksInfo()
        self.__wFltrPeaks  = self.__filterPeaks() # __wFltrPeaks is dictionary with keys
                                                  # peakPosition, peakScale, peakMaxCoeff,
                                                  # peakArea, peakNoise, peakSNR, peakIntensity
                                                  # values are lists of corresponding values for 
                                                  # detected peaks
        self.__pickedPeaksTuple = self.__makePickedPeaksTuple() # tuple with dictionary corresponding to
                                                                # each picked peak. The keys of dictionary are
                                                                # cycleNumber, bestWaveletScale, maxCWTCoeff,
                                                                # area, noise, snr, intensity and rankByIntensity.
                                                                # The value are actual values.
                                                                # Tuple is sorted in descending order by intensity
        



    def __multipleScalesCWT(self):
        # so the wCoefs matrix is going to be 2D zero matrix 
        # number of scales by lenght of the data
        wCoefs = zeros((len(self.__wScales), self.__data.shape[0]))
        for scale_i in self.__wScales:
            # get index of the scale_i in wScales
            ind = where(self.__wScales == scale_i)[0][0]
            # computes the CWT coefficients and puts them into
            # the raw corresponding the scale_i
            wCoefs[ind,:] = self.__sinleScaleCWT( self.__data, scale_i)
        return wCoefs
    
    
    def __sinleScaleCWT(self, data, scale):
        '''
        same type of CWT used in MassSpecWavelet package of Bioconductor
        '''
        psi_xval = linspace(-8, +8, 1024) #length = 1024    
        psi = (2/3**0.5 * pi**(-0.25)) * (1 - psi_xval**2) * exp(-psi_xval**2/2)
        originalLength  = data.shape[0]
        newLength       = int(2**ceil(log2(originalLength)))
        addedLength     = newLength -originalLength
        addedData       = data[-addedLength:][::-1]
        newData         = hstack((data,addedData))
        #
        psi_xval    = psi_xval - psi_xval[0]
        dxval       = psi_xval[1] #step
        xmax        = psi_xval[-1] #max x
        #
        f = repeat(0.0, newLength)
        j = floor(arange(0,scale*xmax+1)/(scale*dxval))    
        lenWave = len(j)
        f[0:lenWave] = psi.take(j.tolist())[::-1] - mean(psi.take(j.tolist()))
        if len(f) > newLength:
            print 'scale %s is too large!' % scale
            raise Exception
        wCoefs = (1/scale**0.5) * self.__convolution(newData, f)
        wCoefs = hstack(( wCoefs[(newLength - floor(lenWave/2.0)):newLength], 
                          wCoefs[0:(newLength - floor(lenWave/2.0))] ))
        return wCoefs[0:originalLength]


    def __convolution(self, x, y):
        '''
        This function gives the same answer as in R
        convolve(x,y)
        '''
        length = x.shape[0]
        z = array([0]*length)
        x = hstack((x,x,x))
        y = hstack((z,y,z))
        res = convolve(x,y[::-1],mode='same')[floor(length/2.0):length+floor(length/2.0)]
        return res


    def __findLocalMaximaMultipleScales(self):
        wLocMax = zeros((self.__wCoefs.shape[0], self.__wCoefs.shape[1]))    
        for i in range(len(self.__wScales)):
            scale_i = self.__wScales[i]
            winSize_i = scale_i * 2 + 1
            if winSize_i < 5:
                winSize_i = 5
            wLocMax[i,:] = self.__findLocalMaximaSingleScale(self.__wCoefs[i,:], winSize_i)
        return wLocMax

    
    def __findLocalMaximaSingleScale(self, coefs, winSize):
        '''
        get local maxima for wavelet coefficients for a single scale
        '''
        #how much do I need to add to be denominated in winSize
        length =  coefs.shape[0]
        fullWindows = 1 + floor(length/winSize)
        toAddLength = fullWindows*winSize - length
        coefsExtended = hstack( ( coefs, repeat( coefs[-1], toAddLength)))
        coefsExtended.shape = (-1,winSize)
        coefsExtended = coefsExtended.T
        #index of row with msximal values
        coefsMaxima = coefsExtended.argmax(axis=0) #it just give indeces!!! axis = 0 runs over columns
        #index of columns with local maxima
        coefsColsWithLocalMaxima = where(apply_along_axis(lambda x: (max(x)>x[0]) & (max(x)>x[-1]),0,coefsExtended))[0]
        localMax = repeat(0, length)
    
        col = (coefsColsWithLocalMaxima)*winSize
        row = coefsMaxima.take(coefsColsWithLocalMaxima)
        
        localMax.put((col + row),1)
        #
        # shifting and repeating
        shift = int(floor(winSize/2))
        fullWindows = 1 + floor((length+shift)/winSize)
        toAddLength = fullWindows*winSize - (length + shift)
        coefsExtended = hstack((
                                repeat( coefs[0],  shift),
                                coefs,
                                repeat( coefs[-1], toAddLength)
                                ))
        coefsExtended.shape = (-1,winSize)
        coefsExtended = coefsExtended.T
        #index of row with msximal values
        coefsMaxima = coefsExtended.argmax(axis=0) #it just give indeces!!! axis = 0 runs over columns
        #index of columns with local maxima
        coefsColsWithLocalMaxima = where( apply_along_axis(lambda x: (max(x)>x[0]) & (max(x)>x[-1]),0,coefsExtended))[0]
        localMax.put(((coefsColsWithLocalMaxima)*winSize + coefsMaxima.take(coefsColsWithLocalMaxima)-shift),1)
        #    
        # Check whether there is some local maxima have in between distance less than winSize
        maxInd = where(localMax == 1)[0]
        selInd = where(diff(maxInd) < winSize)[0]
        if len(selInd) > 0:
            selMaxInd1 = maxInd[selInd]
            selMaxInd2 = maxInd[selInd + 1]
            temp = coefs[selMaxInd1] - coefs[selMaxInd2]
            localMax.put(selMaxInd1[temp <= 0], 0)
            localMax.put(selMaxInd2[temp >  0], 0)
    
        amp_Th = 0
        localMax[coefs < amp_Th] = 0
        return localMax


    def __findRidges( self):
        wLocMax = self.__wLocMax
        wScales = self.__wScales
        maxInd_curr = where(wLocMax[0,:] > 0)[0].tolist()
        colInd = range(1,wLocMax.shape[0]) #actually row index
        ridgeList = dict( zip( maxInd_curr, [[i] for i in maxInd_curr]))
        nLevel = len(colInd)
        for j in range(nLevel):
            col_j     = colInd[j]
            scale_j   = wScales[j+1] #scales are right from 51 to 1
            winSize_j = scale_j * 2 + 1
            selPeak_j = []
            for k in range(len(maxInd_curr)):
                ind_k   = maxInd_curr[k] #number, not an array
                start_k = max(( ind_k - winSize_j,                    0))
                end_k   = min(( ind_k + winSize_j, wLocMax.shape[1] - 1))
                ind_curr = where(wLocMax[col_j, start_k:end_k] > 0)[0] + start_k# - 1
                if ind_curr.shape[0] == 0: #actually that means break in the ridge. right?
                    ind_curr = ind_k 
                elif ind_curr.shape[0] > 1:
                    ind_curr = ind_curr[argmin(abs(ind_curr - ind_k))]
                else:
                    ind_curr = ind_curr[0] #convert from array to number
                oldEntry = ridgeList[ind_k]
                del ridgeList[ind_k]
                ridgeList[ind_curr] =  [ind_curr] + oldEntry
                selPeak_j.append(ind_curr)
            if scale_j > 1: # add non selected peaks from current level.
                '''
                I do not really care about those peaks, that observed at scale 1 only
                '''
                maxInd_next = where(wLocMax[col_j,:] > 0)[0]
                idx = [i not in selPeak_j for i in maxInd_next]
                unSelPeak_j = maxInd_next[where(idx)[0]].tolist()
                newPeak_j = dict(zip(unSelPeak_j,[[i] for i in unSelPeak_j]))
                for item in newPeak_j.items():
                    ridgeList[item[0]] = item[1]
                maxInd_curr = selPeak_j + unSelPeak_j
                maxInd_curr = dict.fromkeys(maxInd_curr).keys()
        return ridgeList


    def __extractPeaksInfo(self):
        # data, wRidges, wScales, wCoefs, localNoiseWindow, noiseQuantile

        ridges = [list(i) for i in self.__wRidges.items()]
        ridges.sort( key=lambda x: x[0]) # sorting by position. simple ridges.sort() would work as well
    
        # !!! hardcoded !!! Probably should switch to % of data points
        winSize_noise = int( self.__wCoefs.shape[1] * self.__localNoiseWindow) 
        noise = abs( self.__wCoefs[-1,])
        x = noise.copy()
        x.sort()
        globalNoise = x[int(x.shape[0] * self.__noiseQuantile)]
    
        dataLength = self.__wCoefs.shape[1]
    
        peakPosition    = []
        peakScale       = []
        peakMaxCoeff    = []
        peakArea        = []
        peakNoise       = []
        peakSNR         = []
        peakIntensity   = []
        
        for i in range(len(ridges)):
            ridge_i  = ridges[i]
            useRdg_i = ridge_i[1] #these are basically col indeces
            #
            #levels?? what rows to take?
            # 1st scalePosition as index
            # last len(ridge_i[1]) as index
            # useRdg_i are the levels, right? No I need just indeces.
            indices_i   = zip( range( 0, len(ridge_i[1])), useRdg_i)
            ridgeCoeffs_i = [self.__wCoefs[::-1,:][j] for j in indices_i]
            maxInd_i    = ridgeCoeffs_i.index(max(ridgeCoeffs_i))
            bestScale_i = self.__wScales[::-1][indices_i[maxInd_i][0]]
            peakScale.append( bestScale_i)
            peakPosition.append( indices_i[maxInd_i][1])
            peakMaxCoeff.append( max( ridgeCoeffs_i))
            peakArea.append( max( ridgeCoeffs_i) * ( bestScale_i**0.5) * ( 63.9/45.0)) # (63.9/45.0) is ratio between the total and positive areas of mexican hat wavelet
            #
            ind_i   = ridges[i][0]
            start_i = max( ind_i - winSize_noise,              0)
            end_i   = min( ind_i + winSize_noise, dataLength - 1)
            x = noise[start_i:end_i].copy()
            x.sort()
            noiseLevel_i = x[int(round(x.shape[0] * self.__noiseQuantile))-1] # quantile
            #noiseLevel_i = std(x)
            #noiseLevel_i = median(abs(x - median(x))) * 1.4826
            #
            ## Limit the minNoiseLevel to avoid the case of very low noise level, e.g., smoothed spectrum. HOW??? By using global value
            noiseLevel_i = max(globalNoise, noiseLevel_i)
            peakNoise.append(noiseLevel_i)        
            peakSNR.append(peakMaxCoeff[i]/noiseLevel_i)
            #
            # peak intensity calculation
            start_i = max( ind_i - bestScale_i,              0)
            end_i   = min( ind_i + bestScale_i, dataLength - 1)
            x = self.__data[start_i:end_i].copy()
            x.sort()
            peakIntensity_i = x[int(round(x.shape[0] * 0.99))-1] # 99th quantile
            peakIntensity.append( peakIntensity_i)
    
        wPeaks = {'peakPosition':  peakPosition,
                  'peakScale':     peakScale,
                  'peakMaxCoeff':  peakMaxCoeff,                            
                  'peakArea':      peakArea,
                  'peakNoise':     peakNoise,
                  'peakSNR':       peakSNR,
                  'peakIntensity': peakIntensity}
        return wPeaks


    def __filterPeaks(self):
        loc_wPeaks = cp.deepcopy( self.__wAllPeaks)
        loc_wPeaks = self.__filterPeaksBySNR( loc_wPeaks)
        loc_wPeaks = self.__filterPeaksByProximityToTheEdge( loc_wPeaks)
        loc_wPeaks = self.__filterPeaksByScale( loc_wPeaks)
        return loc_wPeaks
        

    def __filterPeaksBySNR( self, loc_wPeaks):
        '''
        filter our by SNR
        '''
        peakSNR  = loc_wPeaks['peakSNR']
        idx = []
        for i in range(len(peakSNR)):
            if (peakSNR[i] > self.__snrMin):
                idx.append(i)    
        for k in loc_wPeaks.keys():
            temp = loc_wPeaks[k]
            temp = [temp[i] for i in idx]
            loc_wPeaks[k] = temp
        return loc_wPeaks    


    def __filterPeaksByProximityToTheEdge( self, loc_wPeaks):
        '''
        get rid off boundary peaks
        WARNING !!! edgeSize is hardcoded I would change to %
        '''
        edgeSize = 100
        maxPosition = len(self.__data)
        peakPosition  = loc_wPeaks['peakPosition']
        idx = []
        for i in range(len(peakPosition)):
            if (peakPosition[i] > edgeSize) & (peakPosition[i] < (maxPosition - edgeSize)):
                idx.append(i)
        for k in loc_wPeaks.keys():
            temp = loc_wPeaks[k]
            temp = [temp[i] for i in idx]
            loc_wPeaks[k] = temp
        return loc_wPeaks    


    def __filterPeaksByScale( self, loc_wPeaks):
        '''
        filter our by scale range
        '''
        peakScale  = loc_wPeaks['peakScale']
        idx = []
        for i in range(len(peakScale)):
            if (peakScale[i] >= self.__scaleMin) & (peakScale[i] <= self.__scaleMax):
                idx.append(i)
        for k in loc_wPeaks.keys():
            temp = loc_wPeaks[k]
            temp = [temp[i] for i in idx]
            loc_wPeaks[k] = temp
        return loc_wPeaks    


    def __makePickedPeaksTuple(self):
        """
        returns a tuple with dictionary corresponding to
        each picked peak. The keys of dictionary are
        cycleNumber, bestWaveletScale, maxCWTCoeff, 
        area, noise, snr, intensity and rankByIntesity.
        The value are actual values.
        The tuple items are in descending order by intensity.
        """
        pickedPeaksTuple = ()
        for peak_index in range(len(self.__wFltrPeaks["peakPosition"])):
            currentPeak = {
                           "cycleNumber":      self.__wFltrPeaks["peakPosition"][peak_index],
                           "bestWaveletScale": self.__wFltrPeaks["peakScale"][peak_index],
                           "maxCWTCoeff":      self.__wFltrPeaks["peakMaxCoeff"][peak_index],
                           "area":             self.__wFltrPeaks["peakArea"][peak_index],
                           "noise":            self.__wFltrPeaks["peakNoise"][peak_index],
                           "snr":              self.__wFltrPeaks["peakSNR"][peak_index],
                           "intensity":        self.__wFltrPeaks["peakIntensity"][peak_index]
                          }
            pickedPeaksTuple += (currentPeak,)
            
        # sort peaks by intensity
        # manually encoded selection sort
        size = len( pickedPeaksTuple)
        pickedPeaksList = list( pickedPeaksTuple)
        for i in range(size-1):
            maxPeak_index = i
            for j in range(i, size):
                if (pickedPeaksList[maxPeak_index]["intensity"] < pickedPeaksList[j]["intensity"]):
                    maxPeak_index = j
            temp = pickedPeaksList[i]
            pickedPeaksList[i] = pickedPeaksList[maxPeak_index]
            pickedPeaksList[maxPeak_index] = temp
        pickedPeaksTuple = tuple(pickedPeaksList)

        # add rank by intensity value
        for i in range(len(pickedPeaksTuple)):
            pickedPeaksTuple[i]["rankByIntensity"] = i + 1
        
        return pickedPeaksTuple



    #---public fucntions
    
    
    def getPickedPeaksTuple(self):
        """
        returns a tuple with dictionary corresponding to
        each picked peak. The keys of dictionary are
        cycleNumber, bestWaveletScale, maxCWTCoeff, 
        area, noise, snr, intensity and rankByIntesity.
        The value are actual values.
        The tuple items are in descending order by intensity.
        """
        return self.__pickedPeaksTuple


    def getWaveletCoefficients(self):
        return self.__wCoefs


    
    def plot_cwt_peak_picking_results( self):
    
        dataLength = len(self.__data)
    
        # does this just to have a shortcut to access the data??
        for key in self.__wFltrPeaks.keys():
            globals()[key] = self.__wFltrPeaks[key]
            # peakPosition
            # peakScale
            # peakMaxCoeff
            # peakArea
            # peakNoise
            # peakSNR
    
        sbplot = pylab.subplot(221)
        pylab.plot(data,'b.-')
        pylab.xlim((0, dataLength))
        pylab.hold(True)
        pylab.plot(peakPosition,self.__data[peakPosition],'ro')
        pylab.plot(peakPosition, peakIntensity,'g^')    
        pylab.xlim((0, dataLength))
    
        if peakArea != []:
            sbplot = pylab.subplot(222)
    ##    pylab.vlines( mzInd,     ymin=0, ymax=peakValue,    color='blue')
            pylab.vlines( peakPosition, ymin=0, ymax=peakArea, color='red')
            pylab.xlim((0, dataLength))
    
        if peakSNR != []:
            sbplot = pylab.subplot(224)
    ##    pylab.vlines( mzInd,     ymin=0, ymax=peakSNR,    color='blue')
            pylab.vlines( peakPosition, ymin=0, ymax=peakSNR, color='red')
            pylab.xlim((0, dataLength))
    
        sbplot = pylab.subplot(223)
        pylab.imshow( self.__wCoefs, aspect='auto',interpolation='nearest')
        pylab.xlim((0, dataLength))
        pylab.yticks(arange(len(self.__wScales)),[str(i) for i in self.__wScales.tolist()])
        pylab.ylabel('scales')
        pylab.hold(True)
        peakScalesForImage = [where(self.__wScales == i)[0][0] for i in peakScale]
        pylab.plot( peakPosition, peakScalesForImage, 'kx')
        pylab.xlim((0, dataLength))
    
        pylab.show()









    """
    wPeaks = {'peakPosition':  peakPosition,
              'peakScale':     peakScale,
              'peakMaxCoeff':  peakMaxCoeff,                            
              'peakArea':      peakArea,
              'peakNoise':     peakNoise,
              'peakSNR':       peakSNR,
              'peakIntensity': peakIntensity}
    note: values are lists
    """
    """
    wScales array of scales it has gone through
    """
    """
    wCoefs is a 2D array. Size num scales X num data points.
    Values are wavelet coefficients for given scale and translation.
    """








if __name__ == '__main__':
    import time
    arc = shelve.open('cwt_test.shlv')
    dataDict = arc['upd']
    ##data = dataDict['373.329']['639.346']
    ##data = dataDict['373.329']['389.168']
    ##data = dataDict['949.959']['227.744']
    data = dataDict['949.959']['370.094']
    ##data = dataDict['546.362']['836.410']
    ##data = dataDict['598.318']['1007.468']
    data = array(data)
    data = data[:,2]

    time1 = time.clock()
    cwtPeakPicker = CWTPeakPicker( data)
    time2 = time.clock()
    print "cwt time: ", round(time2 - time1,3), " sec"
    cwtPeakPicker.plot_cwt_peak_picking_results()




