import os
import pylab as p
from numpy import *
from cls_DaughterIon import DaughterIon
from cls_ExhaustivePeakClustering import ExhaustivePeakClustering


"""
Supposed to know 
    data objects:
        1) daughters
        2) clusters
    functions:
        1) create daughters (maintain tuple/dictionary of daughters)
        2) ask daughters to find peaks
        3) create clusters (using info from daughters) (maintain tuple/dictionary of daughters)
        4) Print all daughters for parent as a single image
        5) Print clusters. Cluster cannot print itself, since it requires context of daughter info
    
"""



class ParentIon:
    def __init__(self, 
                 parentIon, 
                 daughterIonsTuple,
                 mrmDataForParentIon,
                 maxScanNumber,
                 archiver ):
        
        self.__parentIon = parentIon #unicode
        self.__daughterIonsTuple = daughterIonsTuple #tuple of unicodes
        self.__mrmDataForParentIon = mrmDataForParentIon # dict keys-daughter m/z in unicode
                                                         # values - tuples of tuples
                                                         # of <int>scan number, 
                                                         # <float>retention time, 
                                                         # <float>intensity
        self.__maxScanNumber = maxScanNumber # integer
        self.__archiver = archiver 
        
        # initialization of new data objects
        self.__daughterIonsClss = {}
        self.__daughterIonsPeaks = {}
        
        self.__clustersTuple = ()
        self.__clustersClss = {}
        
        # run the process
        self.__daughterIonsClss  = self.__makeDaugtherIonsClss()
        self.__daughterIonsPeaks = {} # dictionary 
                                         # with daughter ions as keys
                                         # values are tuples of picked peaks
                                         # in descending order by intensity
                                         # tuple element is a dictionary 
                                         # with keys as peak attirbutes
                                         # and values as actual values


        
    def __makeDaugtherIonsClss(self):
        daughterIonsClss = {}
        for daughterIon in self.__daughterIonsTuple:
            daughterIonsClss[daughterIon] = DaughterIon( daughterIon,
                                                         self.__parentIon,
                                                    self.__mrmDataForParentIon[daughterIon],
                                                    self.__maxScanNumber,
                                                    self.__archiver)
        return daughterIonsClss
    
    
    def __extractDaughterIonsPeaks(self):
        daughterIonsPeaks = {}
        for daughterIon in self.__daughterIonsTuple:
            daughterIonsPeaks[daughterIon] = self.__daughterIonsClss[daughterIon].getPickedPeaksTuple()
        return daughterIonsPeaks

    
    
    
    #---mutators---------------------------------
    def pickPeaks(self):
        for daughterIon in self.__daughterIonsTuple:
            # so now the daughter ion classes has the peaks info ready
            self.__daughterIonsClss[daughterIon].pickPeaks()

        
    def clusterPeaks(self, method = "exhaustive"):
        # self.__daughterIonsPeaks
        # is the right variable to feed to clustering routine
        self.__daughterIonsPeaks = self.__extractDaughterIonsPeaks()
        if method == "exhaustive":
            cluteringClass = ExhaustivePeakClustering(self.__daughterIonsPeaks,
                                                      self.__daughterIonsTuple)
            cluteringClass.clusterPeaks()
            # returns clusters ordered by score
            self.__clustersTuple = cluteringClass.getClusterTuple()
                







        
    def printParentEntireTransitionsAsImage(self):
        '''
        def plotAllTransitions( mrmData, mrmPeaks, maxScan, dsetOutPath="."):
        '''
        '''
        self.__daughterIonsTuple
        '''
        # setting directory for dumping
        dsetOutPath = "."
        outDirPath = os.path.join( dsetOutPath, 'all_transitions')
#        if os.path.exists( outDirPath):
#            shutil.rmtree( outDirPath)
#        os.mkdir( outDirPath)
            
        # plotting all transitions and detected peaks
        numSubPlots = len( self.__daughterIonsTuple)
#        numSubPlots = len(mrmData[k1].keys())
        i = 1 #???
#        keys2 = mrmData[k1].keys()
#        keys2.sort()
#        self.__daughterIonsTuple is the same as keys2. Right?
#        daughterNumber = len( keys2)
#        daughterNumber is the same as numSubPlots
        for daughterIon in self.__daughterIonsTuple:
            x = array(self.__mrmDataForParentIon[daughterIon])
            sbPlt = p.subplot(str(numSubPlots)+'1'+str(i))
            p.subplots_adjust(top=0.900, bottom=0.100)
            p.grid()
            if i == 1:
                p.title('parent with %s m/z' % self.__parentIon)
            if i != numSubPlots:
                sbPlt.set_xticklabels([])
            if i == numSubPlots:
                sbPlt.set_xlabel('scan number')
            sbPlt.yaxis.major.formatter.set_powerlimits((-1,+1))

            # data itself
            p.plot(x[:,0],x[:,2])
            p.hold(True)

            # plot found peaks
            peaksTuple  = self.__daughterIonsPeaks[daughterIon]
            scans       = array([p_i["scanNumber"] for p_i in peaksTuple])
            intensities = array([p_i["intensity"]  for p_i in peaksTuple])
            p.plot(scans, intensities, 'ro', alpha=0.5)

            #
            p.xlim((0,self.__maxScanNumber))
            
            # daughter m/z        
            p.text(0.98, 0.94, '%s m/z' % daughterIon,\
                   fontsize=10,\
                   horizontalalignment='right',\
                   verticalalignment='top',\
                   transform = sbPlt.transAxes,\
                   bbox = dict(facecolor='#EEEEEE',\
                               edgecolor='#EEEEEE',\
                               pad = 5,\
                               alpha=0.8))
            i = i + 1

        imName = 'peptide_'+self.__parentIon+'.png'
        imPath = os.path.join( outDirPath, imName)
        imPath = imName
        p.savefig(str(imPath))
        p.clf()







    
    def printParentPeakClustersAsImage(self):
        # setting directory for dumping
        dsetOutPath = "."
        outDirPath = os.path.join( dsetOutPath, 'all_clusters')
#        if os.path.exists( outDirPath):
#            shutil.rmtree( outDirPath)
#        os.mkdir( outDirPath)
    
        # number of +/- loops not scans
        plotWindow = 100 
            
        # plotting all transitions and detected peaks
#        for k1 in clusters.keys():
#            if clusters[k1] is None:
#                continue
        daughters = self.__daughterIonsTuple
        clustNum = 0
        for cluster_i in self.__clustersTuple:
            # cluster_i is dictionary
            clustNum += 1
            numSubPlots = len( self.__daughterIonsTuple)
            repScanNum = cluster_i["representativeClusterScan"]
            for i in range(numSubPlots):
                #
                # extracting chromatogram data
                data = array(self.__mrmDataForParentIon[ self.__daughterIonsTuple[i]])
                loopIdx = where(data[:,0] == repScanNum)[0][0]
                xMin = data[(loopIdx-100),0]
                xMax = data[(loopIdx+100),0]                
                #
                # setting surrounding scan range to plot. defined by 1st fragment for all
                data_i = data[ (loopIdx-plotWindow) : (loopIdx+plotWindow) ,:] # surrounding data
                #
                # getting the clustered peak
                peak_i = cluster_i["data"][self.__daughterIonsTuple[i]]
                # getting the peak intensity and setting the Y limit
                if peak_i == {}:
                    # just take the data from the average scan
                    peakIntensity = data[loopIdx,2]                    
                else:
                    peakIntensity = cluster_i["data"][self.__daughterIonsTuple[i]]["intensity"]
                peakIntensity = peakIntensity * 1.25 # so there is always at least 25% left
                nextOrderMag = 10**ceil( log10( peakIntensity))
                pieces = 5 # split into 5 pieces
                pieceSize = round( 10 / 5)
                portions = pieceSize * ceil( ( peakIntensity/(nextOrderMag*0.1)) / pieceSize) 
                yMax = portions * (nextOrderMag * 0.1)
                #
                # setting the plot area
                sbPlt = p.subplot( str( numSubPlots)+'1'+str( i+1))
                p.subplots_adjust(top=0.900, bottom=0.100)
                p.grid()
                if i+1 == 1:
                    score    = round(cluster_i["clusterScore"],1)
                    p.title('parent with %s m/z, cluster %s, score %.1f' % (self.__parentIon, clustNum, score))
                if i+1 != numSubPlots:
                    sbPlt.set_xticklabels([])
                if i+1 == numSubPlots:
                    sbPlt.set_xlabel('scan number')
                sbPlt.yaxis.major.formatter.set_powerlimits((-1,+1))
                #
                # data itself
                p.plot( data_i[:,0], data_i[:,2],'b.-')
                p.ylim((0,    yMax))
                p.xlim((xMin, xMax))            
                yLim = sbPlt.get_ylim()
                xLim = sbPlt.get_xlim()
                #
                # plot found peaks in the range
                x_snb = array([j["scanNumber"] for j in self.__daughterIonsPeaks[self.__daughterIonsTuple[i]]])
                y_int = array([j["intensity"] for j in self.__daughterIonsPeaks[self.__daughterIonsTuple[i]]])
                p.plot(x_snb, y_int, 'ro', alpha=0.5)
                #
                # plot cluster line
                if peak_i == {}:
                    pass
                else:
                    p.axvline( cluster_i["data"][self.__daughterIonsTuple[i]]["scanNumber"], c='g', ls='--', alpha=0.66)
                sbPlt.set_xlim( xLim)
                sbPlt.set_ylim( yLim)
                #
                # plot some info
                if peak_i == {}:
                    pass
                else:
                    mzText          = '%s m/z'        % self.__daughterIonsTuple[i]
                    scanText        = 'scan:    %s'   % peak_i['scanNumber']
                    areaText        = 'area:    %.1e' % peak_i["area"]
                    intensityText   = 'intens.: %.1e' % peak_i["intensity"]
                    snrText         = 'SNR:     %.1f' % peak_i["snr"]            
                    legend = '\n'.join([mzText, scanText, areaText, intensityText, snrText])
                    p.text(0.75, 0.90, legend,\
                           fontsize=10,\
                           horizontalalignment='left',\
                           verticalalignment='top',\
                           transform = sbPlt.transAxes,\
                           bbox = dict(facecolor='#EEEEEE',\
                                       edgecolor='#EEEEEE',\
                                       pad = 5,\
                                       alpha=0.8))
                
            imName = 'pep_%s_clust_%s.png' % (self.__parentIon, clustNum)
            imPath = os.path.join( outDirPath, imName)
            imPath = imName
            p.savefig(str(imPath))
            p.clf()

        
        
        
        
        
        