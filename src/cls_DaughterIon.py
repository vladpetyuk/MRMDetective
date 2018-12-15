import sys
import os
import shelve
from cls_CWTPeakPicker import CWTPeakPicker


"""
Supposed to know 
    data objects:
        1) peaks
    functions:
        1) pick peaks 
    
"""



class DaughterIon:
    def __init__(self,
                 daughterIon,
                 parentIon,
                 mrmDataForDaughterIon,
                 maxScanNumber,
                 archiver):
        
        self.__daughterIon = daughterIon # unicode
        self.__parentIon = parentIon #unicode
        self.__mrmDataForDaughterIon = mrmDataForDaughterIon #tuples of tuples
                                                             # of <int>scan number, 
                                                             # <float>retention time, 
                                                             # <float>intensity
        self.__maxScanNumber = maxScanNumber # integer
        self.__archiver = archiver
        
        # new data objects
        self.__pickedPeaksTuple = ()
        self.__pickedPeaksClss = {}
        
  
    def __pickPeaksByCWT(self):
        """
        extracts ion intensities from __mrmDataForDaughterIon
        starts CWTPeakPicker class
        put the peaks into __pickedPeaksTuple
        """
        # OPEN ARCHIVE
        #
        # check if the waveletCoefficients is in there
        #  name mangling
        if self.__archiver.is_wavelet_coefficients_data_available(self.__parentIon, self.__daughterIon):
            waveletCoefficientsArray = self.__archiver.get_wavelet_coefficients_data(self.__parentIon, self.__daughterIon)
        else:
            waveletCoefficientsArray = None
        # THE CORE
        #===============================================================================
        data = tuple([i[2] for i in self.__mrmDataForDaughterIon]) 
        cwtPeakPicker = CWTPeakPicker( data, waveletCoefficientsArray)
        self.__pickedPeaksTuple = cwtPeakPicker.getPickedPeaksTuple()
        #===============================================================================
        # END CORE
        #
        if waveletCoefficientsArray is None:
            waveletCoefficientsArray = cwtPeakPicker.getWaveletCoefficients()            
            self.__archiver.set_wavelet_coefficients_data(self.__parentIon, 
                                                          self.__daughterIon, 
                                                          waveletCoefficientsArray)
               

        
    def __addExtraInfoToPeaks(self):
        """
        adds retention time and scan number to peaks
        """
        for peak in self.__pickedPeaksTuple:
            cycleNumber = peak['cycleNumber']
            peak['scanNumber'] = self.__mrmDataForDaughterIon[cycleNumber][0]
            peak['rt']         = self.__mrmDataForDaughterIon[cycleNumber][1]
        
        
    #---mutators------------
    def pickPeaks(self, method='CWT'):
        """
        redirects to the right peak finding procedure
        Options include:
        1) CWT
        2) though centroids???
        No arguments accepted so far
        """
        if method == 'CWT':
            self.__pickPeaksByCWT()
        else:
            print "unknown method"
            sys.exit(1)
        self.__addExtraInfoToPeaks()

    #---accessors------------
    def getPickedPeaksTuple(self):
        return self.__pickedPeaksTuple
    
        
        