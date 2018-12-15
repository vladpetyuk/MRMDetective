import sys
import copy
import re
import time
import numpy as N
from cls_ThermoRawReaderUtility import ThermoRawReaderUtility


'''
Old documentation

Some terms:
    scan number is the group of measurements
    sharing the same parent ion
    One SUPER scan (all transition cycle) 
    have scans (transitions only for a given parent ions)
    having MICRO scans (actual transitions)
    
Class data members are:

maxSpecNumber
    integer value
    maximum spectrum number
    
transitionData
    list of lists with 4 string values: [u'Parent', u'Center', u'Time', u'Profile']
    parent m/z
    daughter m/z
    time width used for signal integration
    Profile - Method ??? The values is not used anywhere
    
uniqParents
    list of parent ions 
    in the order they are in the instrument method
    
mrmData
    The most important member!!
    dict(dict(list(list)))
    mrmData is organized as dictionary of dictionary having lists of lists:
    primary dictionary keys are parent m/z values,
    secondary dictionary keys are daughter m/z values
    main list contain lists of 3 values: [scan, retention time, intensity]
    
dataHolderDict 
    is a template dictionary
    same structure as mrmData, but the main lists are empty
    dict(dict([]))
'''
    



class MRMDataReaderFromThermoRawFile:

    
    
    def __init__(self, filePath):
        """
        filePath is path to Thermo *.RAW file
        """
        # data objects
        self.__filePath = filePath
        self.__maxSpecNumber = None # maximum spectrum/scan number
        self.__transitionData = []  # list of lists with 4 string values: 
                                    # [u'Parent', u'Center', u'Time', u'Profile']
                                    # parent m/z
                                    # daughter m/z
                                    # time width used for signal integration
                                    # Profile - Method ??? The values is not used anywhere
        self.__uniqParentIonsTuple = [] # list of unique parents in the same order as in
                                       # instrument method 
        self.__parentIonsDictDaughterIonsTuple = {} # dictionary with keys as parent ions m/z and 
                                                    # values as list of daughet ion m/z values in the orger
                                                    # used in instrument method
        self.__mrmData_in_3x_tuples = {} # dictionary with keys as parent ions m/z and
                                        # with values as dictionary with keys as fragment ions m/z
                                        # and with values as tuple of (spectrum number, retention time, intensity)
        self.__parentIonsDictDaughterIonsDictScansTuple = {} # dictionary with keys as parent ions m/z and
                                 # with values as dictionary with keys as fragment ions m/z
                                 # and with values as tuple of spectrum numbers
        self.__parentIonsDictDaughterIonsDictRTTuple = {} # dictionary with keys as parent ions m/z and
                              # with values as dictionary with keys as fragment ions m/z
                              # and with values as tuple of retention times
        self.__parentIonsDictDaughterIonsDictIntensitiesTuple = {} # dictionary with keys as parent ions m/z and
                                       # with values as dictionary with keys as fragment ions m/z
                                       # and with values as tuple of intensities
        self.__mrmData_in_dicts = {} # dictionary with keys as parent ions m/z and
                                     # with values as dictionary with keys as fragment ions m/z
                                     # and with values as tuples of dictionaries with keys 
                                     # spectrum number, retention time, intensity
                                     # values are the corresponding actual values
                                 
        
        # functions called during initialization
        self.__checkFile()
        print "...creating dataset class from RAW file name %s" % self.__filePath
        
    
    
    def __checkFile(self):
        try:
            self.__rawFileInst = ThermoRawReaderUtility( self.__filePath)
            numSpectra = 0
            numSpectra = self.__rawFileInst.call_GetNumSpectra()
            if numSpectra == 0:
                raise Exception
        except:
            print "Cannot open the raw file!"
            sys.exit(1)
        self.__rawFileInst.call_Close()
    
    
    def __getTransitionDataFromInstrumentMethod( self):
        """
        returns list of lists with 4 string values: 
        [u'Parent', u'Center', u'Time', u'Profile']
        parent m/z
        daughter m/z
        time width used for signal integration
        Profile - Method ??? The values is not used anywhere
        """
        instrumentMethod = self.__rawFileInst.call_GetInstMethod()
        instrumentMethod = [i.strip(' ') for i in instrumentMethod.split('\r\n')]
        #TODO: Needs to be fixed! Vlad
        startString = u'Parent    Center    Width    Time   CE   Q1 PW  Q3 PW  Tube Lens'
        startString = u'Parent    Center    Width     Time    CE   Q1 PW  Q3 PW  Tube Lens'
##        startString = u'Parent    Center    Width     Time    CE   Q1 PW  Q3 PW  S-Lens'
        endString = u''
        startIdx = instrumentMethod.index( startString)
        endIdx = instrumentMethod[startIdx:].index( endString) + startIdx
        transitionData = instrumentMethod[startIdx:endIdx]
        transitionData = [i.split('  ') for i in transitionData]
        transitionData = [[j.strip(' ') for j in i if j != u''] for i in transitionData]
        transitionData = [[i[0],i[1],i[3],u'Profile'] for i in transitionData]
        return transitionData


    def __makeDataHolderDict( self):
        # I need unique ordered list of parents
        parents = [i[0] for i in self.__transitionData[1:]]
        uniqParentIonsList = []
        for p_i in range( len( parents)):
            if parents[p_i] not in parents[p_i+1:]:
                uniqParentIonsList.append( parents[p_i])
        numberOfPeptides = len( uniqParentIonsList)
        uniqParentsIonsDict = dict.fromkeys( uniqParentIonsList, [])
        for t in self.__transitionData[1:]:
            temp = uniqParentsIonsDict[t[0]]
            uniqParentsIonsDict[t[0]] = temp + [t[1]]
        for k in uniqParentsIonsDict.keys():
            temp = uniqParentsIonsDict[k]
            uniqParentsIonsDict[k] = dict.fromkeys( temp, [])
        return uniqParentIonsList, uniqParentsIonsDict



    def __makeUniqParentIonsTuple( self):
        """
        returns tuple of unique parent ion m/z values
        in the same order as in instrument method 
        """
        parents = [i[0] for i in self.__transitionData[1:]]
        uniqParentIonsTuple = ()
        for p_i in range( len( parents)):
            if parents[p_i] not in parents[p_i+1:]:
                uniqParentIonsTuple += ( parents[p_i],)
        return uniqParentIonsTuple


    def __makeParentIonsDictDaughterIonsTuple(self):
        """
        returns dictionary with keys as parent ions m/z and 
        values as tuple of daughter ion m/z values in the order
        used in instrument method
        """
        parentIonsDictDaughterIonsTuple = dict.fromkeys( self.__uniqParentIonsTuple, ())
        for transitionEntry in self.__transitionData[1:]:
            # appends the daughter ions from the read transition to the tuple 
            parentIonsDictDaughterIonsTuple[ transitionEntry[0]] += ( transitionEntry[1],) 
        return parentIonsDictDaughterIonsTuple



    def __makeParentIonsDictDaughterIonsDict(self):
        """
        returns dictionary with keys as parent ions m/z and 
        values as tuple of daughter ion m/z values in the order
        used in instrument method
        """
        parentIonsDictDaughterIonsDict = copy.deepcopy(
            self.__parentIonsDictDaughterIonsTuple)
        for parentIon in parentIonsDictDaughterIonsDict.keys():
            parentIonsDictDaughterIonsDict[parentIon] = dict.fromkeys(
                parentIonsDictDaughterIonsDict[parentIon], ())
        return parentIonsDictDaughterIonsDict



##    def __readMRMData( self):
##        """
##        The major function for extracting the data
##        """
##        data = copy.deepcopy(self.__parentIonsDictDaughterIonsDict)
####        filterPttrn = re.compile('.+\s(\d+\.\d+)@.+\s\[.+\]')
##        filterPttrn = re.compile('.+\s(\d+\.\d+)(@.+)?\s\[.+\]')
##        # get the number of used parent ions
##        numberOfPeptides = len(data.keys())
##        # run through each scan number
##        print "prints '.' per every 1000 scans"
##        for s in range( 1, self.__maxSpecNumber + 1):
##            if s % 1000 == 0:
##                print '.',
##            fltr = self.__rawFileInst.call_GetFilterForScanNum( s)
##            # get the parent m/z from filter
##            parentFromScanFilter = filterPttrn.match( fltr).group(1)
##            # the index in the list of peptide
##            parentIdx = s % numberOfPeptides
##            # if the peptide from the list is not the same as from filter
##            #     raise error!!
##            if self.__uniqParentIonsTuple[parentIdx-1] != parentFromScanFilter:
##                raise Exception
##            # retention time for the scan number
##            rt = self.__rawFileInst.call_RTFromScanNum( s)
##            # mass list as tuple of 2 tuples
##            #    1st tuple mz values of daughter ions
##            #    2nd tuple intensities
##            ml = self.__rawFileInst.call_GetMassListFromScanNum( s)
##            print ml
##            # extract and round the m/z values
##            mzs = [u'%.3f' % round(float(i),3) for i in ml[0]]
##            # fill the data list for individual transition, 
##            # that is parent/daughter m/z values combination
##            # with lists having s values: [scan, retention time, intensity]
##            for mz_index in range(len(mzs)):
##                daughter = mzs[mz_index] # just to make it clear
##                data[parentFromScanFilter][daughter] += ( ( s, rt, ml[1][mz_index]),)
##        print '\n'
##        return data









    def __readMRMData( self):
        """
        The major function for extracting the data
        """
        data = copy.deepcopy(self.__parentIonsDictDaughterIonsDict)
##        filterPttrn = re.compile('.+\s(\d+\.\d+)@.+\s\[.+\]')
        filterPttrn = re.compile('.+\s(\d+\.\d+)(@.+)?\s\[.+\]')
        # get the number of used parent ions
        numberOfPeptides = len(data.keys())
        # run through each scan number
        print "prints '.' per every 1000 scans"
        for s in range( 1, self.__maxSpecNumber + 1):
            if s % 1000 == 0:
                print '.',
            fltr = self.__rawFileInst.call_GetFilterForScanNum( s)
            # get the parent m/z from filter
            parentFromScanFilter = filterPttrn.match( fltr).group(1)
            # the index in the list of peptide
            parentIdx = s % numberOfPeptides
            # if the peptide from the list is not the same as from filter
            #     raise error!!
            if self.__uniqParentIonsTuple[parentIdx-1] != parentFromScanFilter:
                raise Exception
            # retention time for the scan number
            rt = self.__rawFileInst.call_RTFromScanNum( s)
            # mass list as tuple of 2 tuples
            #    1st tuple mz values of daughter ions in PROFILE
            #    2nd tuple intensities
            ml = self.__rawFileInst.call_GetMassListFromScanNum( s)

            # getting daughter mz values
            mzs = data[ parentFromScanFilter].keys()
            
            # fill the data list for individual transition, 
            # that is parent/daughter m/z values combination
            # with lists having s values: [scan, retention time, intensity]
            for mz_index in range(len(mzs)):
                daughter = mzs[mz_index] # just to make it clear
                intensity = self.__obtain_intensity( ml, float(daughter))
                data[parentFromScanFilter][daughter] += ( ( s, rt, intensity),)
        print '\n'
        return data




    def __obtain_intensity( self, ml, daughter):
        # 50 mDa tolerance
        # difference between mz values in profile 70 mDa
        arrMl = N.array(ml)
        arrMlInd = N.abs(arrMl[0,] - daughter) < 0.050
        return N.mean( arrMl[1,arrMlInd])
        
        








    
    def __makeParentIonsDictDaughterIonsDictScansTuple(self):
        """
        returns a dictionary with parent ions as keys and values as
        dictionary with daughter ions as keys and scan number tuple as values
        """
        returnValue = copy.deepcopy(self.__parentIonsDictDaughterIonsDict)
        for parentIon in returnValue.keys():
            for daughterIon in returnValue[parentIon].keys():
                # scan or spectrum number index is 0
                scansTuple = tuple([int(i[0]) for i in self.__mrmData_in_3x_tuples[parentIon][daughterIon]])
                returnValue[parentIon][daughterIon] = scansTuple
        return returnValue 


    def __makeParentIonsDictDaughterIonsDictRTTuple(self):
        """
        returns a dictionary with parent ions as keys and values as
        dictionary with daughter ions as keys and retention times tuple as values
        """
        returnValue = copy.deepcopy(self.__parentIonsDictDaughterIonsDict)
        for parentIon in returnValue.keys():
            for daughterIon in returnValue[parentIon].keys():
                # retention time index is 1
                rtTuple = tuple([float(i[1]) for i in self.__mrmData_in_3x_tuples[parentIon][daughterIon]])
                returnValue[parentIon][daughterIon] = rtTuple
        return returnValue 
    

    def __makeParentIonsDictDaughterIonsDictIntensityTuple(self):
        """
        returns a dictionary with parent ions as keys and values as
        dictionary with daughter ions as keys and intensities tuple as values
        """
        returnValue = copy.deepcopy(self.__parentIonsDictDaughterIonsDict)
        for parentIon in returnValue.keys():
            for daughterIon in returnValue[parentIon].keys():
                # intensity index is 2 ???? WTF
                intensityTuple = tuple([float(i[2]) for i in self.__mrmData_in_3x_tuples[parentIon][daughterIon]])
                returnValue[parentIon][daughterIon] = intensityTuple
        return returnValue 
        

    def __makeMrmData_in_dicts(self):
        """
        returns a dictionary with parent ions as keys and values as
        dictionary with daughter ions as keys and tuples of dict with keys 
        as scan number, retention time and intensity and the corresponding values
        as dictionaty values
        """
        returnValue = copy.deepcopy(self.__parentIonsDictDaughterIonsDict)
        for parentIon in returnValue.keys():
            for daughterIon in returnValue[parentIon].keys():
                for entry in self.__mrmData_in_3x_tuples[parentIon][daughterIon]:
                    currentDict = {'scan':      entry[0],
                                   'rt':        entry[1],
                                   'intensity': entry[2]}
                    returnValue[parentIon][daughterIon] += (currentDict,)
        return returnValue 



    
    #---public functions------------------------------------------
    def readFile(self):
        """
        reads the Thermo *.RAW file
        and extracts the relevant MRM data information
        """
        self.__rawFileInst = ThermoRawReaderUtility( self.__filePath)
        self.__maxSpecNumber = self.__rawFileInst.call_GetNumSpectra()
        self.__transitionData = self.__getTransitionDataFromInstrumentMethod()
        self.__uniqParentIonsTuple = self.__makeUniqParentIonsTuple()
        self.__parentIonsDictDaughterIonsTuple = self.__makeParentIonsDictDaughterIonsTuple()
        self.__parentIonsDictDaughterIonsDict  = self.__makeParentIonsDictDaughterIonsDict()
        self.__mrmData_in_3x_tuples = self.__readMRMData()
        self.__parentIonsDictDaughterIonsDictScansTuple = self.__makeParentIonsDictDaughterIonsDictScansTuple() 
        self.__parentIonsDictDaughterIonsDictRTTuple = self.__makeParentIonsDictDaughterIonsDictRTTuple()
        self.__parentIonsDictDaughterIonsDictIntensityTuple = self.__makeParentIonsDictDaughterIonsDictIntensityTuple()
        self.__mrmData_in_dicts = self.__makeMrmData_in_dicts()
        # and finally close the *.raw file
        self.__rawFileInst.call_Close()
    

    
    # accessors
    def getMaxSpecNumber(self):
        """
        returns maximum scan/spectrum number
        """
        return self.__maxSpecNumber

    
    def getParentIonTuple(self):
        """
        returns tuple of unique parent ion m/z values in the order 
        they used in the instrument method
        """
        return self.__uniqParentIonsTuple
    
    
    def getParentIonsDictDaughterIonsTuple(self):
        """
        returns dictionary with keys as parent ions m/z and 
        values as list of daughet ion m/z values in the orger
        used in instrument method
        """
        return self.__parentIonsDictDaughterIonsTuple
    

    def getMrmDataDict(self):
        """
        returns dictionary with keys as parent ions m/z and
        with values as dictionary with keys as fragment ions m/z
        and with values as tuple of [spectrum number, retention time, intensity]
        """
        return self.__mrmData_in_3x_tuples
    

#    def getMrmDataDicts(self):
#        """
#        returns a dictionary with parent ions as keys and values as
#        dictionary with daughter ions as keys and tuples of dict with keys 
#        as scan number, retention time and intensity and 
#        the corresponding values as dictionaty values
#        """
#        return self.__mrmData_in_dicts
    

    def getMrmDataScansTuple(self):
        """
        returns a dictionary with parent ions as keys and values as
        dictionary with daughter ions as keys and scan number tuple as values
        """
        return self.__parentIonsDictDaughterIonsDictScansTuple    


    def getMrmDataRTTuple(self):
        """
        returns a dictionary with parent ions as keys and values as
        dictionary with daughter ions as keys and retention time tuple as values
        """
        return self.__parentIonsDictDaughterIonsDictRTTuple    
    

    def getMrmDataIntensityTuple(self):
        """
        returns a dictionary with parent ions as keys and values as
        dictionary with daughter ions as keys and intensity tuple as values
        """
        return self.__parentIonsDictDaughterIonsDictIntensityTuple  

       


if __name__ == "__main__":
    # definetely need to desing a test
    print "start of the test"
    import os
    print os.getcwd()
    instX = MRMDataReaderFromThermoRawFile(r"9pept_1nM_inQCShew_pt25_IF_MRM_C1.RAW")
    instX.readFile()
    print "end of the test"

    
    
    
    
    
