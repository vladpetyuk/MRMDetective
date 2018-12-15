import os
import shelve
from cls_MRMDataReaderFromMzXmlFile import MRMDataReaderFromMzXmlFile
from cls_MRMDataReaderFromThermoRawFile import MRMDataReaderFromThermoRawFile





class MRMDataReader:
    """
    The class is just for reading the MRM data in general format.
    The process is passed further to specific file type readers.
    """
    
    def __init__(self, mrmDatasetFilePath):
        """
        mrmDatasetFilePath is the path to the file with original MRM data
        """
        # dataset file path
        self.__mrmDatasetFilePath = mrmDatasetFilePath
        dataSetName = os.path.splitext(self.__mrmDatasetFilePath)[0]
        extension = os.path.splitext(self.__mrmDatasetFilePath)[1].upper()
        arcFilePath = dataSetName + ".MRMORI"        
        #
        # data to extract:
        self.__maxSpecNumber = 0
        self.__parentIonsTuple = ()
        self.__parentIonsDictDaughterIonsTuple = {}
        self.__mrmDataDict = {}
        

        if extension == ".RAW":
            # checking if file is good or not happens at the initialization step
            self.__specificReader = MRMDataReaderFromThermoRawFile( self.__mrmDatasetFilePath) 
        elif extension == ".MZXML":
            self.__specificReader = MRMDataReaderFromMzXmlFile( self.__mrmDatasetFilePath)
        else:
            print "Unknown extention of MRM dataset file!"
            sys.exit(1)
        #
        # go ahead if everything is fine
        self.__specificReader.readFile()
        #
        # get values for the above data objects
        self.__extractDataFromSpecificReader()
        #
#        # save into archive so it is easy to get next time
#        self.__saveOriginalDataIntoArchive(arcFilePath)
                







    
    def __extractDataFromSpecificReader(self):
        self.__maxSpecNumber = self.__specificReader.getMaxSpecNumber()
        self.__parentIonsTuple = self.__specificReader.getParentIonTuple()
        self.__parentIonsDictDaughterIonsTuple =\
            self.__specificReader.getParentIonsDictDaughterIonsTuple()
        self.__mrmDataDict = self.__specificReader.getMrmDataDict()
        

    #---public functions-----   
    # accessors
    def get_max_spec_number(self):
        """
        returns maximum scan/spectrum number
        """
        return self.__maxSpecNumber

    
    def get_parent_ions_tuple(self):
        """
        returns tuple of unique parent ion m/z values in the order 
        they used in the instrument method
        """
        return self.__parentIonsTuple
    
    
    def get_parent_ions_dict_daughter_ions_tuple(self):
        """
        returns dictionary with keys as parent ions m/z and 
        values as list of daughet ion m/z values in the orger
        used in instrument method
        """
        return self.__parentIonsDictDaughterIonsTuple
    

    def get_mrm_data_dict(self):
        """
        returns dictionary with keys as parent ions m/z and
        with values as dictionary with keys as fragment ions m/z
        and with values as tuple of tuples (spectrum number, retention time, intensity)
        """
        return self.__mrmDataDict
    

#    def getMrmDataDicts(self):
#        """
#        returns a dictionary with parent ions as keys and values as
#        dictionary with daughter ions as keys and tuples of dict with keys 
#        as scan number, retention time and intensity and 
#        the corresponding values as dictionaty values
#        """
#        return self.__specificReader.getMrmDataDicts()
    

#    def getMrmDataScansTuple(self):
#        """
#        returns a dictionary with parent ions as keys and values as
#        dictionary with daughter ions as keys and scan number tuple as values
#        """
#        return self.__specificReader.getMrmDataScansTuple()
#
#
#    def getMrmDataRTTuple(self):
#        """
#        returns a dictionary with parent ions as keys and values as
#        dictionary with daughter ions as keys and retention time tuple as values
#        """
#        return self.__specificReader.getMrmDataRTTuple()    
#    
#
#    def getMrmDataIntensityTuple(self):
#        """
#        returns a dictionary with parent ions as keys and values as
#        dictionary with daughter ions as keys and intensity tuple as values
#        """
#        return self.__specificReader.getMrmDataIntensityTuple()



