import sys
import os
from cls_MRMDataReader import MRMDataReader
from cls_ParentIon import ParentIon
from cls_Archiver import Archiver



class DatasetControl:
    """
    What class is needed for:
    Receives the raw MRM data and processing settings file paths.
    Creates MRMDataset class instance.
    Calls peak picking.
    Calls peak grouping.
    Calls reporting of the resutls

    Arguments for initialization:
        1) MRM dataset file path
        2) Processing options file name
        
    Knows data objects:
        1) Path to MRM data file.
        2) Path to processing options file.
        3) MRMDataset object


    Knows how to:
        1) Call creation of the MRMDataset class instance.
        2) Call peak picking.
        3) Call peak grouping
        4) Call result reporting

    """
    
    
    def __init__( self, mrmDatasetFilePath, processingOptionsFilePath):
        """
        Arguments:
            1) the MRM raw data file path
            2) processing options file path
        Side Effect:
            Sets the dataFileName and processingOptionsFileName
            variables.
        """
        # initialization arguments
        self.__mrmDatasetFilePath = mrmDatasetFilePath
        self.__processingOptionsFilePath = processingOptionsFilePath
        #
        # auxiliary class for saving time when same data processed multiple times
        self.__archiver = Archiver(self.__mrmDatasetFilePath)
        #        
        # initializing some class variables
        self.__parentIonsClss = {}
        #--
        self.__maxSpecNumber = 0
        self.__parentIonsTuple = ()
        self.__parentIonsDictDaughterIonsTuple = {}
        self.__mrmDataDict = {} 
        
        
        
        print "...processing options %s" % self.__processingOptionsFilePath
        


    
    
    def __makeParentIonsClss(self):
        """
        returns a dictionary the keys are parent ion m/z 
        and values are classes for parent ions
        """
        parentIonsClss = {}
        for parentIon in self.__parentIonsTuple:
            parentIonsClss[parentIon] = ParentIon( parentIon,
                                                    self.__parentIonsDictDaughterIonsTuple[parentIon],
                                                    self.__mrmDataDict[parentIon],
                                                    self.__maxSpecNumber,
                                                    self.__archiver)
        return parentIonsClss


    def __readOrinalDataFromSource(self, source):
        self.__maxSpecNumber = source.get_max_spec_number()
        self.__parentIonsTuple = source.get_parent_ions_tuple()
        self.__parentIonsDictDaughterIonsTuple = source.get_parent_ions_dict_daughter_ions_tuple()
        self.__mrmDataDict = source.get_mrm_data_dict()



    #---calls for action---------------------------------
    def readMRMData(self):
        """
        Calls the reader class of the MRM original dataset file.
        This class will exctract all the needed information from the provided file.
        """
        if self.__archiver.is_original_data_available():
            print "...reading original data from archive"
            # read from archive
            source = self.__archiver
            self.__readOrinalDataFromSource(source)
        else:
            #use a reader       
            source = MRMDataReader( self.__mrmDatasetFilePath)
            self.__readOrinalDataFromSource(source)
            # then save into archive
            self.__archiver.save_original_data(self.__maxSpecNumber,
                                               self.__parentIonsTuple,
                                               self.__parentIonsDictDaughterIonsTuple,
                                               self.__mrmDataDict
                                               )
        #
        # create the dictionary of parent ion classes
        # I'm not sure if this is a right location
        self.__parentIonsClss  = self.__makeParentIonsClss()
        
    
    
    def pickPeaks(self):
        print "...picking peaks"
        for parentIon in self.__parentIonsTuple:
            print "w",
            self.__parentIonsClss[parentIon].pickPeaks()
        print ""


    def clusterPeaks(self):
        for parentIon in self.__parentIonsTuple:
            print "g",
            self.__parentIonsClss[parentIon].clusterPeaks()
        print ""
    
    
    def reportResults(self):
        print "...reporting results"
        for parentIon in self.__parentIonsTuple:
            print "p",
            self.__parentIonsClss[parentIon].printParentEntireTransitionsAsImage()
            self.__parentIonsClss[parentIon].printParentPeakClustersAsImage()
        print ""

    


