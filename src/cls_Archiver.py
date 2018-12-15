import os
import shelve




class Archiver(object):
    def __init__(self, mrmDatasetFilePath):
        self.__archiveExt = ".MRMDTC"
        dataSetName = os.path.splitext(mrmDatasetFilePath)[0]
        self.__arcFilePath = dataSetName + self.__archiveExt
        #
        # original data variables
        self.__maxSpecNumber = 0
        self.__parentIonsTuple = ()
        self.__parentIonsDictDaughterIonsTuple = ()
        self.__mrmDataDict = {}
        


    def is_original_data_available(self):
        arcFileExist = os.path.exists(self.__arcFilePath)
        if arcFileExist:
            archive = shelve.open(self.__arcFilePath)
            if archive.has_key("originalData"):
                originalDataInFile = True
            else:
                originalDataInFile = False
            archive.close()
        # make final conlusion    
        if arcFileExist and originalDataInFile:
            isAvailable = True
            self.__read_original_data()
        else:
            isAvailable = False
        return isAvailable
        
            
        

        
    def __read_original_data(self):
        """
        read original data from archive
        """
        archive = shelve.open(self.__arcFilePath)
        originalData = archive["originalData"]
        self.__maxSpecNumber = originalData["maxSpecNumber"]
        self.__parentIonsTuple = originalData["parentIonsTuple"]
        self.__parentIonsDictDaughterIonsTuple = originalData["parentIonsDictDaughterIonsTuple"]
        self.__mrmDataDict = originalData["mrmDataDict"]
        archive.close()


    def save_original_data(self, maxSpecNumber,
                                 parentIonsTuple,
                                 parentIonsDictDaughterIonsTuple,
                                 mrmDataDict
                           ):
        """
        save original data into archive
        """
        originalData = {"maxSpecNumber": maxSpecNumber,
                        "parentIonsTuple": parentIonsTuple,
                        "parentIonsDictDaughterIonsTuple": parentIonsDictDaughterIonsTuple,
                        "mrmDataDict": mrmDataDict}
        archive = shelve.open(self.__arcFilePath)
        archive["originalData"] = originalData 
        archive.close()        
        


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


        
    def is_wavelet_coefficients_data_available(self, parent, daughter):
        archive = shelve.open(self.__arcFilePath)
        key = self.__make_key_for_wavelet_coefficients_access(parent, daughter)
        if archive.has_key(key):
            isAvailable = True
        else:
            isAvailable = False
        archive.close()
        return isAvailable
    
    
    def get_wavelet_coefficients_data(self, parent, daughter):
        archive = shelve.open(self.__arcFilePath)
        key = self.__make_key_for_wavelet_coefficients_access(parent, daughter)
        data = archive[key]
        archive.close()
        return data
    

    def set_wavelet_coefficients_data(self, parent, daughter, data):
        archive = shelve.open(self.__arcFilePath)
        key = self.__make_key_for_wavelet_coefficients_access(parent, daughter)
        archive[key] = data
        archive.close()

    
    def __make_key_for_wavelet_coefficients_access(self, parent, daughter):
        return str("waveletCoefficients" + parent + daughter)
    