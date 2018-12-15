import sys
from cls_DatasetControl import DatasetControl



class MainApp:
    """
    Arguments for constructor:
        MRM data file path
        (some other options in the future)
    Knows data objects:
        Path to MRM data file.
        Path to processing options file.
    Knows how to:
        Check the validity and assign arguments.
        Print wellcome header.
        Create DatasetProcessControl class instance.
        Ask DatasetProcessControl to find peaks.
        Ask DatasetProcessControl to group peaks.
        Ask DatasetProcessControl to report results.
    """
    def __init__(self, args):
        self.__args = args
        self.__mrmDatasetFilePath = None
        self.__processingOptionsFilePath = "MRMDetectiveDefaultOptions.xml"
        self.__checkAndAssignArgs()
        
    def __checkAndAssignArgs(self):
        if (len(self.__args) == 0):
            self.__printWellcomeHeader()
            sys.exit(1)            
        else:
            self.__mrmDatasetFilePath = self.__args[0]
            
    def __printWellcomeHeader(self):
        print "\n"\
              "No arguments received!\n"\
              "Command is runMRMDetective.py file\n"

    def run(self):
        # all the dataset reading and data structe creation happens during initialization
        __datasetProcessingControl = DatasetControl( self.__mrmDatasetFilePath,
                                                     self.__processingOptionsFilePath)
        __datasetProcessingControl.readMRMData()
        __datasetProcessingControl.pickPeaks()
        __datasetProcessingControl.clusterPeaks()
        __datasetProcessingControl.reportResults()

    def printReport(self):
        print "...processing is done"




def main():
    """
    Creates the MainApp class instance.
    Run process.
    """
    app = MainApp( sys.argv[1:])
    app.run()
    app.printReport()


if __name__ == "__main__":
    main()
    
