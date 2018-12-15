from pprint import pprint as pp
#from numpy import *
from copy import deepcopy






class ExhaustivePeakClustering:
    
    REQURED_MIN_AGREEMENT = 0.9 # 0.4 another testing option
    
    def __init__(self, daughterIonsPeaks, daughterIonsTuple):
        self.__daughterIonsPeaks = daughterIonsPeaks # dictionary 
                                                     # with daughter ions as keys
                                                     # values are tuples of picked peaks
                                                     # in descending order by intensity
                                                     # tuple element is a dictionary 
                                                     # with keys as peak attirbutes
                                                     # and values as actual values
                                                     # {
                                                     #  'daughter1':({},{},),
                                                     # }
                                                     
        self.__daughterIonsTuple = daughterIonsTuple # tuple of daughter ion m/z in 
                                                     # the same order as on intrument
                                                     # ('daughter1',,,)
        
        # clustering results is a tuple of dictionaries
        # like: 'score':val, 'size':val, 'daughtersPresence':(0|1,), 'peaks':{'':{}}
        # where peaks value is dictionary with keys of present daughter ions
        # and values are dictionaries representing peaks
        self.__clusterTuple = ()


    def __generateCombinations( self, x):
        """
        Argument:
            x is list1 of lists2 (or tuple of tuples)
        Return:
            tuple of length of multiple of
            lengths of lists2. Each element of the tuple is
            tuple of lenght(list1) and is combination 
            of elements from lists2
        """
        combs = [] # this list contains all the combinations
        def nestedLoops( combination, remainedLists):
            """
            beatiful recursive nested for loop implementation
            """
            if len(remainedLists) == 1:
                for item in remainedLists[0]:
                    combs.append( combination + [item])
            else:
                for item in remainedLists[0]:
                    nestedLoops( combination + [item], remainedLists[1:])
        nestedLoops( [], x)
        # convert to tuples (for safety reasons)
        combs = tuple([tuple(i) for i in combs])
        return combs


    def __prefilterCombinationCodesBySize(self, __combinationCodes):
#        combs = deepcopy(list(self.__combinationCodes))
        minSize = len(__combinationCodes[0]) * self.REQURED_MIN_AGREEMENT
        filteredCombs = [i for i in __combinationCodes if sum(i) > minSize]
        return tuple( filteredCombs)


    def __generateCombinationsForCode(self, code):
        peaksForCombinations = ()
        for i in zip(self.__daughterIonsTuple, code):
            if i[1] == 1:
                peaksForCombinations += (self.__daughterIonsPeaks[i[0]],)
        tplCombinationTuplesForCode = self.__generateCombinations(peaksForCombinations)
        tplCombinationDictsForCode = ()
        for peakCombDataTuple in tplCombinationTuplesForCode:
            peakCombDataList = list(peakCombDataTuple)
            peakCombDataList.reverse()
            peakCombDataDict = {}
            # get individual tuple of dictionary/peaks
            for i in zip(self.__daughterIonsTuple, code):
                if   i[1] == 0:
                    peakCombDataDict[i[0]] = {}
                elif i[1] == 1:
                    peakCombDataDict[i[0]] = peakCombDataList.pop()
                else:
                    raise Exception, "How is that possible? :)"
            tplCombinationDictsForCode += ({"data":peakCombDataDict},)
        return tplCombinationDictsForCode


    def __calcultateMaxCycleNumberDifferenceForPeakCombination(self, comb):
        cycleNumbers = []
        for daughter_i in comb["data"].keys():
            if comb["data"][daughter_i] != {}:
                cycleNumbers.append( comb["data"][daughter_i]["cycleNumber"])
        cycleNumbers.sort()
        maxCycleNumberDifference = cycleNumbers[-1] - cycleNumbers[0]
        return maxCycleNumberDifference
    
    
    def __isGoodCombination(self, comb, code):
        """
        return true is max difference between cycle numbers is
        no more than (number of involved daughters - 1)
        """
        TEST_ONLY = 10 # 1 means neighboring only
        isGood = False
        if comb["maxCycleNumberDifference"] <= (sum(code) - 1) * TEST_ONLY :
            isGood = True
        return isGood
    
    def __findCodeSupersets(self, code, codesTuple):
        superset = ()
        for codesTuple_i in codesTuple:
            isLarger = sum(codesTuple_i) > sum(code)
            allOnesCovered = len([-1 for i in zip(codesTuple_i,code) if i[0]-i[1] < 0]) == 0
            if isLarger and allOnesCovered:
                 superset += (codesTuple_i,)
        return superset
             

    def __testIfSuperset(self, smallPeakCluster, largePeakCluster):
        """
        smallPeakSet, largePeakSet are peak cluster dictionaries
        We also assume, cycleNumber is unique identifies of a peak in a
        given daughter
        """
        isSuperset = True
        daughters = smallPeakCluster["data"].keys()
        for daughter_i in daughters:
            peakFromSmall = smallPeakCluster["data"][daughter_i]
            if peakFromSmall != {}:
                peakFromLarge = largePeakCluster["data"][daughter_i]
                # check if the peaks are the same
                if peakFromSmall["cycleNumber"] != peakFromLarge["cycleNumber"]:
                    isSuperset = False
                    break
        return isSuperset 

        
    def __computerClusterScore(self, cluster):
        """
        returns average SNR for cluster
        """
        numAlignedTransitions = len(cluster["presence"])
        sumOfSNR = 0
        for daughter in cluster["data"].keys():
            peak = cluster["data"][daughter]
            if peak != {}:
                sumOfSNR += peak["snr"]
        averageSNR = sumOfSNR/numAlignedTransitions
        return averageSNR 
    
    
    def __getRepresentativeClusterScan(self, cluster):
        """
        returns representative scan number for cluster
        """
        repPeak = {"intensity":-1}
        for daughter in cluster["data"].keys():
            peak = cluster["data"][daughter]
            if (peak != {}) and (peak["intensity"] > repPeak["intensity"]):
                    repPeak = peak
        return repPeak["scanNumber"]
        

    
    #---mutators---
    def clusterPeaks(self):
        """
        start clustering.
        put some values into __clusterTuple
        """
        # == 1 ==
        # encode allowed combinations
        # as 0 and 1 for daughter ions
        # somehow has to take into account those 
        # situations when no peaks found for a transition
        self.__combinationCodes = () # stores tuples of combinations (0|1)
        codesPossibilities = ((1,0),) * len(self.__daughterIonsTuple)
        self.__combinationCodes = self.__generateCombinations( codesPossibilities)
        self.__combinationCodes = self.__prefilterCombinationCodesBySize( self.__combinationCodes)
        
        
        # == 2 ==
        # generate all actual combinations 
        # for the given combination code
        self.__combinations = {} # key = combination code, value = tuple of peaks
                                 # tuple elements are dictionaries representing peak groups/combinations
                                 # keys are daughter values are peaks(dictionaries) if
                                 # peaks are missing it is {}
        for code_i in self.__combinationCodes:
             self.__combinations[code_i] = self.__generateCombinationsForCode(code_i)
        
        # == 3 ==
        # generate some distance metric for all actual combinations
        # for example as maximum difference between cycle numbers
        #
        # code_i is a tuple
        for code_i in self.__combinations.keys():
            # comb_i is dictionary
            for comb_i in self.__combinations[code_i]:
                comb_i["presence"] = code_i
                comb_i["maxCycleNumberDifference"] = self.__calcultateMaxCycleNumberDifferenceForPeakCombination(comb_i)
                
        
        # == 4 ==
        # filter based on generated score 
        # (may take into account cluster size itself)
        self.__filteredCombinations = {} # key = combination code, value = tuple of peaks
                                         # tuple elements are dictionaries representing peak groups/combinations
                                         # keys are daughter values are peaks(dictionaries) if
                                         # peaks are missing it is {}
        # code_i is a tuple                                                 
        for code_i in self.__combinations.keys():
            self.__filteredCombinations[code_i] = ()
            # comb_i is dictionary
            for comb_i in self.__combinations[code_i]:
                if self.__isGoodCombination(comb_i, code_i):
                    self.__filteredCombinations[code_i] += (comb_i,)


        
        
        # == 5 ==
        # evaluate if some clusters that can be viewed as 
        # subsets of larger clusters
        for code_i in self.__filteredCombinations.keys():
            # gotta figure out which codes make sense to search
            codesToSearch = self.__findCodeSupersets(code_i, self.__combinationCodes)
            clusterTuple = self.__filteredCombinations[code_i]
            for cluster_i in clusterTuple:
                supersetFound = False
                cluster_i["unique"] = True
                for codesToSearch_i in codesToSearch:
                    if supersetFound:
                        break
                    supersetsClusterTuple = self.__filteredCombinations[codesToSearch_i]
                    for supersetsClusterTuple_i in supersetsClusterTuple:
                        isSuperset = self.__testIfSuperset(cluster_i, supersetsClusterTuple_i)
                        if isSuperset:
                            supersetFound = True
                            cluster_i["unique"] = False
                            break

        # == 6 ===
        # remove sub-clusters
        self.__clustersDict = {} # key = combination code, value = tuple of peaks
                                         # tuple elements are dictionaries representing peak groups/combinations
                                         # keys are daughter values are peaks(dictionaries) if
                                         # peaks are missing it is {}
        # code_i is a tuple                                                 
        for code_i in self.__filteredCombinations.keys():
            self.__clustersDict[code_i] = ()
            # comb_i is dictionary
            for comb_i in self.__filteredCombinations[code_i]:
                if comb_i["unique"]:
                    del comb_i["unique"]
                    self.__clustersDict[code_i] += (comb_i,)         
        
        
        # == 7 ==
        # generate some confidence score for cluster
        # like average SNR, multiple of SNR, multiple of ranks.
        # sum of ranks or sum of ranks divided by factorial of number of
        # involved transitions
        #
        # code_i is a tuple                                                 
        for code_i in self.__clustersDict.keys():
            # comb_i is dictionary
            for comb_i in self.__clustersDict[code_i]:
                comb_i["clusterScore"] = self.__computerClusterScore(comb_i)          
        
        
        # == 8 ==
        # append passed clusters to cluster tuple
        # code_i is a tuple                                                 
        for code_i in self.__clustersDict.keys():
            # comb_i is dictionary
            for comb_i in self.__clustersDict[code_i]:
                comb_i["representativeClusterScan"] = self.__getRepresentativeClusterScan(comb_i)
                self.__clusterTuple += (comb_i,)   

    
    #---accessors---
    def getClusterTuple(self):
        """
        retrieve the results
        """
        return self.__clusterTuple



if __name__ == "__main__":
    import shelve
    arc = shelve.open("for_clustering.shlv")
    daughterIonsTuple = arc["daughterIonsTuple"]
    daughterIonsPeaks = arc["daughterIonsPeaks"]
    clusterer = ExhaustivePeakClustering(daughterIonsPeaks, daughterIonsTuple)
    clusterer.groupPeaks()
    pp(clusterer.getClusterTuple())
    
    
    