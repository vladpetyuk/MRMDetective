from comtypes.client import CreateObject
from comtypes.automation import VARIANT
from comtypes import BSTR
from ctypes import *


class ThermoRawReaderUtility( object):

    def __init__( self, pathToFile):
        self.xrf = CreateObject( 'Xrawfile.Xrawfile')
        self.xrf.Open( szFileName = pathToFile)
        self.xrf.SetCurrentController( nControllerType = 0, nControllerNumber = 1)

    def call_Close( self):
        self.xrf.Close()

    def call_GetNumSpectra( self):
        numSpec = c_long(0)
        p_numSpec = pointer( numSpec)
        silence = self.xrf.GetNumSpectra( p_numSpec)
        return numSpec.value

    def call_RTFromScanNum( self, scanNumber):
        rt = c_double(0)
        p_rt = pointer(rt)
        silence = self.xrf.RTFromScanNum( nScanNumber=scanNumber, pdRT = p_rt)
        return rt.value
        
    def call_GetFilterForScanNum( self, scanNumber):
        filterEntry = BSTR()
        p_filterEntry = pointer( filterEntry)
        silence = self.xrf.GetFilterForScanNum( scanNumber, p_filterEntry)
        return filterEntry.value

    def call_GetInstMethod( self):
        instrMethod = BSTR()
        p_instrMethod = pointer( instrMethod)
        silence = self.xrf.GetInstMethod( 0, p_instrMethod)
        return instrMethod.value

    def call_GetFilters( self):
        filters = VARIANT()
        p_filters = pointer( filters)
        silence = self.xrf.GetFilters( p_filters, pointer(c_int()))
        return filters.value    
        
    def call_GetMassListFromScanNum( self, scanNum):
        scanNum   = c_long( scanNum)
        p_scanNum = pointer( scanNum)
        centroidPeakWidth   = c_double( 0.0) # unused ??
        p_centroidPeakWidth = pointer( centroidPeakWidth)
        massList = VARIANT()
        p_massList = pointer( massList)
        peakFlags  = VARIANT()
        p_peakFlags = pointer( peakFlags)
        arraySize = c_int(0)
        p_arraySize = pointer( arraySize)
        silence = self.xrf.GetMassListFromScanNum(
                                                   pnScanNumber = p_scanNum,
                                                   bstrFilter = '',
                                                   nIntensityCutoffType = 0,
                                                   nIntensityCutoffValue = 0,
                                                   nMaxNumberOfPeaks = 0,
                                                   bCentroidResult = 0,
                                                   pdCentroidPeakWidth = p_centroidPeakWidth,
                                                   pvarMassList = p_massList,
                                                   pvarPeakFlags = p_peakFlags,
                                                   pnArraySize = p_arraySize)
        return massList.value
        





if __name__ == '__main__':

    # path to test file
    path1 = "D:\\_Python_Projects\\tripple_quad_data_processing\\__whole_thing_bleeding_edge\\25jan2008_MP_EIF_0p1_1b.RAW"

    rawInst = ThermoRawReaderUtility( path1)    
    x = rawInst.call_GetNumSpectra()
    rt = rawInst.call_RTFromScanNum(x)
    im = rawInst.call_GetInstMethod()
    ml = rawInst.call_GetMassListFromScanNum( x)
    ft = rawInst.call_GetFilterForScanNum( x)
    ft = rawInst.call_GetFilters()
    rawInst.call_Close()

    print x, rt
    print ft
    print ml
    print im
    print ft
    print "----------"
    print ml
    








