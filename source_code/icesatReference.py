# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 11:18:58 2019

@author: malonzo
"""

# # Import modules
# import os
# import numpy as np
# import time as runTime
# import ntpath
# from osgeo import gdal, ogr
# from getAtlMeasuredSwath_auto import getAtlMeasuredSwath
# from icesatIO import (readLas, writeLas, readGeoidFile, readDEMepsg, formatDEM)
# from icesatUtils import (getCoordRotFwd, getCoordRotRev, getLatLon2UTM,
#                          getGeoidHeight, identifyEPSG, transform, getRaster)
# from icesatUtils import getNameParts

import numpy as np
import time

from getAtlTruthSwath_auto import getAtlTruthSwath, getTruthFilePaths

from getMeasurementError_auto import getMeasurementError, offsetsStruct
# from icesatReader import read_atl03_geolocation
# from icesatReader import match_atl_to_atl03
from icesatUtils import indexMatch, getGeoidHeight
from icesatIO import getTruthHeaders, readGeoidFile

def estimate_segment_id_legacy(geolocation, gt, atlTruthStructLegacy):

    # geolocation = read_atl03_geolocation(atl03filepath, gt)
    
    # geolocation, rotation_data, epsg = match_atl_to_atl03(geolocation, atl03) 
    
    # Find closest (within 10 meters) 
    seg_id = np.array(geolocation.segment_id)
    seg_at = np.array(geolocation.alongtrack)
    truth_at = atlTruthStructLegacy.alongTrack.flatten()
    index = indexMatch(seg_at, truth_at)
    datalen = len(seg_at)
    index[index >= datalen] = (datalen - 1)
    include = np.zeros(len(truth_at))
    seg_at_comp = np.array([seg_at[x] for x in index])
    seg_id_truth = np.array([seg_id[x] for x in index])
    diff = np.abs(seg_at_comp - truth_at)
    
    include[diff < 12] = 1
    
    return seg_id_truth, include
    
def legacy_get_truth_swath(atl03legacy, rotationData, truthSwathDir, truthFileType, 
                           outFilePath, buffer = 50, useExistingTruth = False,
                           createTruthFile = True):

    # Get input truth file(s)
    truthFilePaths = getTruthFilePaths(truthSwathDir, truthFileType, logFileID=False)
              
    # Get truth file header info
    if(not(useExistingTruth)):
        truthHeaderDF = getTruthHeaders(truthFilePaths, truthFileType, logFileID=False)
    else:
        truthHeaderDF = False
    # endIf
    
    # Call getAtlTruthSwath
    print('RUNNING getAtlTruthSwath...\n')
    atlTruthData = getAtlTruthSwath(atl03legacy, rotationData, 
                                    truthHeaderDF, truthFilePaths,
                                    buffer, outFilePath, createTruthFile, 
                                    truthFileType, useExistingTruth, 
                                    logFileID=False)

    # Reclassify anything 
    atlTruthData.classification[atlTruthData.classification == 3] = 4
    atlTruthData.classification[atlTruthData.classification == 5] = 4
    
    return atlTruthData
    
def legacy_get_meas_error(atl03legacy, atlTruthData, rotationData, outFilePath, truthgroundclass=2):
    # truthgroundclass = 2
    offsetsCrossTrackBounds = np.array([-48, 48])      # Cross-track limits
    offsetsAlongTrackBounds = np.array([-48, 48 ])      # Along-track limits 
    offsetsRasterResolutions = np.array([8, 4, 2, 1])  # Step-down resolutions
    offsetsUseVerticalShift = False    # Option to use a vertical shift
    offsetsVerticalShift = 0   # Vertical shift to use if above set to True
    #measClassFilter = 1 # Meas Classes (0 = Unclass, 1 = Ground)
    filterData = truthgroundclass 
    useMeasSigConf = False # Use measured signal confidence
    offsets = offsetsStruct(offsetsCrossTrackBounds, offsetsAlongTrackBounds, 
                            offsetsRasterResolutions, offsetsUseVerticalShift, 
                            offsetsVerticalShift )
    createMeasCorrFile = True # Option to create ouput measured corrected .las
    makePlots = False         # Option to make output plots
    showPlots = False         # Option to show output plot windows
    refHeightType = 'HAE'

    atlCorrections = getMeasurementError(atl03legacy, atlTruthData, 
                                         refHeightType, 
                                         rotationData, outFilePath, 
                                         useMeasSigConf, filterData, offsets, 
                                         createMeasCorrFile, makePlots, 
                                         showPlots)

    return atlCorrections
    
def apply_offset_legacy(atl03struct, atlCorrections):
    atl03struct.df.alongtrack = atl03struct.df.alongtrack +\
        atlCorrections.alongTrackBounds
    atl03struct.df.crosstrack = atl03struct.df.crosstrack +\
        atlCorrections.crossTrackBounds
    atl03struct.df.h_ph = atl03struct.df.h_ph + atlCorrections.verticalShift



# Unit test
if __name__ == "__main__":
    import os
    # from getAtlTruthSwath_auto import getAtlTruthSwath
    # from getMeasurementError_auto import getMeasurementError, offsetsStruct
    from icesatReader import get_atl03_struct
    from icesatReader import convert_atl03_to_legacy
    from icesatReader import get_atl_alongtrack
    import pandas as pd

    if os.name == 'nt':
        basepath03 = 'Z:/data/release/002/ATL03_r002/Finland/'
        basepath08 = 'Z:/data/release/002/ATL08_r002/Finland/'
    else:
        basepath03 = '/laserpewpew/data/release/002/ATL03_r002/Finland/'
        basepath08 = '/laserpewpew/data/release/002/ATL08_r002/Finland/'

    atl03file = 'ATL03_20181118120428_07770103_002_01.h5'
    atl08file = 'ATL08_20181118120428_07770103_002_01.h5'

    # Inputs
    atl03filepath = basepath03 + atl03file
    atl08filepath = basepath08 + atl08file
    gt = 'gt1r'
    
    header_file_path =\
        '/LIDAR/server/USERS/eric/1_experiment/Finland_HeaderData.mat'
        
    kml_bounds_txt1 = '/LIDAR/server/USERS/eric/2_production/kmlBounds.txt'
   
    print('Generate ATL03 Struct')
    atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath, 
                             epsg = '32635', kml_bounds_txt = kml_bounds_txt1, 
                             header_file_path = header_file_path)    
    
    atl03.df = atl03.df[atl03.df['time'] < 12]
    
    df, rotation_data = get_atl_alongtrack(atl03.df)
    
    atl03.df = df
    atl03.rotationData = rotation_data

    print('Convert Struct to Legacy')    
    atl03legacy, rotationData = convert_atl03_to_legacy(atl03)
    
    # Legacy Truth Swath Inputs
    buffer = 50                 # Distance in cross-track (meters) around ATL03 track to look for truth data 
    useExistingTruth = False     # Option to use existing truth data if it exists
    # truthSwathDir = 'N:/data/TanDEMX//MississippiValley_TDR_DSM/Indiana/TDT_N37W086_02_DEM.TIF'
    truthSwathDir = '/laserpewpew/data/validation/data/Finland/LAS_UTM'
    outFilePath = '/LIDAR/server/USERS/eric/1_experiment/ecosystem'
    createTruthFile = True      # Option to create output truth .las file

    
    timeStart = time.time()
    
    # Call getAtlMeasuredSwath
    # print('RUNNING getAtlMeasuredSwath...\n')
    # atl03Data, atl08Data, headerData, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, gtNum, trimInfo, createAtl03LasFile, createAtl03KmlFile, createAtl08KmlFile, createAtl03CsvFile, createAtl08CsvFile)
        
    # Call getAtlTruthSwath
    print('Run Legacy Truth Swath')
    truthFileType = '.las'
    atlTruthData = legacy_get_truth_swath(atl03legacy, atl03.rotationData, 
                                          truthSwathDir, truthFileType, 
                           outFilePath, buffer = 50, useExistingTruth = False,
                           createTruthFile = True)
    
    # End timer
    
    geoidDataFile = '/laserpewpew/data/validation/geoid/FIN2005N00/geoidFin2005N00_latlon.mat'
    geoidData = readGeoidFile(geoidDataFile)
    atlTruthData = getGeoidHeight(geoidData,atlTruthData)
    timeEnd = time.time()
    timeElapsedTotal = timeEnd - timeStart
    timeElapsedMin = np.floor(timeElapsedTotal / 60)
    timeElapsedSec = timeElapsedTotal % 60
        
    # Print completion message
    print('   Script Completed in %d min %d sec.' % (timeElapsedMin, 
                                                      timeElapsedSec))
    print('\n')

    atlCorrections = legacy_get_meas_error(atl03legacy, atlTruthData, 
                                           rotationData, outFilePath)
    
    ##Filter Truth Data by Include
    alongtrack = atlTruthData.alongTrack.flatten()
    crosstrack = atlTruthData.crossTrack.flatten()
    z = atlTruthData.z.flatten()
    easting = atlTruthData.easting.flatten()
    northing = atlTruthData.northing.flatten()
    classification = atlTruthData.classification.flatten()
    intensity = atlTruthData.intensity.flatten()
    

    df_truth = pd.DataFrame(z,columns=['z'])
    df_truth = pd.concat([df_truth,pd.DataFrame(
        crosstrack,columns=['crosstrack'])],axis=1)
    df_truth = pd.concat([df_truth,pd.DataFrame(
        alongtrack,columns=['alongtrack'])],axis=1)
    df_truth = pd.concat([df_truth,pd.DataFrame(
        easting,columns=['easting'])],axis=1)
    df_truth = pd.concat([df_truth,pd.DataFrame(
        northing,columns=['northing'])],axis=1)
    df_truth = pd.concat([df_truth,pd.DataFrame(
        classification,columns=['classification'])],axis=1)

    
