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
import pandas as pd
import datetime

# import time
# from icesatUtils import identifyEPSG
from icesatUtils import getCoordRotFwd
from icesatUtils import transform

# from getAtlTruthSwath_auto import getAtlTruthSwath
from getAtlTruthSwath_auto import getTruthFilePaths

# from getMeasurementError_auto import getMeasurementError, offsetsStruct
# from icesatReader import read_atl03_geolocation
# from icesatReader import match_atl_to_atl03
# from icesatUtils import indexMatch, getGeoidHeight
# from icesatIO import getTruthHeaders, readGeoidFilem, readLas
from icesatIO import getTruthHeaders
from icesatIO import readLasHeader
from icesatIO import formatDEM
# from icesatIO import loadTifFile
# from icesatIO import readLas
# from icesatUtils import identifyEPSG
# from icesatUtils import getCoordRotFwd
# from icesatUtils import transform
# from icesatUtils import getGeoidHeight
from icesatUtils import getCoordRotRev
# from icesatUtils import superFilter
# from icesatUtils import getRaster
# from icesatUtils import getUTM2LatLon
from laspy.file import File
from icesatUtils import indexMatch
# from icesatUtils import indexMatch
from ace import ace
from getMeasurementError import getMeasurementError
from icesatReader import get_atl_alongtrack
from icesatCalVal import perfect_classifier


# def estimate_segment_id_legacy(geolocation, gt, atlTruthStructLegacy):

#     # geolocation = read_atl03_geolocation(atl03filepath, gt)
    
#     # geolocation, rotation_data, epsg = match_atl_to_atl03(geolocation, atl03) 
    
#     # Find closest (within 10 meters) 
#     seg_id = np.array(geolocation.segment_id)
#     seg_at = np.array(geolocation.alongtrack)
#     truth_at = atlTruthStructLegacy.alongTrack.flatten()
#     index = indexMatch(seg_at, truth_at)
#     datalen = len(seg_at)
#     index[index >= datalen] = (datalen - 1)
#     include = np.zeros(len(truth_at))
#     seg_at_comp = np.array([seg_at[x] for x in index])
#     seg_id_truth = np.array([seg_id[x] for x in index])
#     diff = np.abs(seg_at_comp - truth_at)
    
#     # include[diff < 12] = 1
    
#     return seg_id_truth, include
    
# def legacy_get_truth_swath(atl03legacy, rotationData, truthSwathDir, truthFileType, 
#                            outFilePath, buffer = 50, useExistingTruth = False,
#                            createTruthFile = True):

#     # Get input truth file(s)
#     truthFilePaths = getTruthFilePaths(truthSwathDir, truthFileType, logFileID=False)
              
#     # Get truth file header info
#     if(not(useExistingTruth)):
#         truthHeaderDF = getTruthHeaders(truthFilePaths, truthFileType, logFileID=False)
#     else:
#         truthHeaderDF = False
#     # endIf
    
#     # Call getAtlTruthSwath
#     print('RUNNING getAtlTruthSwath...\n')
#     atlTruthData = getAtlTruthSwath(atl03legacy, rotationData, 
#                                     truthHeaderDF, truthFilePaths,
#                                     buffer, outFilePath, createTruthFile, 
#                                     truthFileType, useExistingTruth, 
#                                     logFileID=False)

#     # Reclassify anything 
#     atlTruthData.classification[atlTruthData.classification == 3] = 4
#     atlTruthData.classification[atlTruthData.classification == 5] = 4
    
#     return atlTruthData
    
# def legacy_get_meas_error(atl03legacy, atlTruthData, rotationData, outFilePath, truthgroundclass=2):
#     # truthgroundclass = 2
#     offsetsCrossTrackBounds = np.array([-48, 48])      # Cross-track limits
#     offsetsAlongTrackBounds = np.array([-48, 48 ])      # Along-track limits 
#     offsetsRasterResolutions = np.array([8, 4, 2, 1])  # Step-down resolutions
#     offsetsUseVerticalShift = False    # Option to use a vertical shift
#     offsetsVerticalShift = 0   # Vertical shift to use if above set to True
#     #measClassFilter = 1 # Meas Classes (0 = Unclass, 1 = Ground)
#     filterData = truthgroundclass 
#     useMeasSigConf = False # Use measured signal confidence
#     offsets = offsetsStruct(offsetsCrossTrackBounds, offsetsAlongTrackBounds, 
#                             offsetsRasterResolutions, offsetsUseVerticalShift, 
#                             offsetsVerticalShift )
#     createMeasCorrFile = True # Option to create ouput measured corrected .las
#     makePlots = False         # Option to make output plots
#     showPlots = False         # Option to show output plot windows
#     refHeightType = 'HAE'

#     atlCorrections = getMeasurementError(atl03legacy, atlTruthData, 
#                                          refHeightType, 
#                                          rotationData, outFilePath, 
#                                          useMeasSigConf, filterData, offsets, 
#                                          createMeasCorrFile, makePlots, 
#                                          showPlots)

#     return atlCorrections
    
# def apply_offset_legacy(atl03struct, atlCorrections):
#     atl03struct.df.alongtrack = atl03struct.df.alongtrack +\
#         atlCorrections.alongTrackBounds
#     atl03struct.df.crosstrack = atl03struct.df.crosstrack +\
#         atlCorrections.crossTrackBounds
#     atl03struct.df.h_ph = atl03struct.df.h_ph + atlCorrections.verticalShift

def loadTruthFile(truthFilePath, atlMeasuredData, rotationData, truthFileType, outFilePath, logFileID=False):
    
    # Initialize output
    atlTruthData = []
    
    # Determine which file type to load
    if(('las' in truthFileType.lower()) or ('laz' in truthFileType.lower())):  
        # Get .las header info
        atlTruthData = loadLasFile(truthFilePath, atlMeasuredData, rotationData, logFileID)             
    elif('tif' in truthFileType):
        # Load .tif file
        atlTruthData = loadTifFile(truthFilePath, atlMeasuredData, rotationData, outFilePath, logFileID)
    # endIf
    
    return atlTruthData

def las_to_df(lasFilePath):
    
    # if metadata is None:
        # Read contents of .las file
    with File(lasFilePath, mode = 'r') as lasFile:
        
        # Create Pandas DF
        las_df = pd.DataFrame()
        # Store output from .las file
        las_df['x'] = lasFile.x
        las_df['y'] = lasFile.y
        las_df['z'] = lasFile.z
        las_df['classification'] = lasFile.classification
        las_df['intensity'] = lasFile.intensity
        las_df['date'] = lasFile.header.get_date()
        
        # Store output into class structure
        # lasData = lasStruct(x, y, z, classification, intensity, headerData)
        
    # EndWith
    
    return las_df

def loadLasFile(truthFilePath, epsg_atl, rotationData, decimate_n = 3):
    
    # Read .las file
    lasTruthData = las_to_df(truthFilePath)
    lasTruthData = lasTruthData.iloc[::decimate_n, :]
    
    # Find EPSG Code from truth file
    truthHeader = readLasHeader(truthFilePath)
    epsg_truth = truthHeader['epsg'][0]
    
    # Find EPSG Code from input file
    # epsg_atl = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
    
        
    # Reproject if necessary
    if(epsg_truth == 'None'):
        
        print('      *WARNING: Invalid reference EPSG code, skipping file.')
        atlTruthData = False
        
    else:
        
        if(epsg_truth != epsg_atl):
            
            # If EPSG code does not match, reproject to input EPSG code
            print('      *Reference file EPSG code does not match ICESat-2, reprojecting reference file...')
            lasTruthData['easting'], lasTruthData['northing'] = transform(epsg_truth, epsg_atl, np.array(lasTruthData.x), np.array(lasTruthData.y))
            
        # endIf
        
        # Rotate TRUTH data to CT/AT plane
        lasTruthData['crosstrack'], lasTruthData['alongtrack'], _, _, _, _ = \
        getCoordRotFwd(np.array(lasTruthData.easting), np.array(lasTruthData.northing), 
                       rotationData.R_mat, rotationData.xRotPt, 
                       rotationData.yRotPt, rotationData.desiredAngle)
        
        # Get reference lat/lon
        # lasTruth_lat, lasTruth_lon = getUTM2LatLon(lasTruthData_x, lasTruthData_y,
        #                                            atlMeasuredData.zone, atlMeasuredData.hemi)
        lasTruthData['lon'], lasTruthData['lat'] = transform(epsg_truth, 'epsg:4326', 
                                                   np.array(lasTruthData.x), 
                                                   np.array(lasTruthData.y))

        # Store data as object
        # atlTruthData = atlTruthStruct(lasTruthData_x, lasTruthData_y, 
        #                               lasTruthData_x_newRot, lasTruthData_y_newRot, 
        #                               lasTruth_lon, lasTruth_lat,
        #                               lasTruthData.z, 
        #                               lasTruthData.classification, 
        #                               lasTruthData.intensity,
        #                               lasTruthData_year,
        #                               lasTruthData_month,
        #                               lasTruthData_day,                                      
        #                               atlMeasuredData.zone, 
        #                               atlMeasuredData.hemi,
        #                               epsg_atl)
    # endIf
    
    return lasTruthData

### Function to load .tif file
def loadTifFile(truthFilePath, epsg_atl, rotationData):
    
    # Read Tif file
    xarr0, yarr0, zarr, intensity, classification, epsg = formatDEM(truthFilePath)
    
    # Convert ints to floats
    xarr0 = xarr0.astype(float)
    yarr0 = yarr0.astype(float)
    zarr = zarr.astype(float)
    
    # Find EPSG Code from tif
    epsg_truth = 'epsg:' + epsg
    
    # Determine if EPSG Code is the same for the ATL03 Measured
    # epsg_atl = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
    
    # Store values
    xarr = xarr0
    yarr = yarr0
    
    # Reproject if necessary
    if(epsg_truth == 'None'):
        
        print('      *WARNING: Invalid reference EPSG code, skipping file.')
        # atlTruthData = False
        tif_df = pd.DataFrame()
        
    else:
        
        if(epsg_truth != epsg_atl):
            
            # If EPSG code does not match, use GDAL Warp to create new Tif
            print('      *Reference file EPSG code does not match ICESat-2, reprojecting reference file...')
            xarr, yarr = transform(epsg_truth, epsg_atl, xarr0, yarr0)
                    
        # Rotate Data for along-track/cross-track
        x_newRot, y_newRot, _, _, _, _ = getCoordRotFwd(xarr, yarr, 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
        
        # Get reference lat/lon
        lat, lon = transform(epsg_truth, 'epsg:4326', xarr, yarr)
        
        # Store Data as Object
        tif_df = pd.DataFrame()
        tif_df['x'] = xarr
        tif_df['y'] = xarr
        tif_df['z'] = xarr
        tif_df['lon'] = xarr
        tif_df['lat'] = xarr
        tif_df['alongtrack'] = xarr
        tif_df['crosstrack'] = xarr
        tif_df['classification'] = xarr
        tif_df['date'] = datetime.datetime.now() # TODO: find actual date of tif

    return tif_df

### Function to find and reproject min/max extents in header data
def reprojectHeaderData(truthHeaderDF, epsg_atl):
    
    # Copy dataframe
    truthHeaderNewDF = truthHeaderDF.copy()
    
    # Input EPSG code
    # epsg_atl = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
    
    # Loop through all truth file EPSG codes
    for i in range(0,len(truthHeaderDF)):
#        writeLog(str(i) + ' out of ' + str(len(truthHeaderDF)), logFileID)
        # Reproject EPSG code to match input EPSG if unmatching
        epsg_truth = truthHeaderDF['epsg'][i]
#        writeLog(epsg_truth, logFileID)
#        writeLog(epsg_atl, logFileID)
        if(epsg_truth!=epsg_atl):
            
            # Header extents
            xmin = truthHeaderDF['xmin'][i]
            xmax = truthHeaderDF['xmax'][i]
            ymin = truthHeaderDF['ymin'][i]
            ymax = truthHeaderDF['ymax'][i]
            x = [xmin, xmax]
            y = [ymin, ymax]
            
            # Reproject extents to input EPSG code
            try:

                xout, yout = transform(epsg_truth, epsg_atl, x, y)
                
                # Store new extents into output dataframe
                truthHeaderNewDF['xmin'][i] = xout[0]
                truthHeaderNewDF['xmax'][i] = xout[1]
                truthHeaderNewDF['ymin'][i] = yout[0]
                truthHeaderNewDF['ymax'][i] = yout[1]
                truthHeaderNewDF['epsg'][i] = epsg_atl
            
            except:
                
                # Store new extents into output dataframe
                print('WARNING: Cannot reproject data, skipping file: %s' %truthHeaderNewDF['fileName'][i])
                truthHeaderNewDF['xmin'][i] = 'None'
                truthHeaderNewDF['xmax'][i] = 'None'
                truthHeaderNewDF['ymin'][i] = 'None'
                truthHeaderNewDF['ymax'][i] = 'None'
                truthHeaderNewDF['epsg'][i] = 'None'
            
            # endTry
                
            
        # endIf
    # endFor
    
    return truthHeaderNewDF

### Function to find which truth tiles ICESat-2 crosses over
def findMatchingTruthFiles(truthHeaderNewDF, atlMeasuredData, rotationData, buffer):
       
    # Get MEASURED rotated buffer bounds data
    xRotL = atlMeasuredData.crosstrack - buffer
    xRotR = atlMeasuredData.crosstrack + buffer
    yRot  = atlMeasuredData.alongtrack
    
    # Get MEASURED buffer bounds data in easting/northing plane
    xL, yL,  _, _, _ = getCoordRotRev(xRotL, yRot, rotationData.R_mat, rotationData.xRotPt, rotationData.yRotPt)
    xR, yR,  _, _, _ = getCoordRotRev(xRotR, yRot, rotationData.R_mat, rotationData.xRotPt, rotationData.yRotPt)
    
    # Rotate truth header file min/max x,y points to CT/AT plane
    xMinyMinHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(((truthHeaderNewDF['xmin']).to_numpy()).astype('float'), 
                                                        ((truthHeaderNewDF['ymin']).to_numpy()).astype('float'), 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
    
    xMinyMaxHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(((truthHeaderNewDF['xmin']).to_numpy()).astype('float'), 
                                                        ((truthHeaderNewDF['ymax']).to_numpy()).astype('float'), 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
    xMaxyMinHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(((truthHeaderNewDF['xmax']).to_numpy()).astype('float'), 
                                                        ((truthHeaderNewDF['ymin']).to_numpy()).astype('float'), 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
    xMaxyMaxHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(((truthHeaderNewDF['xmax']).to_numpy()).astype('float'), 
                                                        ((truthHeaderNewDF['ymax']).to_numpy()).astype('float'), 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
    
    # Find min/max x/y header points inside min/max x buffer points
    xMinyMinXPtsInBuffer = (xMinyMinHeaderRotX >= xRotL.min()) & \
    (xMinyMinHeaderRotX <= xRotR.max())
    xMinyMaxXPtsInBuffer = (xMinyMaxHeaderRotX >= xRotL.min()) & \
    (xMinyMaxHeaderRotX <= xRotR.max())
    xMaxyMinXPtsInBuffer = (xMaxyMinHeaderRotX >= xRotL.min()) & \
    (xMaxyMinHeaderRotX <= xRotR.max())
    xMaxyMaxXPtsInBuffer = (xMaxyMaxHeaderRotX >= xRotL.min()) & \
    (xMaxyMaxHeaderRotX <= xRotR.max())
    
    # Get header points inside MEASURED buffer points
    xHeaderPtsInBuffer = np.c_[xMinyMinXPtsInBuffer | xMinyMaxXPtsInBuffer| xMaxyMinXPtsInBuffer | xMaxyMaxXPtsInBuffer]
    
    # Find min/max x buffer points inside min/max x/y header points
    xyMinMaxHeaderRot = np.column_stack((xMinyMinHeaderRotX,xMinyMaxHeaderRotX,xMaxyMinHeaderRotX,xMaxyMaxHeaderRotX))
    xyMinHeaderRot = np.c_[np.amin(xyMinMaxHeaderRot,axis=1)]
    xyMaxHeaderRot = np.c_[np.amax(xyMinMaxHeaderRot,axis=1)]
    xMinBufferPtsInFile = (xRotL.min() >= xyMinHeaderRot) & (xRotL.min() <= xyMaxHeaderRot)
    xMaxBufferPtsInFile = (xRotR.max() >= xyMinHeaderRot) & (xRotR.max() <= xyMaxHeaderRot)
    
    # Get MEASURED buffer points inside header points
    xBufferPtsInHeader = xMinBufferPtsInFile | xMaxBufferPtsInFile
    
    # Get any points where buffer points are inside header and vice versa
    xPtsInFile = np.ravel(np.logical_or(xHeaderPtsInBuffer, xBufferPtsInHeader))
        
    # Get matching truth file names and x/y min/max points
    matchingFilesPre = ((truthHeaderNewDF['fileName'][xPtsInFile]).to_numpy()).astype('str')
    matchingHeaderXmin = ((truthHeaderNewDF['xmin'][xPtsInFile]).to_numpy()).astype('float')
    matchingHeaderXmax = ((truthHeaderNewDF['xmax'][xPtsInFile]).to_numpy()).astype('float')
    matchingHeaderYmin = ((truthHeaderNewDF['ymin'][xPtsInFile]).to_numpy()).astype('float')
    matchingHeaderYmax = ((truthHeaderNewDF['ymax'][xPtsInFile]).to_numpy()).astype('float')
    
    fileNumsAll = np.arange(0,len(truthHeaderNewDF))
    allFileInds = fileNumsAll[xPtsInFile]
    
    # Get all MEASURED x,y buffer points
    xAll = np.concatenate((xL, xR))
    yAll = np.concatenate((yL, yR))
    
    matchTF = []
    if(len(matchingFilesPre)>0):
        for i in range(0,len(matchingFilesPre)):
            
            # Determine TRUTH files where MEASURED data actually crosses over
            xPtsInFile = (xAll >= matchingHeaderXmin[i]) & (xAll <= matchingHeaderXmax[i])
            yPtsInFile = (yAll >= matchingHeaderYmin[i]) & (yAll <= matchingHeaderYmax[i])
            anyPtsInFile = any(xPtsInFile & yPtsInFile)
    
            # If a TRUTH file is a match, use it
            if(anyPtsInFile):
                matchTF.append(True)
            else:
                matchTF.append(False)
            # endIf
            
        # endFor
    # endIf
    
    if(len(matchTF)>0):
        matchTF = np.array(matchTF)
        matchingFiles = matchingFilesPre[matchTF]
    else:
        matchingFiles = []
    # endIf
    
    matchingFileInds = allFileInds[matchTF]
    
    return matchingFiles, matchingFileInds

def make_buffer(atl03, truth_df, buffer):

    # Index Matching
    indices = indexMatch(atl03.df.alongtrack, truth_df.alongtrack)
    
    # Filter indices
    indices[indices >= len(atl03.df.crosstrack)] = (len(atl03.df.crosstrack) - 1)
    x_check = atl03.df.crosstrack[indices]
    x_diff = np.array(truth_df.crosstrack) - np.array(x_check)
    filter_data = np.where((x_diff < buffer) & (x_diff > -buffer))[0]
    
    # TODO: iloc is slow, this could be sped up by separating into PD array
    return truth_df.iloc[filter_data]

# Unit test
if __name__ == "__main__":
    import os
    # from getAtlTruthSwath_auto import getAtlTruthSwath
    # from getMeasurementError_auto import getMeasurementError, offsetsStruct
    from icesatReader import get_atl03_struct
    from icesatReader import convert_atl03_to_legacy
    # from icesatReader import get_atl_alongtrack
    import pandas as pd

    out_folder = 'E:/0_data/is2/prf/' 
    if os.name == 'nt':
        # basepath03 = 'Z:/data/release/002/ATL03_r002/Finland/'
        # basepath08 = 'Z:/data/release/002/ATL08_r002/Finland/'
        basepath03 = 'E:/0_data/is2/prf/ATL03/'
        basepath08 = 'E:/0_data/is2/prf/ATL08/'
    else:
        basepath03 = '/laserpewpew/data/release/002/ATL03_r002/Finland/'
        basepath08 = '/laserpewpew/data/release/002/ATL08_r002/Finland/'

    # atl03file = 'ATL03_20181118120428_07770103_002_01.h5'
    # atl08file = 'ATL08_20181118120428_07770103_002_01.h5'

    atl03_list = os.listdir(basepath03)
    atl08_list = os.listdir(basepath08)
    atl03file = atl03_list[2]
    atl08file = atl08_list[2]

    # Inputs
    atl03filepath = basepath03 + atl03file
    atl08filepath = basepath08 + atl08file
    gt = 'gt1r'
    
    header_file_path =\
        '/LIDAR/server/USERS/eric/1_experiment/Finland_HeaderData.mat'
        
    kml_bounds_txt1 = '/LIDAR/server/USERS/eric/2_production/kmlBounds.txt'
   
    print('Generate ATL03 Struct')
    atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath, 
                             epsg = '32618', kml_bounds_txt = kml_bounds_txt1, 
                             header_file_path = header_file_path)    
    
    # atl03.df = atl03.df[atl03.df['time'] < 12]

    ###### Truth Header Section #######
    # If existing truth header, read it, if there is none, create it
    truthSwathDir = 'E:/data/2018spl_2959_6647'
    truthFileType = 'laz'
    truthFilePaths = getTruthFilePaths(truthSwathDir, truthFileType, logFileID=False)
    
    truthHeaderDF = getTruthHeaders(truthFilePaths, truthFileType, logFileID=False)

    ##### Reproject Truth Header ######
    epsg_atl = 'epsg:32618'
    truthHeaderNewDF = reprojectHeaderData(truthHeaderDF, epsg_atl)
    
    # Filter atl03.df by size of reproject
    atl03.df = atl03.df[atl03.df.northing > np.min(truthHeaderNewDF.ymin)]
    atl03.df = atl03.df[atl03.df.northing < np.max(truthHeaderNewDF.ymax)]
    atl03.df = atl03.df.reset_index()
    atl03.df, atl03.rotationData = get_atl_alongtrack(atl03.df)
        
    # Find truth files that intersect ICESat-2 track
    #### atlMeasuredData = atlMeasuredData.df
    buffer = 25
    _, matchingTruthFileInds = findMatchingTruthFiles(truthHeaderNewDF, 
                                                       atl03.df,  atl03.rotationData, buffer)
    matchingTruthFiles = np.array(truthFilePaths)[matchingTruthFileInds]
    
    print(matchingTruthFiles)
    
    # truthFilePath = matchingTruthFiles[1]
    
    truth_swath = pd.DataFrame()
    
    for i in range(0,len(matchingTruthFiles)):
        truth_df = loadLasFile(matchingTruthFiles[i], epsg_atl, atl03.rotationData, decimate_n = 5)
        truth_df['classification'] = ace(np.array(truth_df.easting), np.array(truth_df.northing), 
                     np.array(truth_df.z), np.array(truth_df.classification))
        truth_df = make_buffer(atl03, truth_df, buffer)
        truth_swath = truth_swath.append(truth_df)

    # Calculate off-set corrections for ATL03        
    atlCorrections = getMeasurementError(atl03, truth_swath)

    # Apply along-track/cross-track, and height corrections to ATL03    
    # atl03.df.alongtrack = atl03.df.alongtrack + atlCorrections.alongTrack 
    # atl03.df.crosstrack = atl03.df.crosstrack + atlCorrections.crossTrack 
    # atl03.df.h_ph = atl03.df.h_ph + atlCorrections.z 
    
    truth_swath.alongtrack = truth_swath.alongtrack - atlCorrections.alongTrack 
    truth_swath.crosstrack = truth_swath.crosstrack - atlCorrections.crossTrack 
    truth_swath.z = truth_swath.z - atlCorrections.z 
    
    # # Apply beam width measurement
    truth_swath = make_buffer(atl03, truth_swath, 4)
    # truth_swath.z  = truth_swath.z  +
    # Do Perfect Classifier
    measpc, measoc = perfect_classifier(atl03, truth_swath,ground = [2],canopy = [3,4,5], 
                      unclassed = [1, 6, 7, 18], keepsize = True)
    
    atl03.df['truth_label'] = measpc
    
    # Do Binner
    from icesatBinner import calc_truth_binning_radiometry
    from icesatBinner import calc_atl03_binning_radiometry
    
    # calc_truth_binning_radiometry
    
    # Graphing
        
    # import matplotlib.pyplot as plt

    # plt.figure()
    # ax1 = plt.subplot(211)
    # plt.plot(truth_swath.alongtrack[truth_swath.classification > 2],truth_swath.z[truth_swath.classification > 2],'.',color=[0.6,0.6,0.6])    
    # plt.plot(truth_swath.alongtrack[truth_swath.classification == 2],truth_swath.z[truth_swath.classification == 2],'.',color=[0.3,0.3,0.3])

    # plt.plot(atl03.df.alongtrack[atl03.df.truth_label == 0],atl03.df.h_ph[atl03.df.truth_label == 0],'.',color=[0.0,0.95,0.95])
    # plt.plot(atl03.df.alongtrack[atl03.df.classification == 3],atl03.df.h_ph[atl03.df.classification == 3],'.',color=[0,0.9,0])    
    # plt.plot(atl03.df.alongtrack[atl03.df.classification == 2],atl03.df.h_ph[atl03.df.classification == 2],'.',color=[0.05,0.5,0.12])
    # plt.plot(atl03.df.alongtrack[atl03.df.classification == 1],atl03.df.h_ph[atl03.df.classification == 1],'.',color=[0.95,0.5,0.0])

    # plt.subplot(212, sharex=ax1, sharey=ax1)
    # plt.plot(truth_swath.alongtrack[truth_swath.classification > 2],truth_swath.z[truth_swath.classification > 2],'.',color=[0.6,0.6,0.6])    
    # plt.plot(truth_swath.alongtrack[truth_swath.classification == 2],truth_swath.z[truth_swath.classification == 2],'.',color=[0.3,0.3,0.3])

    # plt.plot(atl03.df.alongtrack[atl03.df.truth_label == 0],atl03.df.h_ph[atl03.df.truth_label == 0],'.',color=[0.0,0.95,0.95])
    # plt.plot(atl03.df.alongtrack[atl03.df.truth_label == 3],atl03.df.h_ph[atl03.df.truth_label == 3],'.',color=[0,0.9,0])    
    # plt.plot(atl03.df.alongtrack[atl03.df.truth_label == 2],atl03.df.h_ph[atl03.df.truth_label == 2],'.',color=[0.05,0.5,0.12])
    # plt.plot(atl03.df.alongtrack[atl03.df.truth_label == 1],atl03.df.h_ph[atl03.df.truth_label == 1],'.',color=[0.95,0.5,0.0])

    
    # plt.plot(test.alongtrack[truth_swath.classification > 2],truth_swath.norm_h[test.classification > 2],'.',color=[0.6,0.6,0.6])    
    # plt.plot(test.alongtrack[truth_swath.classification == 2],truth_swath.norm_h[test.classification == 2],'.',color=[0.3,0.3,0.3])
