# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 11:18:58 2019

@author: malonzo
"""

import numpy as np
import pandas as pd
import datetime


from phoreal.utils import getCoordRotFwd, transform, getCoordRotRev, indexMatch
from phoreal.io import getTruthHeaders, readLasHeader, formatDEM, getCoordRotRev
import laspy
import pyproj
from phoreal.ace import ace
from phoreal.getMeasurementError import getMeasurementError
from phoreal.reader import get_atl_alongtrack
from phoreal.CalVal import perfect_classifier

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
    
    return atlTruthData

def las_to_df(lasFilePath):
    
    las_file = laspy.read(lasFilePath)
    las_df = pd.DataFrame()
    # Store output from .las file
    las_df['x'] = np.array(las_file.x)
    las_df['y'] = np.array(las_file.y)
    las_df['z'] = np.array(las_file.z)
    las_df['classification'] = np.array(las_file.classification)
    las_df['intensity'] = np.array(las_file.intensity)
    las_df['date'] = las_file.header.creation_date
    
    return las_df

def loadLasFile(truthFilePath, epsg_atl, rotationData, decimate_n = 3, epsg_truth = 'EPSG:0'):
    
    # Read .las file
    lasTruthData = las_to_df(truthFilePath)
    lasTruthData = lasTruthData.iloc[::decimate_n, :]
    
    # Find EPSG Code from truth file
    # truthHeader = readLasHeader(truthFilePath)
    las_file = laspy.read(truthFilePath)
    if epsg_truth == 0:
        epsg_truth = 'EPSG:' + str(las_file.header.parse_crs().to_epsg(min_confidence=1))

    # epsg_truth = truthHeader['epsg'][0]
    
    # Find EPSG Code from input file
    # epsg_atl = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
            
    # Reproject if necessary
    if(epsg_truth == 'None'):
        
        print('      *WARNING: Invalid reference EPSG code, skipping file.')
        # atlTruthData = False
        
    else:
        
        if(epsg_truth != epsg_atl):
            
            # If EPSG code does not match, reproject to input EPSG code
            print('      *Reference file EPSG code does not match ICESat-2, reprojecting reference file...')
            lasTruthData['easting'], lasTruthData['northing'] = transform(epsg_truth, epsg_atl, np.array(lasTruthData.x), np.array(lasTruthData.y))
                    
        # Rotate TRUTH data to CT/AT plane
        lasTruthData['crosstrack'], lasTruthData['alongtrack'], _, _, _, _ = \
        getCoordRotFwd(np.array(lasTruthData.easting), np.array(lasTruthData.northing), 
                       rotationData.R_mat, rotationData.xRotPt, 
                       rotationData.yRotPt, rotationData.desiredAngle)
        
        # Get reference lat/lon
        lasTruthData['lon'], lasTruthData['lat'] = transform(epsg_truth, 'epsg:4326', 
                                                   np.array(lasTruthData.x), 
                                                   np.array(lasTruthData.y))
    
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
    
    if isinstance(epsg_atl, str) and 'epsg:' not in  epsg_atl:
        epsg_atl = 'epsg:' + epsg_atl
    
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
        tif_df['y'] = yarr
        tif_df['z'] = zarr
        tif_df['lon'] = lat
        tif_df['lat'] = lon
        tif_df['alongtrack'] = x_newRot
        tif_df['crosstrack'] = y_newRot
        tif_df['classification'] = 2
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
    pass
