# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 11:18:58 2019

@author: malonzo
"""

# Import modules
import os
import numpy as np
import time as runTime
import ntpath
from osgeo import gdal, ogr
from getAtlMeasuredSwath_auto import getAtlMeasuredSwath
from icesatIO import (readLas, writeLas, readGeoidFile, readDEMepsg, formatDEM)
from icesatUtils import (getCoordRotFwd, getCoordRotRev, getLatLon2UTM,
                         getGeoidHeight, identifyEPSG, transform)


class atlTruthStruct:
        
    # Define class with designated fields
    def __init__(self, easting, northing, crossTrack, alongTrack, z, 
                 classification, intensity, zone, hemi):
            
        self.easting = np.c_[easting]
        self.northing = np.c_[northing]
        self.crossTrack = np.c_[crossTrack]
        self.alongTrack = np.c_[alongTrack]
        self.z = np.c_[z]
        self.classification = np.c_[classification]
        self.intensity = np.c_[intensity]
        self.zone = zone
        self.hemi = hemi
        
        
def getAtlTruthSwath(atlMeasuredData, headerData, rotationData, 
                     useExistingTruth, truthSwathDir, buffer, 
                     outFilePath, createTruthFile):
    
    # Start timer
    timeStart = runTime.time()
    
    # Print message
    print('   Ground Track Number: %s' % atlMeasuredData.gtNum)
    
    # Initialze truthFilePath/truthDir
    truthFilePath = []
    truthDir = []
    # atlTruthDataFiltered = []
    
    # Initialize dem/lidar flags
    demfilecheck = False
    lidarfilecheck = False
    
    # Determine if truthSwathDir is folder or file
    if os.path.isfile(truthSwathDir) == True:
        file = ntpath.basename(truthSwathDir) 
        truthFilePath = truthSwathDir
        extension = file.split('.')[-1]
        if (extension == 'las') or (extension == 'laz'):
            print('File is lidar LAS/LAZ File')
            lidarfilecheck = True
        else:
            print('File is DSM/DTM .tif File')
            demfilecheck = True
        # endIf
    else:
        truthDir = truthSwathDir
    # endIf
    
    # Get rotation data
    R_mat = rotationData.R_mat
    xRotPt = rotationData.xRotPt
    yRotPt = rotationData.yRotPt
    desiredAngle = rotationData.desiredAngle
    phi = rotationData.phi
    
    # If truth swath exists, read .las file
    if(useExistingTruth and truthFilePath and demfilecheck == False):
        
        # Read in truth file .las file
        print('   Reading in Lidar Truth File: %s' % truthFilePath)
        lasTruthData = readLas(truthFilePath)

        # Rotate TRUTH data to CT/AT plane
        print('   Rotating Data to CT/AT Frame...')
        lasTruthData_x_newRot, lasTruthData_y_newRot, _, _, _, _ = \
        getCoordRotFwd(lasTruthData.x, lasTruthData.y, 
                       R_mat, xRotPt, yRotPt, desiredAngle)

#        # Store data as object  
#        atlTruthData = atlTruthStruct(lasTruthData.x, lasTruthData.y, 
#                                      lasTruthData_x_newRot, 
#                                      lasTruthData_y_newRot, 
#                                      lasTruthData.z, 
#                                      lasTruthData.classification, 
#                                      lasTruthData.intensity, 
#                                      headerData.zone, headerData.hemi)
        
        # Store data as object (CHANGE AFTER DEMO!!!)
        atlTruthData = atlTruthStruct(lasTruthData.x, lasTruthData.y, 
                                      lasTruthData_x_newRot, 
                                      lasTruthData_y_newRot, 
                                      lasTruthData.z, 
                                      lasTruthData.classification, 
                                      lasTruthData.intensity, 
                                      atlMeasuredData.zone, 
                                      atlMeasuredData.hemi)
    
    elif(demfilecheck == True):
        ###New
        
        # 1. Find EPSG Code from DEM
        epsg_dem = readDEMepsg(truthFilePath)
        
        # 2. Determine if EPSG Code is the same for the ATL03 Measured
        epsg_atl = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
        epsg_atl = epsg_atl.split(':')[-1]
        
        if epsg_atl == epsg_dem:
            print('   DEM Projection matches target')
            # 3b. If EPSG code is the same, set DEM as a target
            truthFilePathDEM = truthFilePath
        else:
            print('   Reprojecting DEM')
            # 3a. If EPSG code does not match, use GDAL Warp to create new DEM
            file = ntpath.basename(truthSwathDir) 
            truthFilePath = truthSwathDir
            originalname = file.split('.')[0]
            newname = originalname + "_projected.tif"
            truthFilePathDEM = os.path.normpath(os.path.join(outFilePath,newname))
            srs = ogr.osr.SpatialReference()
            srs.ImportFromEPSG(np.int(epsg_atl))
            warpoptions = gdal.WarpOptions(dstSRS = srs)
            gdal.Warp(truthFilePathDEM, truthFilePath, options = warpoptions)

        # 4. Call the Format DEM Function
        print('   Reading in DEM Truth File: %s' % truthFilePathDEM)
        xarr, yarr, zarr, intensity, classification, epsg = \
        formatDEM(truthFilePathDEM)           
        
        # 5. Run a quick filter to remove points not even remotely close 
        print('   Running Quick Filter')
        xmin = np.min(atlMeasuredData.easting.flatten()) - 100
        xmax = np.max(atlMeasuredData.easting.flatten()) + 100
        ymin = np.min(atlMeasuredData.northing.flatten()) - 100
        ymax = np.max(atlMeasuredData.northing.flatten()) + 100
        
        quickfilter = np.where((xarr > xmin) & (xarr < xmax) & (yarr > ymin) & (yarr < ymax))
        
        xcoord = xarr[quickfilter]
        ycoord = yarr[quickfilter]
        zarr = zarr[quickfilter]
        intensity = intensity[quickfilter]
        classification = classification[quickfilter]
        
        
        
#        ###Original
#
#        
#        # Apply EPSG code to X,Y
#        # 'epsg:32613'
#        epsg_in = 'epsg:' + str(epsg)
#        epsg_out = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
#        
#        if epsg_in == epsg_out:
#            print('   DEM and Measured Match EPSG Code')
#            xcoord = xarr
#            ycoord = yarr
#        else:
#            print('   Transforming Projection for DEM')
#            xx,yy = zip(*[transform(epsg_in, epsg_out, xarr[i],yarr[i])  \
#                          for i in range(0,len(xarr))])
#            xcoord = np.asarray(xx)
#            ycoord = np.asarray(yy)
        
        
        # Rotate Data for along-track/cross-track
        print('   Rotating Data to CT/AT Frame...')
        x_newRot, y_newRot, _, _, _, _ = getCoordRotFwd(xcoord, ycoord, R_mat, 
                                                        xRotPt, yRotPt, 
                                                        desiredAngle)
        
        # Store Data as Object
        atlTruthData = atlTruthStruct(xcoord, ycoord, x_newRot, y_newRot, 
                                      zarr, classification, intensity, 
                                      atlMeasuredData.zone, 
                                      atlMeasuredData.hemi)
        print('   DEM being used for ATL Truth Data')

    else:
        
        # Create new truth .las file
        print('   Creating new Truth File...')
        
        # Get MEASURED rotated buffer bounds data
        xRotL = atlMeasuredData.crossTrack - buffer
        xRotR = atlMeasuredData.crossTrack + buffer
        yRot = atlMeasuredData.alongTrack
        
        # Get MEASURED buffer bounds data in easting/northing plane
        xL, yL,  _, _, _ = getCoordRotRev(xRotL, yRot, R_mat, xRotPt, yRotPt)
        xR, yR,  _, _, _ = getCoordRotRev(xRotR, yRot, R_mat, xRotPt, yRotPt)
                    
        # Rotate truth header file min/max x,y points to CT/AT plane
        xMinyMinHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(headerData.xmin, 
                                                            headerData.ymin, 
                                                            R_mat, 
                                                            xRotPt, yRotPt, 
                                                            desiredAngle)
        xMinyMaxHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(headerData.xmin, 
                                                            headerData.ymax, 
                                                            R_mat, xRotPt, 
                                                            yRotPt, 
                                                            desiredAngle)
        xMaxyMinHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(headerData.xmax, 
                                                            headerData.ymin, 
                                                            R_mat, 
                                                            xRotPt, yRotPt, 
                                                            desiredAngle)
        xMaxyMaxHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(headerData.xmax, 
                                                            headerData.ymax, 
                                                            R_mat, 
                                                            xRotPt, yRotPt, 
                                                            desiredAngle)
        
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
        xHeaderPtsInBuffer = xMinyMinXPtsInBuffer | xMinyMaxXPtsInBuffer| xMaxyMinXPtsInBuffer | xMaxyMaxXPtsInBuffer
        
        # Find min/max x buffer points inside min/max x/y header points
        xyMinMaxHeaderRot = np.concatenate((xMinyMinHeaderRotX,xMinyMaxHeaderRotX,xMaxyMinHeaderRotX,xMaxyMaxHeaderRotX),axis=1)
        xyMinHeaderRot = np.c_[np.amin(xyMinMaxHeaderRot,axis=1)]
        xyMaxHeaderRot = np.c_[np.amax(xyMinMaxHeaderRot,axis=1)]
        xMinBufferPtsInFile = (xRotL.min() >= xyMinHeaderRot) & (xRotL.min() <= xyMaxHeaderRot)
        xMaxBufferPtsInFile = (xRotR.max() >= xyMinHeaderRot) & (xRotR.max() <= xyMaxHeaderRot)
        
        # Get MEASURED buffer points inside header points
        xBufferPtsInHeader = xMinBufferPtsInFile | xMaxBufferPtsInFile
        
        # Get any points where buffer points are inside header and vice versa
        xPtsInFile = xHeaderPtsInBuffer | xBufferPtsInHeader
        
        # Get matching truth file names and x/y min/max points
        matchingFiles = headerData.tileName[xPtsInFile]
        matchingHeaderXmin = headerData.xmin[xPtsInFile]
        matchingHeaderXmax = headerData.xmax[xPtsInFile]
        matchingHeaderYmin = headerData.ymin[xPtsInFile]
        matchingHeaderYmax = headerData.ymax[xPtsInFile]
        
        # Get all MEASURED x,y buffer points
        xAll = np.concatenate((xL, xR))
        yAll = np.concatenate((yL, yR))
        
        # Initialize loop variables
        counter = 0
        subData_x_old = []
        subData_y_old = []
        subData_z_old = []
        subData_classification_old = []
        subData_intensity_old = []
        
        # Find matching truth files
        print('   Looking in Truth File Path: %s' % truthDir)
        print('   Matching Truth .las Files:')
        for i in range(0,len(matchingFiles)):
            
            # Determine TRUTH files where MEASURED data actually crosses over
            xPtsInFile = (xAll >= matchingHeaderXmin[i]) & (xAll <= matchingHeaderXmax[i])
            yPtsInFile = (yAll >= matchingHeaderYmin[i]) & (yAll <= matchingHeaderYmax[i])
            anyPtsInFile = any(xPtsInFile & yPtsInFile)

            # If a TRUTH file is a match, use it
            if(anyPtsInFile):
                
                # Print TRUTH file name to screen
                fileNum = counter + 1
                print('   %s) %s' % (fileNum, matchingFiles[i]) )
                
                # Build path to TRUTH file 
                truthLasFile = os.path.normpath(truthDir + '/' + matchingFiles[i]) 
    
                # Try to read .las file
                try:
                    
                    # Read TRUTH .las/.laz file
                    # truthLasFile = truthLasFile.split('.')[0][:-3] + 'las.las'
                    # print(truthLasFile)
                    lasData = readLas(truthLasFile)
                    
                except:
                    
                    # Print Error message
                    print('    ERROR: Could not read file.')
                    lasData = False
                    
                # endTry
                
                # If there is .las data, process it
                if(lasData):
                    
                    # Covert to UTM coords if in lat/lon
                    if(headerData.coordType == 'Lat/Lon'):
                        lasDataNewX, lasDataNewY, _, _ = getLatLon2UTM(lasData.x,lasData.y)
                        lasData = atlTruthStruct(lasDataNewX, lasDataNewY, lasData.z, lasData.classification, lasData.intensity, headerData.zone, headerData.hemi)
                    # endif
                    
                    # Rotate TRUTH lasData to CT/AT plane
                    xRotLasData, yRotLasData,  _, _, _, _ = getCoordRotFwd(lasData.x, lasData.y, R_mat, xRotPt, yRotPt, desiredAngle)
    
                    # Find TRUTH lasData points inside TRUTH buffer
                    yRotLocalInds = (yRot >= yRotLasData.min()) & (yRot <= yRotLasData.max()) 
                    xRotLocalL = xRotL[yRotLocalInds]
                    xRotLocalR = xRotR[yRotLocalInds]
                    xyBufferInds = (xRotLasData >= xRotLocalL.min()) & (xRotLasData <= xRotLocalR.max())
                    
                    # Extract TRUTH lasData points inside TRUTH buffer
                    subData_x = lasData.x[xyBufferInds]
                    subData_y = lasData.y[xyBufferInds]
                    subData_z = lasData.z[xyBufferInds]
                    subData_classification = lasData.classification[xyBufferInds]
                    subData_intensity = lasData.intensity[xyBufferInds]
                    
                    # Append TRUTH buffer points onto previous data
                    subData_x_new = np.append(subData_x_old,subData_x)
                    subData_y_new = np.append(subData_y_old,subData_y)
                    subData_z_new = np.append(subData_z_old,subData_z)
                    subData_classification_new = np.append(subData_classification_old,subData_classification)
                    subData_intensity_new = np.append(subData_intensity_old,subData_intensity)
                    
                    # Store old TRUTH buffer data
                    subData_x_old = subData_x_new
                    subData_y_old = subData_y_new
                    subData_z_old = subData_z_new
                    subData_classification_old = subData_classification_new
                    subData_intensity_old = subData_intensity_new
                
                # endIf
                
                # Increment counter
                counter += 1
                
            # endIf
        # endFor
        
        # Store TRUTH swath data if it exists
        if('subData_x_new' in locals()):
        
            # Rotate TRUTH data to CT/AT plane
            print('   Rotating Data to CT/AT Frame...')
            subData_x_newRot, subData_y_newRot, _, _, _, _ = getCoordRotFwd(subData_x_new, subData_y_new, R_mat, xRotPt, yRotPt, desiredAngle)
    
            # Store data in class structure
            atlTruthData = atlTruthStruct(subData_x_new, subData_y_new, subData_x_newRot, subData_y_newRot, subData_z_new, subData_classification_new, subData_intensity_new, headerData.zone, headerData.hemi)
        
            # Apply geoid model Z corrections if necessary
            regionName = atlMeasuredData.kmlRegionName
            
            # Initialize geoidFile variable
            geoidFile = False
            
            # Determine geoid model to use
            if('sonoma' in regionName.lower() or 'indiana' in regionName.lower()):     
                geoidFile = 'geoid12b_latlon.mat'
            elif('finland' in regionName.lower()):
                geoidFile = 'geoidFin2005N00_latlon.mat'
            elif('rangiora' in regionName.lower() or 'wellington' in regionName.lower()):
                geoidFile = 'geoidNZ_latlon.mat'
            # endif
            
            if(geoidFile):
            
                # Print status message
                print('')
                print('   STATUS: Applying Truth Z Correction for %s (Geoid File = %s).' % (regionName, geoidFile) )
                    
                # Load Geoid file
                geoidData = readGeoidFile(geoidFile)
        
                # Get geoidal heights and add to orthometric heights
                atlTruthData = getGeoidHeight(geoidData,atlTruthData)
                
            # EndIf
        
        else:
            
            # Set to false
            print('   STATUS: No Truth Data Stored for this Track.')
            atlTruthData = False
            
        # endIf
        
    # endIf
    
    # Create output file
    if(createTruthFile and atlTruthData): 
        
        print('   Writing truth .las file...')
        outName = atlMeasuredData.atl03FileName + '_' + atlMeasuredData.gtNum + '_REFERENCE_' + str(buffer) + 'L' + str(buffer) + 'Rm_buffer.las'
        outPath = os.path.normpath(outFilePath + '/' + outName)
        
        # If output directory does not exist, create it
        if(not os.path.exists(os.path.normpath(outFilePath))):
            os.mkdir(os.path.normpath(outFilePath))
        # EndIf
        
        # Write .las file
        writeLas(np.ravel(atlTruthData.easting),np.ravel(atlTruthData.northing),np.ravel(atlTruthData.z),'utm',outPath,np.ravel(atlTruthData.classification),np.ravel(atlTruthData.intensity),None,atlTruthData.hemi,atlTruthData.zone)
    
    # endIf
    
    # End timer
    timeEnd = runTime.time()
    timeElapsedTotal = timeEnd - timeStart
    timeElapsedMin = np.floor(timeElapsedTotal / 60)
    timeElapsedSec = timeElapsedTotal % 60
    
    # Print completion message
    print('   Module Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
    print('\n')
    
    # Return object
    return atlTruthData
    
    
# Unit test
if __name__ == "__main__":
    
    print('Test')
