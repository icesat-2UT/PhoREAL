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
from icesatUtils import getRaster

import glob
import vox
import acer
from scipy import interpolate
from operator import itemgetter
import time



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
    
def get_ace(x,y,z,c,i):
    np.warnings.filterwarnings('ignore')
#    print('    Start')    
    try:
        start_time = time.time()
        xminground = np.min(x[c == 2])
        xmaxground = np.max(x[c == 2])
        yminground = np.min(y[c == 2])
        ymaxground = np.max(y[c == 2])
        
        checkpoint = time.time()
#        print('    Split in-range data')    
        infilt = (((x >= xminground) & (x <= xmaxground)) & ((y >= yminground) & (y <= ymaxground)))
        xin = x[infilt]
        yin = y[infilt]
        zin = z[infilt]
        cin = c[infilt]
        iin = i[infilt]
    #    ein = e[infilt]
    #    nin = n[infilt]
        
        outfilt = np.logical_not(infilt)
        xout = x[outfilt]
        yout = y[outfilt]
        zout = z[outfilt]
        cout = c[outfilt]
        iout = i[outfilt]
        
        
        groundclass = 2
        res = 1
        
#        print("    %s seconds" % (time.time() - checkpoint))
        
        checkpoint = time.time()
        # Grid the ground points, 1m grids
#        print('    Create DEM')
        gx = x[c == groundclass]
        gy = y[c == groundclass]
        gz = z[c == groundclass]
        
        grid = getRaster(gx, gy, gz, res, 'median', fillValue = -999, time = [])
        
        # Set -999 to nan
        dem = grid.grid
        dem[dem == -999] = np.nan
    
    #    print("    %s seconds" % (time.time() - checkpoint))
        
    #    checkpoint = time.time
#        print('    Prepare for Interpolation')    
        # Derive x and y indexes from mesh-grid
        x_indx, y_indx = np.meshgrid(np.arange(0, dem.shape[1]),
                                 np.arange(0, dem.shape[0]))
        
        # Mask all invalid values
        zs_masked = np.ma.masked_invalid(dem)
        
        # Retrieve the valid, non-Nan, defined values
        valid_xs = x_indx[~zs_masked.mask]
        valid_ys = y_indx[~zs_masked.mask]
        valid_zs = zs_masked[~zs_masked.mask]
        
#        print("    %s seconds" % (time.time() - checkpoint))
        
        checkpoint = time.time()
        # Generate interpolated array of z-values
#        print('    Interpolate Missing Values')
        dem = interpolate.griddata((valid_xs, valid_ys), valid_zs.ravel(),
                                         (x_indx, y_indx), method='linear')
        del valid_xs
        del valid_ys
        del valid_zs
        del x_indx
        del y_indx
        del zs_masked
        
        
        xg = grid.x
        yg = grid.y
        zg = dem
        
        xg = xg.flatten()
        yg = yg.flatten()
        zg = zg.flatten()
        
        xg = np.reshape(xg,[len(xg),1])
        yg = np.reshape(yg,[len(xg),1])
        zg = np.reshape(zg,[len(xg),1])
        
        xc = np.reshape(xin,[len(xin),1])
        yc = np.reshape(yin,[len(xin),1])
        zc = np.reshape(zin,[len(xin),1])
        cc = np.reshape(cin,[len(xin),1])
        ic = np.reshape(iin,[len(xin),1])
        
#        print("    %s seconds" % (time.time() - checkpoint))
        
        checkpoint = time.time()    
#        print('    Voxelize, Part 1')
        _, dict1, _, _ = vox.voxelize_dict(xc,yc,zc,1.0)
        gvind, gdict1, _, gdigits  = vox.voxelize_dict(xg,yg,zg,1.0)
#        print("    %s seconds" % (time.time() - checkpoint))
        checkpoint = time.time()    
#        print('    Voxelize, Part 2')
        refl = [vox.ijk_to_voxelind(list(gvind[i]),gdigits) for i in range(0,len(xg))]
        
        indexc = itemgetter(*refl)(dict1)
        indexg = itemgetter(*refl)(gdict1)
#        print("    %s seconds" % (time.time() - checkpoint))
        
        checkpoint = time.time()    
#        print('    Reassign unclassed values')   
        cc[cc == 13] = 0
        cc[cc == 3] = 4
        cc[cc == 5] = 4
        cc = acer.aceReassign(indexc, indexg, zc, cc, zg)
        xc = xc.flatten()
        yc = yc.flatten()
        zc = zc.flatten()
        cc = cc.flatten()
        ic = ic.flatten()
    #    ein = ein.flatten()
    #    ein = ein.flatten()
        xc = np.concatenate((xc, xout))
        yc = np.concatenate((yc, yout))
        zc = np.concatenate((zc, zout))
        cc = np.concatenate((cc, cout))
        ic = np.concatenate((ic, iout)) 
    #    ein = np.concatenate((ein, eout)) 
    #    nin = np.concatenate((nin, nout)) 
        return xc, yc, zc,cc,ic
    except:
        print('ACE Failed, returning original values')
        return x,y,z,c,i
        
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
                    
                    # Compute ACE
                    subData_x, subData_y, subData_z, subData_classification,subData_intensity = get_ace(subData_x,subData_y,subData_z,subData_classification,subData_intensity)
                    
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
    
    print('GET ATL TRUTH SWATH (UNIT TEST):\n')
    
    ##### Start Inputs for getAtlMeasuredSwath

    # Path to ATL03 Input File
    # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20181126114738_08990103_001_01.h5' # FINLAND
    # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20181030110205_04860106_001_01.h5' # SONONMA
    # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20190101195003_00670202_001_01_sreq_2696.h5' # SONOMA
    # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20190228170214_09510202_001_02_sreq_2696.h5' # SONOMA
    # atl03FilePath = '//bigtex/laserpewpew/data/release/001/ATL03_r001/ATL03_20190426213703_04370308_001_01.h5' # Brazil    
    
    # Path to ATL08 Input File
    # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20181126114738_08990103_952_01.h5' # FINLAND
    # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20181030110205_04860106_952_01.h5' # SONOMA
    # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20190101195003_00670202_952_01.h5' # SONOMA
    # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20190228170214_09510202_952_02.h5' # SONOMA
    # atl08FilePath = False    
    
    # Path to Output Directory
    # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_finland_20181126_python' # FINLAND
    # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_sonoma_20181030_python' # SONOMA
    # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_sonoma_20190101_python' # SONOMA
    # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_sonoma_20190228_python' # SONOMA
    # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test' # SONOMA
    
#    if os.name == 'nt':
#        base03 = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/' 
#        base08 = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/'
#        outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test'
#    atl03FilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/ATL03_20181030110205_04860106_001_01_sub_218.h5' # WSMR
#    atl08FilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/ATL08_20181030110205_04860106_001_01_sub_218.h5' # WSMR
#    outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test''

#    atl03FilePath = '/LIDAR/server/USERS/eric/2_production/icesat2_icepy/ATL03_20181030110205_04860106_001_01_sub_218.h5' # WSMR
#    atl08FilePath = '/LIDAR/server/USERS/eric/2_production/icesat2_icepy/ATL08_20181030110205_04860106_001_01_sub_218.h5' # WSMR
#    outFilePath = '/LIDAR/server/USERS/eric/2_production/other/test'
#   
#    atl03FilePath = '/laserpewpew/data/release/002/ATL03_r002/Indiana/ATL03_20181021204535_03550102_002_01.h5'
#    atl08FilePath = '/laserpewpew/data/release/002/ATL08_r002/Indiana/ATL08_20181021204535_03550102_002_01.h5'

# Current brad
#    atl03FilePath = 'Z:/data/release/002/ATL03_r002/Finland/ATL03_20181016133637_02740103_002_01.h5'
#    atl08FilePath = 'Z:/data/release/002/ATL08_r002/Finland/ATL08_20181016133637_02740103_002_01.h5'

    atl03FilePath = 'Z:/data/release/002/ATL03_r002/Indiana/ATL03_20181016090944_02710106_002_01.h5'
    atl08FilePath = 'Z:/data/release/002/ATL08_r002/Indiana/ATL08_20181016090944_02710106_002_01.h5'
    
    # outFilePath = '/LIDAR/server/USERS/eric/2_production/other/test'
    outFilePath = 'C:/Users/malonzo/GLAM/Python/PhoREAL_GUI/GUI_3.0/dist2'

    # Ground track number to analyze
    gtNum = 'gt2r'
    
    # User options
    trimInfo = 'manual,time,0,3'   # OPTIONS: ('none', 'auto', or 'manual')
                        # none - does not trim data at all
                        # auto - trims ATL03 track to extent of truth bounding region
                        # manual - trims ATL03 track by latitude or time
                            # Example: 'manual,lat,38,39'
                            # Only uses data between latitude 38 and 39 deg
                            # Example: 'manual,time,3,4'
                            # Only uses data between time 3 and 4 seconds
    
    createAtl03LasFile = True    # Option to create output measured ATL03 .las file
    createAtl03KmlFile = True    # Option to create output measured ATL03 .kml file
    createAtl08KmlFile = True    # Option to create output measured ATL08 .kml file
    createAtl03CsvFile = True    # Option to create output measured ATL03 .csv file
    createAtl08CsvFile = True    # Option to create output measured ATL08 .csv file
        
    ##### End Inputs for getAtlMeasuredSwath
    
    
    ##### Start Inputs for getAtlTruthSwath
    
    buffer = 50                 # Distance in cross-track (meters) around ATL03 track to look for truth data 
    useExistingTruth = True     # Option to use existing truth data if it exists
#    truthSwathDir = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/sonoma_ATL03_20181030110205_04860106_001_01_gt2r_TRUTH_50L50Rm_buffer_classified_short.las'
#    truthSwathDir = '/LIDAR/server/USERS/eric/2_production/icesat2_icepy/sonoma_ATL03_20181030110205_04860106_001_01_gt2r_TRUTH_50L50Rm_buffer_classified_short_dtm_utm.tif'
#    truthSwathDir = 'LIDAR/data/TanDEMX/MississippiValley_TDR_DSM/TDT_N37W085_02_DEM.TIF'
#    truthSwathDir = 'N:/data/TanDEMX/S_finland_TDR DSM_20190109.1457/ULD_TanDEM-X_TDR_DSM_20190109.1457_N62E022.tif'
    truthSwathDir = 'N:/data/TanDEMX//MississippiValley_TDR_DSM/Indiana/TDT_N37W086_02_DEM.TIF'
    createTruthFile = True      # Option to create output truth .las file

    
    ##### End Inputs for getAtlTruthSwath
    
    
    ##### CODE BELOW -- DO NOT EDIT ###############################################
    
    timeStart = runTime.time()
    
    # Call getAtlMeasuredSwath
    print('RUNNING getAtlMeasuredSwath...\n')
    atl03Data, atl08Data, headerData, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, gtNum, trimInfo, createAtl03LasFile, createAtl03KmlFile, createAtl08KmlFile, createAtl03CsvFile, createAtl08CsvFile)
        
    # Call getAtlTruthSwath
    print('RUNNING getAtlTruthSwath...\n')
    atlTruthData = getAtlTruthSwath(atl03Data, headerData, rotationData, useExistingTruth, truthSwathDir, buffer, outFilePath, createTruthFile)
            
    # End timer
    timeEnd = runTime.time()
    timeElapsedTotal = timeEnd - timeStart
    timeElapsedMin = np.floor(timeElapsedTotal / 60)
    timeElapsedSec = timeElapsedTotal % 60
        
    # Print completion message
    print('   Script Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
    print('\n')

# endIf