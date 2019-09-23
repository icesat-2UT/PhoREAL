# -*- coding: utf-8 -*-
"""
Script to perform most basic functionalities required to get ATL03 swath

Copyright 2019 Applied Research Laboratories, University of Texas at Austin

This package is free software; the copyright holder gives unlimited
permission to copy and/or distribute, with or without modification, as
long as this notice is preserved.

Authors:
    Mike Alonzo
    Eric Guenther
    
Date: September 20, 2019

"""

# Import modules
import os
import numpy as np
import time as runTime
from icesatIO import (readAtl03H5, readAtl03DataMapping, readAtl08DataMapping, 
                      readTruthRegionsTxtFile, readHeaderMatFile, 
                      writeLas, writeKml, writeArrayToCSV)
from icesatUtils import (getNameParts, getAtl08Mapping, getLatLon2UTM, 
                         getCoordRotFwd, getClosest)


class atlMeasuredStruct:
        
    # Define class with designated fields
    def __init__(self, lat, lon, easting, northing, crossTrack, alongTrack, z, 
                 time, signalConf, classification, gtNum, zone, hemi, 
                 kmlRegionName, headerFilePath, truthFilePath, atl03FilePath, 
                 atl03FileName, atl08FilePath, atl08FileName, trackDirection, 
                 alt03h5Info, dataIsMapped):
            
        self.lat = np.c_[lat]
        self.lon = np.c_[lon]
        self.easting = np.c_[easting]
        self.northing = np.c_[northing]
        self.crossTrack = np.c_[crossTrack]
        self.alongTrack = np.c_[alongTrack]
        self.z = np.c_[z]
        self.time = np.c_[time]
        self.signalConf = np.c_[signalConf]
        self.classification = np.c_[classification]
        self.intensity = np.zeros(np.shape(np.c_[classification]))
        self.gtNum = gtNum
        self.zone = zone
        self.hemi = hemi
        self.kmlRegionName = kmlRegionName
        self.headerFilePath = headerFilePath
        self.truthFilePath = truthFilePath
        self.atl03FilePath = atl03FilePath
        self.atl03FileName = atl03FileName
        self.atl08FilePath = atl08FilePath
        self.atl08FileName = atl08FileName
        self.trackDirection = trackDirection
        self.atlVersion = alt03h5Info.atlVersion
        self.year = alt03h5Info.year
        self.month = alt03h5Info.month
        self.day = alt03h5Info.day
        self.hour = alt03h5Info.hour
        self.minute = alt03h5Info.minute
        self.second = alt03h5Info.second
        self.trackNum = alt03h5Info.trackNum
        self.unknown = alt03h5Info.unknown
        self.releaseNum = alt03h5Info.releaseNum
        self.incrementNum = alt03h5Info.incrementNum
        self.dataIsMapped = dataIsMapped
 
        
class atlRotationStruct:
    
    # Define class with designated fields
    def __init__(self, R_mat, xRotPt, yRotPt, desiredAngle, phi):
        
        self.R_mat = R_mat
        self.xRotPt = xRotPt
        self.yRotPt = yRotPt
        self.desiredAngle = desiredAngle
        self.phi = phi
        

def getAtlMeasuredSwath(atl03FilePath = False, atl08FilePath = False, 
                        outFilePath = False, gtNum = 'gt1r', trimInfo = 'auto', 
                        createLasFile = False, createKmlFile = False, 
                        createCsvFile = False):
    
    # Initialize outputs
    atlMeasuredData = False
    headerData = False
    coordType = False
    rotationData = False
    
    # Only execute code if input ATL03 .h5 file and output path declared
    if(atl03FilePath and outFilePath):
    
        # Start timer
        timeStart = runTime.time()
        
        # Print message
        print('   Ground Track Number: %s' % gtNum)
        
        # Get ATL03 file path/name
        atl03FilePath = os.path.normpath(os.path.abspath(atl03FilePath))
        atl03FileName = os.path.splitext(os.path.basename(atl03FilePath))[0]
        
        # Read ATL03 data from h5 file
        print('   Reading ATL03 .h5 file: %s' % atl03FilePath)
        lat_all = readAtl03H5(atl03FilePath, 'lat_ph', gtNum)
        lon_all = readAtl03H5(atl03FilePath, 'lon_ph', gtNum)
        z_all = readAtl03H5(atl03FilePath, 'h_ph', gtNum)
        deltaTime_all = readAtl03H5(atl03FilePath, 'delta_time', gtNum)
        signalConf_all = readAtl03H5(atl03FilePath, 'signal_conf_ph', gtNum)
        atl03_ph_index_beg, atl03_segment_id = readAtl03DataMapping(atl03FilePath, gtNum)
            
        # Get time from delta time
        time_all = deltaTime_all - np.min(deltaTime_all)
        
        # Get track direction
        if(lat_all[-1] > lat_all[0]):
            trackDirection = 'Ascending'
        else:
            trackDirection = 'Descending'
        # endIf
        
        # Extract metadata from ATL03 file name
        atl03h5Info = getNameParts(atl03FileName)
        
        # Map ATL08 to ATL03 ground photons
        if(atl08FilePath):
        
            # Get ATL08 file path/name
            atl08FilePath = os.path.normpath(os.path.abspath(atl08FilePath))
            
            # Message to screen
            print('   Reading ATL08 .h5 file: %s' % atl08FilePath)  
            atl08_classed_pc_indx, atl08_classed_pc_flag, atl08_segment_id = readAtl08DataMapping(atl08FilePath, gtNum)
    
            print('   Mapping ATL08 to ATL03 Ground Photons...')
            classification = getAtl08Mapping(atl03_ph_index_beg, atl03_segment_id, atl08_classed_pc_indx, atl08_classed_pc_flag, atl08_segment_id)
            
            # Get length to trim data by
            class_length = len(classification)
            lat_length = len(lat_all)
            data_length = np.min([class_length, lat_length]) 
            
            # Trim ATL03 data down to size of classification array
            lat = lat_all[0:data_length]
            lon = lon_all[0:data_length]
            z = z_all[0:data_length]
            time = time_all[0:data_length]
            signalConf = signalConf_all[0:data_length]
            classification = classification[0:data_length]
            dataIsMapped = True
        
        else:
            
            # Message to screen
            print('   Not Mapping ATL08 to ATL03 Ground Photons')
            
            # Store data
            lat = lat_all
            lon = lon_all
            z = z_all
            time = time_all
            signalConf = signalConf_all
            classification = np.zeros(np.size(lat_all))
            dataIsMapped = False
            
        # endIf
        
        # Get trim options
        trimParts = trimInfo.split(',')
        trimMode = trimParts[0]
        trimType = 'None'
        if(('manual' in trimMode.lower()) and (len(trimParts) > 1)):
            trimType = trimParts[1]
            trimMin = float(trimParts[2])
            trimMax = float(trimParts[3])
        # endIf
        
        
        # If selected to manually trim data, do this first
        if('manual' in trimMode.lower()):
            
            # Trim by lat or time
            if('lat' in trimType.lower()):
                print('   Manual Trim Mode (Min Lat: %s, Max Lat: %s)' % (trimMin, trimMax))
                indsToKeep = (lat >= trimMin) & (lat <= trimMax)
            elif('time' in trimType.lower()):
                print('   Manual Trim Mode (Min Time: %s, Max Time: %s)' % (trimMin, trimMax))
                indsToKeep = (time >= trimMin) & (time <= trimMax)
            else:
                print('   Manual Trim Mode is Missing Args, Not Trimming Data')
                indsToKeep = np.ones(np.shape(lat), dtype = bool)
            # endif
            
            # Trim data
            lat = lat[indsToKeep]
            lon = lon[indsToKeep]
            z = z[indsToKeep]
            time = time[indsToKeep]
            signalConf = signalConf[indsToKeep]
            classification = classification[indsToKeep]
            
        # EndIf
        
        # Determine if ATL03 track goes over a lidar truth region
        kmlBoundsTextFile = 'kmlBounds.txt'
        kmlRegionName = False
        headerFilePath = False
        truthFilePath = False
        if(os.path.exists(kmlBoundsTextFile)):
            
            # Message to user
            print('   Finding Truth Region...')
            
            try:
                
                # Read kmlBounds.txt file and get contents
                kmlInfo = readTruthRegionsTxtFile(kmlBoundsTextFile)
                
                # Loop through kmlBounds.txt and find matching TRUTH area
                maxCounter = len(kmlInfo.regionName)
                counter = 0
                while(not kmlRegionName):
                    latInFile = (lat >= kmlInfo.latMin[counter]) & (lat <= kmlInfo.latMax[counter])
                    lonInFile = (lon >= kmlInfo.lonMin[counter]) & (lon <= kmlInfo.lonMax[counter])
                    trackInRegion = any(latInFile & lonInFile)
                    if(trackInRegion):
                        
                        # Get truth region info
                        kmlRegionName = kmlInfo.regionName[counter]
                        headerFilePath = os.path.normpath(kmlInfo.headerFilePath[counter])
                        truthFilePath = os.path.normpath(kmlInfo.truthFilePath[counter])
                        kmlLatMin = kmlInfo.latMin[counter]
                        kmlLatMax = kmlInfo.latMax[counter]
                        kmlLonMin = kmlInfo.lonMin[counter]
                        kmlLonMax = kmlInfo.lonMax[counter]
                        
                        # Print truth region
                        print('   Truth File Region: %s' % kmlRegionName)
                    
                        # Read truth header .mat file
                        headerData = readHeaderMatFile(headerFilePath)
                        coordType = headerData.coordType
                    
                    # endIf
                    
                    # Increment counter
                    counter += 1
                    
                    if(counter >= maxCounter):
                        
                        # Send message to user
                        print('   No Truth File Region Found in kmlBounds.txt')
                        break
                    
                    # endIf
                # endWhile
            
            except:
                
                # Could not read kmlBounds.txt file
                print('   Could not read truth header .mat file. Auto-assigning UTM zone.')
            
        # endIf
        
        # If selected to auto trim data, then trim if truth region exists
        if('auto' in trimMode.lower()):
            
            if(kmlRegionName):
                
                # Trim data based on TRUTH region
                print('   Auto Trim Mode (Trimming Data Based on Truth Region)...')
                indsInRegion = (lat >= kmlLatMin) & (lat <= kmlLatMax) & (lon >= kmlLonMin) & (lon <= kmlLonMax)
                lat = lat[indsInRegion]
                lon = lon[indsInRegion]
                z = z[indsInRegion]
                time = time[indsInRegion]
                signalConf = signalConf[indsInRegion]
                classification = classification[indsInRegion]
                
            # endIf
        # endIf
        
        # Convert lat/lon coordinates to UTM
        print('   Converting Lat/Lon to UTM...')
        if(kmlRegionName and coordType == 'UTM'):
                
            # Force UTM zone to match region header file
            zoneToUse = headerData.zone
            hemiToUse = headerData.hemi
            easting, northing, zone, hemi = getLatLon2UTM(lon, lat, zoneToUse, hemiToUse)
            
        else:
                
            # Allow function to determine UTM zone
            easting, northing, zone, hemi = getLatLon2UTM(lon, lat) 
                
        # endIf
        
        # Print UTM zone
        print('   UTM Zone: %s' % zone)
      
        # Transform MEASURED data to CT/AT plane
        print('   Computing CT/AT Frame Rotation...')
        desiredAngle = 90
        crossTrack, alongTrack, R_mat, xRotPt, yRotPt, phi = getCoordRotFwd(easting, northing, [], [], [], desiredAngle)
        
        # Store rotation object
        rotationData = atlRotationStruct(R_mat, xRotPt, yRotPt, desiredAngle, phi)
        
        # Store data in class structure
        atlMeasuredData = atlMeasuredStruct(lat, lon, \
                                            easting, northing, \
                                            crossTrack, alongTrack, \
                                            z, time, \
                                            signalConf, classification, \
                                            gtNum, zone, hemi, \
                                            kmlRegionName, headerFilePath, truthFilePath, \
                                            atl03FilePath, atl03FileName, \
                                            atl08FilePath, atl03FileName, \
                                            trackDirection, \
                                            atl03h5Info, \
                                            dataIsMapped)
        
        # Create output .las file
        if(createLasFile): 
            
            print('   Writing measured .las file...')
            outName = atlMeasuredData.atl03FileName + '_' + atlMeasuredData.gtNum + '_MEASURED.las'
            outPath = os.path.normpath(outFilePath + '/' + outName)
            
            # If output directory does not exist, create it
            if(not os.path.exists(os.path.normpath(outFilePath))):
                os.mkdir(os.path.normpath(outFilePath))
            # EndIf
            
            # Write .las file
            writeLas(np.ravel(atlMeasuredData.easting),np.ravel(atlMeasuredData.northing),np.ravel(atlMeasuredData.z),'utm',outPath,np.ravel(atlMeasuredData.classification),np.ravel(atlMeasuredData.intensity),atlMeasuredData.hemi,atlMeasuredData.zone)
            
        # endIf
        
        # Create output .kml file
        if(createKmlFile):
            
            print('   Writing measured .kml file...')
            outName = atlMeasuredData.atl03FileName + '_' + atlMeasuredData.gtNum + '_MEASURED.kml'
            outPath = os.path.normpath(outFilePath + '/' + outName)
            
            # If output directory does not exist, create it
            if(not os.path.exists(os.path.normpath(outFilePath))):
                os.mkdir(os.path.normpath(outFilePath))
            # EndIf
            
            # Get array of input time values
            timeStep = 1 # seconds
            timeVals = np.arange(np.min(atlMeasuredData.time), np.max(atlMeasuredData.time) + 1, timeStep)
        
            # Get closest time values from IceSat MEASURED data
            timeIn, indsToUse = getClosest(atlMeasuredData.time, timeVals)
        
            # Reduce lat/lon values to user-specified time scale
            lonsIn = atlMeasuredData.lon[indsToUse]
            latsIn = atlMeasuredData.lat[indsToUse]
    
            # Write .kml file
            writeKml(latsIn, lonsIn, timeIn, outPath)
        
        # endIf
        
        # Create output .csv file
        if(createCsvFile):
            
            print('   Writing measured .csv file...')
            outName = atlMeasuredData.atl03FileName + '_' + atlMeasuredData.gtNum + '_MEASURED.csv'
            outPath = os.path.normpath(outFilePath + '/' + outName)
            
            # If output directory does not exist, create it
            if(not os.path.exists(os.path.normpath(outFilePath))):
                os.mkdir(os.path.normpath(outFilePath))
            # EndIf
            
            # Write .csv file
            namelist = ['Time (sec)', 'Latitude (deg)', 'Longitude (deg)', \
                        'Easting (m)', 'Northing (m)', \
                        'Cross-Track (m)', 'Along-Track (m)', \
                        'Height (m)', \
                        'Classification', 'Signal Confidence']
            datalist = [atlMeasuredData.time, atlMeasuredData.lat, atlMeasuredData.lon, \
                        atlMeasuredData.easting, atlMeasuredData.northing, \
                        atlMeasuredData.crossTrack, atlMeasuredData.alongTrack,\
                        atlMeasuredData.z, \
                        atlMeasuredData.classification, atlMeasuredData.signalConf] 
            writeArrayToCSV(outPath, namelist, datalist)
            
        # endIf
        
        # End timer
        timeEnd = runTime.time()
        timeElapsedTotal = timeEnd - timeStart
        timeElapsedMin = np.floor(timeElapsedTotal / 60)
        timeElapsedSec = timeElapsedTotal % 60
        
        # Print completion message
        print('   Module Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
        print('\n')
        
    else:
        
        # Message to user
        print('Input correct ATL03 input .h5 file and/or output path.')
        
    # endIf
    
    # Return object
    return atlMeasuredData, headerData, rotationData
    
# Unit test
if __name__ == "__main__":
    
    ##### Start Inputs for getAtlMeasuredSwath

    # Path to ATL03 Input File
    # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20181126114738_08990103_001_01.h5' # FINLAND
    atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20181030110205_04860106_001_01.h5' # SONONMA
    # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20190101195003_00670202_001_01_sreq_2696.h5' # SONOMA
    # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20190228170214_09510202_001_02_sreq_2696.h5' # SONOMA
    # atl03FilePath = '//bigtex/laserpewpew/data/release/001/ATL03_r001/ATL03_20190426213703_04370308_001_01.h5' # Brazil    
    
    # Path to ATL08 Input File
    # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20181126114738_08990103_952_01.h5' # FINLAND
    atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20181030110205_04860106_952_01.h5' # SONOMA
    # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20190101195003_00670202_952_01.h5' # SONOMA
    # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20190228170214_09510202_952_02.h5' # SONOMA
    # atl08FilePath = False    
    
    # Path to Output Directory
    # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_finland_20181126_python' # FINLAND
    # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_sonoma_20181030_python' # SONOMA
    # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_sonoma_20190101_python' # SONOMA
    # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_sonoma_20190228_python' # SONOMA
    outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test' # SONOMA

    # Ground track number to analyze
    gtNum = 'gt1r'
    
    # User options
    trimInfo = 'auto'   # OPTIONS: ('none', 'auto', or 'manual')
                        # none - does not trim data at all
                        # auto - trims ATL03 track to extent of truth bounding region
                        # manual - trims ATL03 track by latitude or time
                            # Example: 'manual,lat,38,39'
                            # Only uses data between latitude 38 and 39 deg
                            # Example: 'manual,time,3,4'
                            # Only uses data between time 3 and 4 seconds
    
    createLasFile = True    # Option to create output measured .las file
    createKmlFile = True    # Option to create output measured .kml file
    createCsvFile = True    # Option to create output measured .csv file
            
    ##### End Inputs for getAtlMeasuredSwath
    
    
    ##### CODE BELOW -- DO NOT EDIT ###############################################
    
    timeStart = runTime.time()
    
    # Call getAtlMeasuredSwath
    print('RUNNING getAtlMeasuredSwath...\n')
    atlMeasuredData, headerData, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, gtNum, trimInfo, createLasFile, createKmlFile, createCsvFile)
        
    # End timer
    timeEnd = runTime.time()
    timeElapsedTotal = timeEnd - timeStart
    timeElapsedMin = np.floor(timeElapsedTotal / 60)
    timeElapsedSec = timeElapsedTotal % 60
        
    # Print completion message
    print('   Script Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
    print('\n')

# endIf