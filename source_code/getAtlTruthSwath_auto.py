# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 11:18:58 2019

@author: malonzo
"""

# Import Python modules
import os
import numpy as np
import time as runTime
import ntpath

# Import ICESat-2 modules
from getAtlMeasuredSwath_auto import getAtlMeasuredSwath
from icesatIO import (atlTruthStruct, writeLas, writeLog, \
                      getTruthFilePaths, getTruthHeaders, \
                      reprojectHeaderData, findMatchingTruthFiles, \
                      loadTruthFile, makeBuffer)


# Get ATL Truth Swath
def getAtlTruthSwath(atlMeasuredData, rotationData, truthHeaderDF, truthFilePaths,
                     buffer, outFilePath, createTruthFile, truthFileType, 
                     useExistingTruth, logFileID = False):
        
    # Start timer
    timeStart = runTime.time()
    
    # Initialize data
    atlTruthData = []
    
    # Print message
    gtNum = atlMeasuredData.gtNum
    beamNum = atlMeasuredData.beamNum
    beamStrength = atlMeasuredData.beamStrength
    writeLog('   Ground Track Number: %s (Beam #%s, Beam Strength: %s)\n' %(gtNum, beamNum, beamStrength), logFileID)
        
    if(useExistingTruth):
                        
        # Get truth swath .las file
        truthFilePath = truthFilePaths[0]
        
        # Get truth file header info
        writeLog('   Reading in reference buffer file: %s\n' %truthFilePath, logFileID)
        atlTruthData = loadTruthFile(truthFilePath, atlMeasuredData, rotationData, 
                                     truthFileType, outFilePath, logFileID)
                        
    else:
        
        # Reproject unmatching ESPG files to match ICESat-2 EPSG code
        writeLog('   Reprojecting header file data...', logFileID)
        truthHeaderNewDF = reprojectHeaderData(truthHeaderDF, atlMeasuredData, logFileID)
        
        # Find truth files that intersect ICESat-2 track
        writeLog('   Determining which reference files intersect ground track...', logFileID)
        _, matchingTruthFileInds = findMatchingTruthFiles(truthHeaderNewDF, atlMeasuredData, rotationData, buffer)
        matchingTruthFiles = np.array(truthFilePaths)[matchingTruthFileInds]
        
        # Read truth files that intersect ICESat-2 track
        if(len(matchingTruthFiles)>0):
        
            # Print messages
            writeLog('   Ground track crosses over %d (out of %d) reference files' %(len(matchingTruthFiles), len(truthHeaderDF)), logFileID)
            writeLog('   Reading from directory: %s' %ntpath.dirname(matchingTruthFiles[0]), logFileID)
            writeLog('', logFileID)
            writeLog('   Reading File:', logFileID)
            writeLog('   -------------', logFileID)
            
            # Initialize parameters
            fileNum = 1
            atlTruthData = atlTruthStruct([],[],[],[],[],[],[],[],[],[],[])
        
            # Loop over matching files
            for i in range(0,len(matchingTruthFiles)):
                
                # Get truth file to read
                truthFilePath = matchingTruthFiles[i]
                
                # Get truth base file name
                baseName = ntpath.basename(truthFilePath)
                
                # Read truth file
                writeLog('   %d) %s' %(fileNum, baseName), logFileID)
                atlTruthDataSingle = loadTruthFile(truthFilePath, atlMeasuredData, rotationData, truthFileType, outFilePath, logFileID)
                  
                # Get truth file buffer
                if(bool(atlTruthDataSingle)):
                    
                    writeLog('      Buffering data...', logFileID)
                    atlTruthDataBuffer = makeBuffer(atlTruthDataSingle, atlMeasuredData, rotationData, buffer)
                
                    # Append truth files
                    atlTruthData.append(atlTruthDataBuffer)
    
                    # Increment counter
                    fileNum += 1
                
                # endIf
                
            # endFor
        
        else:
            
            # No data to process
            atlTruthData = atlTruthStruct([],[],[],[],[],[],[],[],[],[],[])
            
            writeLog('', logFileID)
            writeLog('   WARNING: No matching reference files intersect ground track.', logFileID)
            
        # endIf
    # endIf
    
    # Test if atlTruthData is empty
    atlTruthEmpty = (len(atlTruthData.easting)==0) and (len(atlTruthData.northing)==0)
    
    # Inform user if data is empty
    if(atlTruthEmpty):
        writeLog('   WARNING: Reference data is empty.', logFileID)
    # endIf
    
    # Create output file
    if(createTruthFile and not(atlTruthEmpty)): 
        
        # Set output file name
        if(useExistingTruth):
            outName = os.path.basename(truthFilePath)
        else:
            outName = atlMeasuredData.atl03FileName + '_' + \
                      atlMeasuredData.gtNum + '_REFERENCE_' + str(buffer) + \
                      'L' + str(buffer) + 'Rm_buffer.las'
        # endIf
        
        # Set full output file name and path
        outPath = os.path.normpath(outFilePath + '/' + outName)

        # If output directory does not exist, create it
        if(not os.path.exists(os.path.normpath(outFilePath))):
            os.makedirs(os.path.normpath(outFilePath))
        # EndIf
        
        # Write to file
        writeLog('', logFileID) 
        writeLog('   Writing reference .las file...', logFileID)
        writeLas(np.ravel(atlTruthData.easting),
                 np.ravel(atlTruthData.northing),
                 np.ravel(atlTruthData.z),
                 'utm', outPath,
                 np.ravel(atlTruthData.classification),
                 np.ravel(atlTruthData.intensity),
                 None, atlTruthData.hemi,atlTruthData.zone)
    # endIf
    
    # End timer
    timeEnd = runTime.time()
    timeElapsedTotal = timeEnd - timeStart
    timeElapsedMin = np.floor(timeElapsedTotal / 60)
    timeElapsedSec = timeElapsedTotal % 60
    
    # Print completion message
    writeLog('', logFileID)
    writeLog('   Module Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec), logFileID)
    writeLog('\n', logFileID)
    
    # Return object
    return atlTruthData

# endDef
    
    
# Unit test
if __name__ == "__main__":
    
    writeLog('GET ATL TRUTH SWATH (UNIT TEST):\n')
    
    ##### Start Inputs for getAtlMeasuredSwath

    atl03FilePath = 'C:/Users/malonzo/GLAM/delete/ATL03_20181030110205_04860106_002_01.h5'
    atl08FilePath = 'C:/Users/malonzo/GLAM/delete/ATL08_20181030110205_04860106_002_01.h5'

    outFilePath = 'C:/Users/malonzo/GLAM/delete'

    # Ground track number to analyze
    gtNum = 'gt1r'
    
    # User options
    trimInfo = 'none'   # OPTIONS: ('none', 'auto', or 'manual')
                        # none - does not trim data at all
                        # auto - trims ATL03 track to extent of truth bounding region
                        # manual - trims ATL03 track by latitude or time
                            # Example: 'manual,lat,38,39'
                            # Only uses data between latitude 38 and 39 deg
                            # Example: 'manual,time,3,4'
                            # Only uses data between time 3 and 4 seconds
    
    createAtl03LasFile = True    # Option to create output measured ATL03 .las file
    createAtl03KmlFile = False    # Option to create output measured ATL03 .kml file
    createAtl08KmlFile = False    # Option to create output measured ATL08 .kml file
    createAtl03CsvFile = False    # Option to create output measured ATL03 .csv file
    createAtl08CsvFile = False    # Option to create output measured ATL08 .csv file
        
    ##### End Inputs for getAtlMeasuredSwath
    
    
    ##### Start Inputs for getAtlTruthSwath
    
    buffer = 50                 # Distance in cross-track (meters) around ATL03 track to look for truth data 
    useExistingTruth = False     # Option to use existing truth data if it exists

    truthSwathDir = 'Z:/data/validation/data/WSMR/ir'
    truthFileType = '.las'
    createTruthFile = True      # Option to create output truth .las file

    ##### End Inputs for getAtlTruthSwath
    
    
    ##### CODE BELOW -- DO NOT EDIT ###############################################
    
    timeStart = runTime.time()
    
    # Call getAtlMeasuredSwath
    print('RUNNING getAtlMeasuredSwath...\n')
    atl03Data, atl08Data, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, gtNum, trimInfo, createAtl03LasFile, createAtl03KmlFile, createAtl08KmlFile, createAtl03CsvFile, createAtl08CsvFile)

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
    atlTruthData = getAtlTruthSwath(atl03Data, rotationData, 
                                    truthHeaderDF, truthFilePaths,
                                    buffer, outFilePath, createTruthFile, 
                                    truthFileType, useExistingTruth, logFileID=False)
                
    # End timer
    timeEnd = runTime.time()
    timeElapsedTotal = timeEnd - timeStart
    timeElapsedMin = np.floor(timeElapsedTotal / 60)
    timeElapsedSec = timeElapsedTotal % 60
        
    # Print completion message
    print('   Script Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
    print('\n')

# endIf