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
                      loadTruthFile, makeBuffer, offsetsStruct)
from icesatUtils import (superFilter)
from icesatCalVal import (perfectClassifier)
from getMeasurementError_auto import getMeasurementError, offsetsStruct


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
            atlTruthData = atlTruthStruct([],[],[],[],[],[],[],[],[],[],[],[],[],[])
        
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
            atlTruthData = atlTruthStruct([],[],[],[],[],[],[],[],[],[],[],[],[],[])
            
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

    atl03FilePath = 'C:/Users/malonzo/GLAM/delete/ATL03_20181030110205_04860106_002_01.h5'
    atl08FilePath = 'C:/Users/malonzo/GLAM/delete/ATL08_20181030110205_04860106_002_01.h5'
    
    atl03FilePath = 'C:/Users/malonzo/GLAM/delete/ATL03_20191129150021_09730506_003_01_sreq_3257.h5'
    atl08FilePath = []

    atl03FilePath = 'E:/0_data/is2/prf/ATL03/ATL03_20190505104857_05680302_005_01.h5'
    atl08FilePath = 'E:/0_data/is2/prf/ATL08/ATL08_20190505104857_05680302_005_01.h5'

    atl03FilePath = 'E:/0_data/is2/prf/ATL03/ATL03_20191005033230_01260502_005_01.h5'
    atl08FilePath = 'E:/0_data/is2/prf/ATL08/ATL08_20191005033230_01260502_005_01.h5'
    
    atl03_folder = 'E:/0_data/is2/prf/ATL03/'
    atl08_folder = 'E:/0_data/is2/prf/ATL08/'
    out_folder = 'E:/0_data/is2/prf/'    
    import os 
    
    atl03_list = os.listdir(atl03_folder)
    atl08_list = os.listdir(atl08_folder)
    
    # for i in range(0,len(atl03_list)):
    # gt_list = ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']
    # for i in range(0,len(atl03_list)):
    i = 2
    atl03FilePath = os.path.join(atl03_folder, atl03_list[i])
    atl08FilePath = os.path.join(atl08_folder, atl08_list[i])


    # outFilePath = '/LIDAR/server/USERS/eric/2_production/other/test'
    outFilePath = 'E:/output'

    # Ground track number to analyze
    gtNum = 'gt1r'
    # for gtNum in gt_list:
        
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
#    truthSwathDir = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/sonoma_ATL03_20181030110205_04860106_001_01_gt2r_TRUTH_50L50Rm_buffer_classified_short.las'
#    truthSwathDir = '/LIDAR/server/USERS/eric/2_production/icesat2_icepy/sonoma_ATL03_20181030110205_04860106_001_01_gt2r_TRUTH_50L50Rm_buffer_classified_short_dtm_utm.tif'
#    truthSwathDir = 'LIDAR/data/TanDEMX/MississippiValley_TDR_DSM/TDT_N37W085_02_DEM.TIF'
#    truthSwathDir = 'N:/data/TanDEMX/S_finland_TDR DSM_20190109.1457/ULD_TanDEM-X_TDR_DSM_20190109.1457_N62E022.tif'
#    truthSwathDir = 'C:/Users/malonzo/GLAM/delete'
    truthSwathDir = 'E:/data/2018spl_2959_6647'
    # truthSwathDir = 'E:/data/2019_spl_3800m'
    # truthSwathDir = 'E:/data/2019spl_2000m'

#    truthSwathDir = 'C:/Users/malonzo/GLAM/delete/ATL03_20191129150021_09730506_003_01_sreq_3257_gt1r_REFERENCE_50L50Rm_buffer.las'
    truthFileType = '.laz'
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
    
    print('RUNNING getMeasurementError...\n')
# atlCorrections = getMeasurementError(atl03Data, atlTruthData, rotationData, outFilePath, useMeasSigConf, filterData, offsets, createMeasCorrFile, makePlots, showPlots)

    refHeightType = 'hae'
    
    offsetsCrossTrackBounds = np.array([-50,50])      # Cross-track limits to search for geolocation error
    offsetsAlongTrackBounds = np.array([-50,50])      # Along-track limits to search for geolocation error
    offsetsRasterResolutions = np.array([8, 4, 2, 1])  # Multi-resolutional step-down raster resolutions (in meters)
    refHeightType = 'HAE'              # 'HAE' or 'MSL'
    offsetsUseVerticalShift = False    # Option to use a vertical shift
    offsetsVerticalShift = 0           # Vertical shift to use if above set to True (in meters)
    useMeasSigConf = False             # Use measured signal confidence (or use ground truth)
                                     # Meas Classes (0 = Unclass, 1 = Ground, 2 = Low Veg, 3 = High Veg), Texpert Truth Classes (0 = Unclass, 2 = Ground, 4 = Veg, 6 = Building)
    filterData = [1]               # Signal Confidence (0, 1, 2, 3, 4)
    offsets = offsetsStruct(offsetsCrossTrackBounds, offsetsAlongTrackBounds, offsetsRasterResolutions, offsetsUseVerticalShift, offsetsVerticalShift )
    createMeasCorrFile = True     # Option to create ouput measured corrected .las file
    makePlots = True              # Option to make output plots
    showPlots = False             # Option to show output plot windows
    
    atlCorrections = getMeasurementError(atl03Data, atlTruthData, refHeightType, 
                            rotationData, outFilePath, 
                            useMeasSigConf, filterData, offsets, createMeasCorrFile, 
                            makePlots, showPlots, logFileID = False)
    
    # Store all objects into one
    # atlMeasuredDataAll.append(atl03Data)
    # atlTruthDataAll.append(atlTruthData)
    # atlCorrectionsAll.append(atlCorrections)
        
    # EndFor
    
    # End timer
    timeEnd = runTime.time()
    timeElapsedTotal = timeEnd - timeStart
    timeElapsedMin = np.floor(timeElapsedTotal / 60)
    timeElapsedSec = timeElapsedTotal % 60
        
    # Print completion message
    print('   Geolocation offset Complete %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
    print('\n')
    
    atlTruthData.alongTrack = atlTruthData.alongTrack - atlCorrections.alongTrack
    atlTruthData.crossTrack = atlTruthData.crossTrack - atlCorrections.crossTrack
    atlTruthData.z = atlTruthData.z - atlCorrections.z

    atlTruthData.year = np.r_['c',(atlTruthData.year)]
    atlTruthData.month = np.r_['c',(atlTruthData.month)]
    atlTruthData.day = np.r_['c',(atlTruthData.day)]
    
    atl03Data, atlTruthData = superFilter(atl03Data, atlTruthData,
                                          xBuf = 6, classCode = [], verbose=False)
    
    measpc, measoc = perfectClassifier(atl03Data, atlTruthData,ground = [1],canopy = [2,3], 
                      unclassed = [0], keepsize = True)
    
    # Make Pandas DF with everything
    ##Filter Truth Data by Include
    alongtrack = atl03Data.alongTrack.flatten()
    crosstrack = atl03Data.crossTrack.flatten()
    z = atl03Data.z.flatten()
    # easting = atl03Data.easting.flatten()
    # northing = atl03Data.northing.flatten()
    classification = atl03Data.classification.flatten()    
    
    # Export Pandas DF as PC
    import pandas as pd
    df_meas = pd.DataFrame(z,columns=['z'])
    df_meas = pd.concat([df_meas,pd.DataFrame(
        crosstrack,columns=['crosstrack'])],axis=1)
    df_meas = pd.concat([df_meas,pd.DataFrame(
        alongtrack,columns=['alongtrack'])],axis=1)
    # df_truth = pd.concat([df_truth,pd.DataFrame(
    #     easting,columns=['easting'])],axis=1)
    # df_truth = pd.concat([df_truth,pd.DataFrame(
    #     northing,columns=['northing'])],axis=1)
    df_meas = pd.concat([df_meas,pd.DataFrame(
        classification,columns=['label'])],axis=1)
    df_meas = pd.concat([df_meas,pd.DataFrame(
        measpc,columns=['truth_label'])],axis=1)
    out_name = atl03_list[0].split('.h5')[0] + '_' + gtNum + '.csv'
    out_file = os.path.join(out_folder,out_name)
    df_meas.to_csv(out_file)

