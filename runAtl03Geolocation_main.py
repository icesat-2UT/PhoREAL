# -*- coding: utf-8 -*-
"""

This is the main script used to run the ATL03 Geolocation Correction algorithm

All user inputs are declared below.

The ATL03 Geolocation Correction algorithm is composed of 3 main functions:
    1) getAtlMeasuredSwath - gets IceSat-2 ATL03 data from .h5 file
    2) getAtlTruthSwath - gets truth swath around IceSat-2 ATL03 track
    3) getMeasurementError - finds geolocation offsets of IceSat-2 data

Created on Fri Aug 30 08:40:15 2019

@author: malonzo
"""

import numpy as np
import time as runTime
from getAtlMeasuredSwath_auto import getAtlMeasuredSwath 
from getAtlTruthSwath_auto import getAtlTruthSwath
from getMeasurementError_auto import getMeasurementError, offsetsStruct


##### Start Inputs for getAtlMeasuredSwath

# Path to ATL03 Input File
# atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20181030110205_04860106_001_01.h5' # SONONMA
# atl03FilePath = 'Z:/data/release/R002/ATL03_rR002/ATL03_20190928175636_00280506_R002_01_sreq_3002.h5' # WSMR 20190928
# atl03FilePath = 'Z:/data/release/001/ATL03_r001/ATL03_20190331023705_00280306_001_02_sreq_2695.h5' # WSMR 20190331
# atl03FilePath = 'z:/data/release/R002/ATL03_rR002/ATL03_20191008052407_01730502_R002_01_sub_344.h5' # WSMR 20191008
# atl03FilePath = 'z:/data/release/R002/ATL03_rR002/ATL03_20191031162419_05310506_R002_01_sub_344.h5' # WSMR 20191008
# atl03FilePath = 'z:/data/release/R002/ATL03_rR002/ATL03_20200107010352_01730602_R002_01_sub_344.h5' # WSMR 20200107
# atl03FilePath = 'z:/data/release/002/ATL03_r002/ATL03_20190928175636_00280506_002_01_sreq_3181.h5' # WSMR 20190928 r002
# atl03FilePath = 'z:/data/release/002/ATL03_r002/ATL03_20191012051547_02340502_002_01_sreq_3181.h5' # WSMR 20191012 r002
# atl03FilePath = 'z:/data/release/002/ATL03_r002/ATL03_20191031162419_05310506_002_01_sreq_3181.h5' # WSMR 20191031 r002
atl03FilePath = 'z:/data/release/002/ATL03_r002/WSMR/ATL03_20190331023705_00280306_002_01.h5' # WSMR 20190331 r002

# Path to ATL08 Input File
# atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20181030110205_04860106_952_01.h5' # SONOMA
# atl08FilePath = 'Z:/data/release/R002/ATL08_rR002/ATL08_20190928175636_00280506_R002_01_sub_344.h5' # WSMR
atl08FilePath = 'z:/data/release/002/ATL08_r002/WSMR/ATL08_20190331023705_00280306_002_01.h5' # WSMR 20190331 r002
# atl08FilePath = []

# Path to Output Directory
# outFilePath = '//lidar-server/lidar/USERS/test'
# outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/rapid002_wsmr_20190928_old' # WSMR
# outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/r001_wsmr_20190331' # WSMR  
# outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/rapid002_wsmr_20191008' # WSMR 
# outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/rapid002_wsmr_20191031' # WSMR
# outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/rapid002_wsmr_20200107' # WSMR
# outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/r002_wsmr_20190928_new' # WSMR r002
# outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/r002_wsmr_20191012_new' # WSMR r002
# outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/r002_wsmr_20191031_new' # WSMR r002
outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/r002_wsmr_20190331' # WSMR r002

# Ground track number(s) to analyze
gtNums = ['gt2r']

# User options
trimInfo = 'auto'   # OPTIONS: ('none', 'auto', or 'manual')
                    # none - does not trim data at all
                    # auto - trims ATL03 track to extent of truth bounding region
                    # manual - trims ATL03 track by latitude or time
                        # Example: 'manual,lat,38,39'
                        # Only uses data between latitude 38 and 39 deg
                        # Example: 'manual,time,3,4'
                        # Only uses data between time 3 and 4 seconds
                                    
createAtl03LasFile = False    # Option to create output measured ATL03 .las file
createAtl03KmlFile = False    # Option to create output measured ATL03 .kml file
createAtl08KmlFile = False    # Option to create output measured ATL08 .kml file
createAtl03CsvFile = False    # Option to create output measured ATL03 .csv file
createAtl08CsvFile = False    # Option to create output measured ATL08 .csv file
    
##### End Inputs for getAtlMeasuredSwath


##### Start Inputs for getAtlTruthSwath

buffer = 50                 # Distance in cross-track (meters) around ATL03 track to look for truth data 
useExistingTruth = True     # Option to use existing truth data if it exists
truthSwathDir = '//lidar-server/lidar/USERS/mike/iceSat2/truth_data/'  # Path to existing truth data
createTruthFile = False      # Option to create output truth .las file

##### End Inputs for getAtlTruthSwath


##### Start Inputs for getMeasurementError

offsetsCrossTrackBounds = np.array([-48, 48])      # Cross-track limits to search for geolocation error
offsetsAlongTrackBounds = np.array([-48, 48])      # Along-track limits to search for geolocation error
offsetsRasterResolutions = np.array([8, 4, 2, 1])  # Multi-resolutional step-down raster resolutions (in meters)
offsetsUseVerticalShift = False    # Option to use a vertical shift
offsetsVerticalShift = 0           # Vertical shift to use if above set to True (in meters)
measClassFilter = 1                # Meas Classes (0 = Unclass, 1 = Ground, 2 = Low Veg, 3 = High Veg)
truthClassFilter = 2               # Texpert Truth Classes (0 = Unclass, 2 = Ground, 4 = Veg, 6 = Building)
offsets = offsetsStruct(offsetsCrossTrackBounds, offsetsAlongTrackBounds, offsetsRasterResolutions, offsetsUseVerticalShift, offsetsVerticalShift )
createMeasCorrFile = True     # Option to create ouput measured corrected .las file
makePlots = True              # Option to make output plots
showPlots = False             # Option to show output plot windows
  
##### End Inputs for getMeasurementError  


##### CODE BELOW -- DO NOT EDIT ###############################################

timeStart = runTime.time()

atlMeasuredDataAll = []
atlTruthDataAll = []
atlCorrectionsAll = []

for i in range(0,len(gtNums)):
    
    gtNum = gtNums[i]

    # Call getAtlMeasuredSwath
    print('RUNNING getAtlMeasuredSwath...\n')
    atl03Data, atl08Data, headerData, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, gtNum, trimInfo, createAtl03LasFile, createAtl03KmlFile, createAtl08KmlFile, createAtl03CsvFile, createAtl08CsvFile)
        
    # Call getAtlTruthSwath
    print('RUNNING getAtlTruthSwath...\n')
    atlTruthData = getAtlTruthSwath(atl03Data, headerData, rotationData, useExistingTruth, truthSwathDir, buffer, outFilePath, createTruthFile)
                
    # Call getMeasurementError
    print('RUNNING getMeasurementError...\n')
    atlCorrections = getMeasurementError(atl03Data, atlTruthData, rotationData, outFilePath, measClassFilter, truthClassFilter, offsets, createMeasCorrFile, makePlots, showPlots)
    
    # Store all objects into one
    atlMeasuredDataAll.append(atl03Data)
    atlTruthDataAll.append(atlTruthData)
    atlCorrectionsAll.append(atlCorrections)
    
# EndFor

# End timer
timeEnd = runTime.time()
timeElapsedTotal = timeEnd - timeStart
timeElapsedMin = np.floor(timeElapsedTotal / 60)
timeElapsedSec = timeElapsedTotal % 60
    
# Print completion message
print('   Script Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
print('\n')