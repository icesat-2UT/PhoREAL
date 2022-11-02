# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 13:14:46 2019

This module will take in ATL03 MEASURED and TRUTH data and find the 
geolocation error in the MEASURED data based on TRUTH data.

@author: malonzo
"""

# Filter Runtime warnings
import warnings

# Import modules
import os
import numpy as np
import time as runTime
import sys
# from getAtlMeasuredSwath_auto import getAtlMeasuredSwath
# from getAtlTruthSwath_auto import getAtlTruthSwath
from phoreal.io import (writeLas, getTruthFilePaths, getTruthHeaders,
                      atlMeasuredDataReducedStruct, atlCorrectionsStruct, offsetsStruct)
from phoreal.utils import (ismember, getRaster, getIntersection2d, getCoordRotRev, getUTM2LatLon)
from phoreal.plot import (plotContour, plotZY, plotZT)

def limBounds(x, y, lim):
    # x = max([x, lim[0]])
    # y = min([lim[1], y])
    return x, y

# Function to get ICESat-2 offsets relative to reference data
# def getMeasurementError(atlMeasuredData, atlTruthData, refHeightType, 
#                         rotationData, outFilePath, 
#                         useMeasSigConf, filterData, offsets, createMeasCorrFile, 
#                         makePlots, showPlots, logFileID = False):

    
def getMeasurementError(atl03, truth_df):
    # Static variables
    refHeightType = 'hae'
    outFilePath = ''
    useMeasSigConf = False
    filterData = [1,2,3,4]
    createMeasCorrFile = False
    makePlots = False
    showPlots = False

    # static variable for offset object
    offsetsCrossTrackBounds = np.array([-25,25])      # Cross-track limits to search for geolocation error
    offsetsAlongTrackBounds = np.array([-25,25])      # Along-track limits to search for geolocation error
    offsetsRasterResolutions = np.array([8, 4, 2, 1])  # Multi-resolutional step-down raster resolutions (in meters)
    offsetsUseVerticalShift = False    # Option to use a vertical shift
    offsetsVerticalShift = 0           # Vertical shift to use if above set to True (in meters)
    offsets = offsetsStruct(offsetsCrossTrackBounds, offsetsAlongTrackBounds, 
                            offsetsRasterResolutions, offsetsUseVerticalShift, 
                            offsetsVerticalShift )
        
    rotationData = atl03.rotationData
    
    # Start timer
    timeStart = runTime.time()
    
    # Get and print ground track number
    gtNum = atl03.gtNum
    beamNum = atl03.beamNum
    beamStrength = atl03.beamStrength
    print('   Ground Track Number: %s (Beam #%s, Beam Strength: %s)\n' %(gtNum, beamNum, beamStrength))

    # Turn off runtime warnings (from np.nanmean)
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    # EndIf
    
    # Get array of offsets
    crossTrackBounds = offsets.crossTrackBounds
    alongTrackBounds = offsets.alongTrackBounds
    rasterResolutions = offsets.rasterResolutions
    
    # Get reduced MEASURED data (limited to bounds of TRUTH data)
    measInTruthInds = (np.array(atl03.df.alongtrack) >= np.array(truth_df.alongtrack).min()) & (np.array(atl03.df.alongtrack) <= np.array(truth_df.alongtrack).max())
    atlMeasuredDataReducedX = np.array(atl03.df.easting)[measInTruthInds]
    atlMeasuredDataReducedY = np.array(atl03.df.northing)[measInTruthInds]
    if('hae' in refHeightType.lower()):
        # Use ellipsoidal heights (HAE - Height Above Ellipsoid)
        print('   Comparing Reference Heights to ICESat-2 Ellipsoidal (HAE) Heights')
        atlMeasuredDataReducedZ = np.array(atl03.df.h_ph)[measInTruthInds]
    else:
        # User orthometric heights (MSL - Mean Sea Level)
        print('   Comparing Reference Heights to ICESat-2 Orthometric (MSL) Heights')
        atlMeasuredDataReducedZ = np.array(atl03.df.h_ph)[measInTruthInds]
    # endIf
    atlMeasuredDataReducedXrot = np.array(atl03.df.crosstrack)[measInTruthInds]
    atlMeasuredDataReducedYrot = np.array(atl03.df.alongtrack)[measInTruthInds]
    atlMeasuredDataReducedLat = np.array(atl03.df.lat_ph)[measInTruthInds]
    atlMeasuredDataReducedLon = np.array(atl03.df.lon_ph)[measInTruthInds]
    atlMeasuredDataReducedTime = np.array(atl03.df.time)[measInTruthInds]
    atlMeasuredDataReducedClass = np.array(atl03.df.classification)[measInTruthInds]
    atlMeasuredDataReducedSignal = np.array(atl03.df.signal_conf_ph)[measInTruthInds]
    
    # Get MEASURED reduced class
    atlMeasuredDataReduced = atlMeasuredDataReducedStruct(atlMeasuredDataReducedX, \
                                                          atlMeasuredDataReducedY, \
                                                          atlMeasuredDataReducedZ, \
                                                          atlMeasuredDataReducedXrot, \
                                                          atlMeasuredDataReducedYrot, \
                                                          atlMeasuredDataReducedLat, \
                                                          atlMeasuredDataReducedLon, \
                                                          atlMeasuredDataReducedTime, \
                                                          atlMeasuredDataReducedClass, \
                                                          atlMeasuredDataReducedSignal)
    
    # Filter Ground data from ATL08 if available
    if(useMeasSigConf):
    
        print('   Using Signal Confidence %s for ICESat-2 Offset Computation' %filterData)
        
        # Filter MEASURED data points based on signal confidence
        measGroundInds = ismember(atlMeasuredDataReduced.signalConf,filterData)[0]
        measXRot = atlMeasuredDataReduced.crossTrack[measGroundInds]
        measYRot = atlMeasuredDataReduced.alongTrack[measGroundInds]
        measZ = atlMeasuredDataReduced.z[measGroundInds]
        measT = atlMeasuredDataReduced.time[measGroundInds]
        
        # Do not filter ground TRUTH data points
        truthXRot = np.array(truth_df.crosstrack)
        truthYRot = np.array(truth_df.alongtrack)
        truthZ = np.array(truth_df.z)
    
    else:
    
        print('   Using Ground Truth for ICESat-2 Offset Computation')
        
        # Filter ground MEASURED data points
        icesatGroundIndex = 1
        measGroundInds = ismember(atlMeasuredDataReduced.classification,icesatGroundIndex)[0]
        measXRot = atlMeasuredDataReduced.crossTrack[measGroundInds]
        measYRot = atlMeasuredDataReduced.alongTrack[measGroundInds]
        measZ = atlMeasuredDataReduced.z[measGroundInds]
        measT = atlMeasuredDataReduced.time[measGroundInds]
        
        # Filter ground TRUTH data points
        truthGroundInds = ismember(np.array(truth_df.classification),filterData)[0]
        truthXRot = np.array(truth_df.crosstrack)[truthGroundInds]
        truthYRot = np.array(truth_df.alongtrack)[truthGroundInds]
        truthZ = np.array(truth_df.z)[truthGroundInds]
        
    # endIf
        
    # Check and make sure there is Measured and Truth data after filtering
    measDataExists = (len(measXRot)>0) and (len(measYRot)>0)
    truthDataExists = (len(truthXRot)>0) and (len(truthYRot)>0)
    if(measDataExists and truthDataExists):
    
        if((len(crossTrackBounds) > 1) and (len(alongTrackBounds) > 1)):
                
            # Initialize parameters
            correctionsCrossTrack = 0
            correctionsAlongTrack = 0
    
            crossTrackEdgeError = 0
            alongTrackEdgeError = 0
            numErr = 0
            crossTrackBoundLimits = [-48, 48]
            alongTrackBoundLimits = [-1000, 1000]
    
            # Get array offsets
            crossTrackMinBound = (np.fix(crossTrackBounds[0]/rasterResolutions[0])*rasterResolutions[0]).astype(int)
            crossTrackMaxBound = (np.fix((crossTrackBounds[-1])/rasterResolutions[0])*rasterResolutions[0]).astype(int)
            alongTrackMinBound = (np.fix(alongTrackBounds[0]/rasterResolutions[0])*rasterResolutions[0]).astype(int)
            alongTrackMaxBound = (np.fix((alongTrackBounds[-1])/rasterResolutions[0])*rasterResolutions[0]).astype(int)
            
            # Loop through all raster resolutions (multi-resolutional approach)
            totRasterResolutions = len(rasterResolutions)
            # for i in range(0,totRasterResolutions):
            i = 0
            while i < totRasterResolutions:
                
                # Get current raster resolution
                rasterResolution = rasterResolutions[i]
    
                # Set cross-track/along-track offsets
                if(i > 0):
                    
                    # Get previous raster resolution
                    rasterResolutionPrev = rasterResolutions[i - 1]
                    
                    # Set cross-track offsets
                    crossTrackBoundStart = correctionsCrossTrack - rasterResolutionPrev - rasterResolution + crossTrackEdgeError
                    crossTrackBoundEnd = correctionsCrossTrack + rasterResolutionPrev + rasterResolution + crossTrackEdgeError
                    crossTrackBoundStart, crossTrackBoundEnd = limBounds(crossTrackBoundStart, crossTrackBoundEnd, crossTrackBoundLimits)
                    crossTrackMinBound = (np.fix(crossTrackBoundStart/rasterResolution)*rasterResolution).astype(int)
                    crossTrackMaxBound = (np.fix(crossTrackBoundEnd/rasterResolution)*rasterResolution).astype(int)
                    crossTrackOffsets = np.arange(crossTrackMinBound, crossTrackMaxBound + 1, rasterResolution)
                
                    # Set along-track offsets
                    alongTrackBoundStart = correctionsAlongTrack - rasterResolutionPrev - rasterResolution + alongTrackEdgeError
                    alongTrackBoundEnd = correctionsAlongTrack + rasterResolutionPrev + rasterResolution + alongTrackEdgeError
                    alongTrackBoundStart, alongTrackBoundEnd = limBounds(alongTrackBoundStart, alongTrackBoundEnd, alongTrackBoundLimits)
                    alongTrackMinBound = (np.fix(alongTrackBoundStart/rasterResolution)*rasterResolution).astype(int)
                    alongTrackMaxBound = (np.fix(alongTrackBoundEnd/rasterResolution)*rasterResolution).astype(int)
                    alongTrackOffsets = np.arange(alongTrackMinBound, alongTrackMaxBound + 1, rasterResolution)
        
                else:
    
                    # Set cross-track offsets
                    crossTrackBoundStart = crossTrackBounds[0] + crossTrackEdgeError
                    crossTrackBoundEnd = crossTrackBounds[-1] + crossTrackEdgeError
                    crossTrackBoundStart, crossTrackBoundEnd = limBounds(crossTrackBoundStart, crossTrackBoundEnd, crossTrackBoundLimits)
                    crossTrackMinBound = (np.fix(crossTrackBoundStart/rasterResolutions[0])*rasterResolutions[0]).astype(int)
                    crossTrackMaxBound = (np.fix(crossTrackBoundEnd/rasterResolutions[0])*rasterResolutions[0]).astype(int)
                    crossTrackOffsets = np.arange(crossTrackMinBound, crossTrackMaxBound + 1, rasterResolutions[0])
    
                    # Set along-track offsets
                    alongTrackBoundStart = alongTrackBounds[0] + alongTrackEdgeError
                    alongTrackBoundEnd = alongTrackBounds[-1] + alongTrackEdgeError
                    alongTrackBoundStart, alongTrackBoundEnd = limBounds(alongTrackBoundStart, alongTrackBoundEnd, alongTrackBoundLimits)
                    alongTrackMinBound = (np.fix(alongTrackBoundStart/rasterResolutions[0])*rasterResolutions[0]).astype(int)
                    alongTrackMaxBound = (np.fix(alongTrackBoundEnd/rasterResolutions[0])*rasterResolutions[0]).astype(int)
                    alongTrackOffsets = np.arange(alongTrackMinBound, alongTrackMaxBound + 1, rasterResolutions[0])
    
                # endIf
                
                # Print message to screen
                print('')
                print('   Test Case %d of %d' %(i+1, totRasterResolutions))
                if not (crossTrackEdgeError == 0 and alongTrackEdgeError == 0):
                    print('   Relocated search center to [%s, %s]' % (crossTrackEdgeError, alongTrackEdgeError))
    
                print('   Raster Resolution = %s m' % rasterResolution)
                print('   Cross Track Offsets = [%s, %s] (Adjusted to align with %s m resolution)' %(crossTrackOffsets[0], crossTrackOffsets[-1], rasterResolution))
                print('   Along Track Offsets = [%s, %s] (Adjusted to align with %s m resolution)' %(alongTrackOffsets[0], alongTrackOffsets[-1], rasterResolution))
                
                # Make Offsets Int's and flip along-track offsets array so larger numbers at top
                crossTrackOffsets = crossTrackOffsets.astype(int)
                alongTrackOffsets = (np.flipud(alongTrackOffsets)).astype(int)
                
                # Format TRUTH data
                truthXRotReduced = np.ravel(truthXRot)
                truthYRotReduced = np.ravel(truthYRot)
                truthZReduced = np.ravel(truthZ)
                
                # Rasterize MEASURED data
                gridMethod = 'Mean'
                fillValue = np.nan
                print('   Gridding ICESat-2 Data at %s m Resolution using %s Values...' % (rasterResolution, gridMethod))
                measRasterRot = getRaster(measXRot, measYRot, measZ, rasterResolution, gridMethod, fillValue, measT)
                
                # Rasterize TRUTH data
                print('   Gridding Reference Data at %s m Resolution using %s Values...' % (rasterResolution, gridMethod))
                truthRasterRot = getRaster(truthXRotReduced, truthYRotReduced, truthZReduced, rasterResolution, gridMethod, fillValue)
                
                # Find offsets with minimum MAE
                print('   Finding Offsets with Minimum MAE...')
                
                # Store MEASURED data into array variables
                measRasterXRot = np.c_[np.ravel(measRasterRot.x)]
                measRasterYRot = np.c_[np.ravel(measRasterRot.y)]
                measRasterZ    = np.c_[np.ravel(measRasterRot.grid)]
                # measRasterT    = np.c_[np.ravel(measRasterRot.t)]
                
                # Store TRUTH data into array variables
                truthRasterXRot = np.c_[np.ravel(truthRasterRot.x)]
                truthRasterYRot = np.c_[np.ravel(truthRasterRot.y)]
                truthRasterZ    = np.c_[np.ravel(truthRasterRot.grid)]
    
                # Create 2 column array of all MEASURED X,Y data
                measXYall = np.column_stack([measRasterXRot, measRasterYRot])
                
                # Create 2 column array of all TRUTH X,Y data
                truthXYall = np.column_stack([truthRasterXRot, truthRasterYRot])
    
                # Find common TRUTH and MEASURED indices
                _, truthIndsCommon, measIndsCommon = getIntersection2d(truthXYall, measXYall, assume_unique=True)
                
                # Only use MEASURED indices in common with TRUTH indices            
                measRasterXRot_common = measRasterXRot[measIndsCommon]
                measRasterYRot_common = measRasterYRot[measIndsCommon]
                measRasterZ_common = measRasterZ[measIndsCommon]
                # measRasterT_common = measRasterT[measIndsCommon]
                
                # Find row, col for indices in TRUTH raster (that are in common with MEASURED raster)
                truthColsInit = (truthIndsCommon%np.shape(truthRasterRot.x)[1]).astype(int)
                truthRowsInit = (np.floor(truthIndsCommon/np.shape(truthRasterRot.x)[1])).astype(int)
                
                # Find min/max row, col for indicies in TRUTH raster
                truthMinXBounds = 0
                truthMaxXBounds = np.shape(truthRasterRot.grid)[1]
                truthMinYBounds = 0
                truthMaxYBounds = np.shape(truthRasterRot.grid)[0]
                        
                # Pre-allocate results data
                numHorzCombos = np.size(crossTrackOffsets)
                numVertCombos = np.size(alongTrackOffsets)
                resultsCrossTrackShift = np.zeros((numVertCombos,numHorzCombos))
                resultsAlongTrackShift = np.zeros((numVertCombos,numHorzCombos))
                resultsVerticalShift = np.zeros((numVertCombos,numHorzCombos))
                resultsMAE = np.zeros((numVertCombos,numHorzCombos))
                resultsRMSE = np.zeros((numVertCombos,numHorzCombos))
                resultsME = np.zeros((numVertCombos,numHorzCombos))
                
                # Loop through all combos and find min MAE
                for I in range(0,numHorzCombos):
                    for J in range(0,numVertCombos):
                        
                        # Get cross-track/along-track offset values
                        crossTrackShift = crossTrackOffsets[I]
                        alongTrackShift = alongTrackOffsets[J]
                    
                        # Get amount to shift X,Y indices
                        moveIndsX = (crossTrackShift/rasterResolution).astype(int)
                        moveIndsY = (alongTrackShift/rasterResolution).astype(int)
                        
                        # Get new TRUTH data from shifted X,Y indices
                        truthColsCurrent = truthColsInit + moveIndsX
                        truthRowsCurrent = truthRowsInit - moveIndsY # y array is inverted
    
                        # Determine indices in bounds
                        colsToKeep = (truthColsCurrent>=truthMinXBounds) & (truthColsCurrent<truthMaxXBounds)
                        rowsToKeep = (truthRowsCurrent>=truthMinYBounds) & (truthRowsCurrent<truthMaxYBounds)
                        indsToKeep = colsToKeep & rowsToKeep
                        
                        # Get current TRUTH Z points
                        truthColsCurrent = truthColsCurrent[indsToKeep]
                        truthRowsCurrent = truthRowsCurrent[indsToKeep]
                        truthRasterZCurr = np.c_[np.ravel(truthRasterRot.grid[truthRowsCurrent, truthColsCurrent])]
                        
                        # Get current MEASURED Z points
                        measRasterZCurr = measRasterZ_common[indsToKeep]
                        
                        # Find Z error
                        zError = truthRasterZCurr - measRasterZCurr
                        
                        # Determine vertical shift
                        if(offsets.useVerticalShift == 1):
                                
                            # Use defined vertical shift
                            verticalShift = offsets.verticalShift
                                
                        else:
                                
                            # MAE is tied to median error (RMSE to mean error)
                            verticalShift = np.nanmedian(zError)
                                
                        # Endif
                        
                        if(np.logical_not(np.isnan(verticalShift))):
                        
                            # Apply vertical shift
                            measRasterZShifted = measRasterZCurr + verticalShift
                            
                            # Get new Z error
                            zError = truthRasterZCurr - measRasterZShifted
                            
                            # Get Mean Absolute Error and RMSE
                            MAE = np.nanmean(abs(zError))
                            RMSE = np.sqrt(np.nanmean(zError**2))
                            ME = np.nanmean(zError)
                        
                        else:
                        
                            # Populate with default values
                            verticalShift = np.nan
                            MAE = np.nan
                            RMSE = np.nan
                            ME = np.nan
                        
                        # endIf
                        
                        # Store results
                        resultsCrossTrackShift[J,I] = crossTrackShift
                        resultsAlongTrackShift[J,I] = alongTrackShift
                        resultsVerticalShift[J,I] = verticalShift
                        resultsMAE[J,I] = MAE
                        resultsRMSE[J,I] = RMSE
                        resultsME[J,I] = ME
                    
                    # endFor
                    
                # endFor

                # Find row, col for min MAE
                minRow = np.where(resultsMAE == np.nanmin(resultsMAE))[0]
                minCol = np.where(resultsMAE == np.nanmin(resultsMAE))[1]
                if( (len(minRow)==0) or (len(minCol)==0) ):
                    print('Error: Gridding all NaN, ensure ICESat-2/Reference files are correct')
                    return []
    
                if(len(minRow)>1):
                    minRow = np.array([np.median(minRow)], dtype='int64')
                # endIf
                
                if(len(minCol)>1):
                    minCol = np.array([np.median(minCol)], dtype='int64')
                # endIf
                
                # Check that estimate is not close to an edge
                nr, nc = np.shape(resultsMAE)
                dr, dc = 1, 1 # edge "closeness", index units
                inBounds = minRow in range(0+dr, nr-dr) and minCol in range(0+dc, nc-dc)
                if not inBounds:
                    # estimated an edge
    
                    p = np.array([float(resultsCrossTrackShift[minRow,minCol]), float(resultsAlongTrackShift[minRow,minCol])])
                    mid = np.array([np.mean([crossTrackBoundStart, crossTrackBoundEnd]), np.mean([alongTrackBoundStart, alongTrackBoundEnd])])
                    vec_diff = p - mid
    
                    numErr += 1
                    if numErr == 5:
                        print('   Error: Could not determine valid estimate, consider increasing search bounds')
                        return []
    
                    else:
                        crossTrackEdgeError += vec_diff[0] #p[0]
                        alongTrackEdgeError += vec_diff[1] #p[1]
                        print('   Warning: Estimate on edge [%4.1f, %4.1f] m, relocating search and re-running test case %d...' % (p[0], p[1], i+1))
                        continue # do not update i
    
                # Set edge error to zero if inBounds, this way
                # we center on estimate and zoom in since estimate is good
                numErr = 0
                crossTrackEdgeError = 0
                alongTrackEdgeError = 0
    
                # Get output associated with min MAE            
                correctionsCrossTrack = resultsCrossTrackShift[minRow, minCol]
                correctionsAlongTrack = resultsAlongTrackShift[minRow, minCol]
                correctionsVertical = resultsVerticalShift[minRow, minCol]
                correctionsMAE = resultsMAE[minRow, minCol]
                correctionsRMSE = resultsRMSE[minRow, minCol]
                correctionsME = resultsME[minRow, minCol]
            
                # Get CT/AT offsets in Easting/Northing plane
                R_mat = rotationData.R_mat
                xRotPt = rotationData.xRotPt
                yRotPt = rotationData.yRotPt
                correctionsEasting, correctionsNorthing, _, _, _ = getCoordRotRev(correctionsCrossTrack, correctionsAlongTrack, R_mat, xRotPt, yRotPt)
                correctionsEasting -= xRotPt
                correctionsNorthing -= yRotPt
                
                # Print results to screen
                print('   Cross-Track Correction = %d m (Easting = %0.1f m)'  % (correctionsCrossTrack, correctionsEasting))
                print('   Along-Track Correction = %d m (Northing = %0.1f m)' % (correctionsAlongTrack, correctionsNorthing))
                print('   Vertical Correction = %0.1f m' % correctionsVertical)
                print('   MAE = %0.2f m' % correctionsMAE)
                print('   RMSE = %0.2f m' % correctionsRMSE)
                
                # Get final MEASURED raster data
                measRasterXRotFinal= measRasterXRot + correctionsCrossTrack
                measRasterYRotFinal = measRasterYRot + correctionsAlongTrack
                measRasterZFinal = measRasterZ + correctionsVertical
                # measRasterTFinal = measRasterT
        
                # Get amount to shift X,Y indices
                moveIndsX = (correctionsCrossTrack/rasterResolution).astype(int)
                moveIndsY = (correctionsAlongTrack/rasterResolution).astype(int)
                
                # Get new TRUTH data from shifted X,Y indices
                truthColsCurrent = truthColsInit + moveIndsX
                truthRowsCurrent = truthRowsInit - moveIndsY # y array is inverted
    
                # Determine indices in bounds
                colsToKeep = (truthColsCurrent>=truthMinXBounds) & (truthColsCurrent<truthMaxXBounds)
                rowsToKeep = (truthRowsCurrent>=truthMinYBounds) & (truthRowsCurrent<truthMaxYBounds)
                indsToKeep = colsToKeep & rowsToKeep
                
                # Get final common TRUTH data
                truthColsCurrent = truthColsCurrent[indsToKeep]
                truthRowsCurrent = truthRowsCurrent[indsToKeep]
                truthRasterXRotCommonFinal = np.c_[np.ravel(truthRasterRot.x[truthRowsCurrent, truthColsCurrent])]
                truthRasterYRotCommonFinal = np.c_[np.ravel(truthRasterRot.y[truthRowsCurrent, truthColsCurrent])]
                truthRasterZCommonFinal = np.c_[np.ravel(truthRasterRot.grid[truthRowsCurrent, truthColsCurrent])]
                        
                # Get final common MEASURED data
                measRasterXRotCommonFinal = measRasterXRot_common[indsToKeep] + correctionsCrossTrack
                measRasterYRotCommonFinal = measRasterYRot_common[indsToKeep] + correctionsAlongTrack
                measRasterZCommonFinal = measRasterZ_common[indsToKeep] + correctionsVertical
                # measRasterTCommonFinal = measRasterT_common[indsToKeep]
                
                # Get Z Error raster (MEASURED - TRUTH)
                zErrorCommonFinal = measRasterZCommonFinal - truthRasterZCommonFinal
                
                # Find Z values that are not NaN in both MEASURED and TRUTH data
                measNotNan = np.logical_not(np.isnan(measRasterZCommonFinal))
                truthNotNan = np.logical_not(np.isnan(truthRasterZCommonFinal))
                allNotNan = measNotNan & truthNotNan
                
                # Filter out non-NaN values in TRUTH data
                truthRasterXRotCommonFinal = truthRasterXRotCommonFinal[allNotNan]
                truthRasterYRotCommonFinal = truthRasterYRotCommonFinal[allNotNan]
                truthRasterZCommonFinal = truthRasterZCommonFinal[allNotNan]
                
                # Filter out non-NaN values in MEASURED data
                measRasterXRotCommonFinal = measRasterXRotCommonFinal[allNotNan]
                measRasterYRotCommonFinal = measRasterYRotCommonFinal[allNotNan]
                measRasterZCommonFinal = measRasterZCommonFinal[allNotNan]
                # measRasterTCommonFinal = measRasterTCommonFinal[allNotNan]
                  
                # Filter out non-NaN values in error data
                zErrorCommonFinal = zErrorCommonFinal[allNotNan]
                
                # Get 100 m segment Z Mean Error data
                errorResolution = 100
                gridMethod = 'Mean'
                fillValue = np.nan
                try:
                    segmentError = getRaster(np.ravel(measRasterXRotCommonFinal), np.ravel(measRasterYRotCommonFinal), np.ravel(zErrorCommonFinal), errorResolution, gridMethod, fillValue)
                    segmentErrorX = np.mean(segmentError.x, axis = 1)
                    segmentErrorY = np.mean(segmentError.y, axis = 1)
                    segmentErrorZ = np.nanmean(segmentError.grid, axis = 1)
                except:
                    segmentErrorX = np.NaN
                    segmentErrorY = np.NaN
                    segmentErrorZ = np.NaN
                # endIf 
                
                # Format Y,Z values for plotting
                segmentErrorXYZ = np.column_stack([segmentErrorX, segmentErrorY, segmentErrorZ])
                segmentErrorXYZsorted = segmentErrorXYZ[segmentErrorXYZ[:,1].argsort(),]
                segmentErrorLo = segmentErrorXYZsorted[:,1] - errorResolution*0.5
                segmentErrorHi = segmentErrorXYZsorted[:,1] + errorResolution*0.5
                segmentErrorXin = (np.column_stack([segmentErrorXYZsorted[:,0], segmentErrorXYZsorted[:,0]])).flatten()
                segmentErrorYin = (np.column_stack([segmentErrorLo, segmentErrorHi])).flatten()
                segmentErrorZin = (np.column_stack([segmentErrorXYZsorted[:,2], segmentErrorXYZsorted[:,2]])).flatten()
        
                # Convert to easting/northing plane
                _, segmentErrorY_plot, _, _, _ = getCoordRotRev(segmentErrorXin, segmentErrorYin, R_mat, xRotPt, yRotPt)        
                segmentErrorZ_plot = segmentErrorZin
        
                # Convert common measured and truth data back to Easting/Northing plane
                try:
                    measRasterXCommonFinal, measRasterYCommonFinal, _, _, _ = getCoordRotRev(measRasterXRotCommonFinal, measRasterYRotCommonFinal, R_mat, xRotPt, yRotPt)
                    truthRasterXCommonFinal, truthRasterYCommonFinal, _, _, _ = getCoordRotRev(truthRasterXRotCommonFinal, truthRasterYRotCommonFinal, R_mat, xRotPt, yRotPt)
                except:
                    measRasterXCommonFinal = np.array([])
                    measRasterYCommonFinal = np.array([])
                    truthRasterXCommonFinal = np.array([])
                    truthRasterYCommonFinal = np.array([])
                # endIf 
                
                # Get shifted data
                shiftedEasting = atl03.df.easting + correctionsEasting
                shiftedNorthing = atl03.df.northing + correctionsNorthing
                shiftedCrossTrack = atl03.df.crosstrack + correctionsCrossTrack
                shiftedAlongTrack = atl03.df.alongtrack + correctionsAlongTrack
                shiftedVertical = atl03.df.h_ph + correctionsVertical
                # shiftedVerticalMsl = atlMeasuredData.zMsl + correctionsVertical
                
                # Get lat/lon
                # shiftedLat, shiftedLon = getUTM2LatLon(shiftedEasting, 
                #                                        shiftedNorthing, 
                #                                        atlMeasuredData.zone, 
                #                                        atlMeasuredData.hemi)
            
                # Store data in class structure
                atlCorrections = atlCorrectionsStruct(shiftedEasting, shiftedNorthing, \
                                                      shiftedCrossTrack, shiftedAlongTrack, \
                                                      [], [], \
                                                      [], [], \
                                                      [], [], \
                                                      [], \
                                                      [], \
                                                      [], [], \
                                                      correctionsCrossTrack, correctionsAlongTrack, \
                                                      correctionsEasting, correctionsNorthing, \
                                                      correctionsVertical, \
                                                      correctionsMAE, correctionsRMSE, correctionsME, \
                                                      measRasterXCommonFinal, \
                                                      measRasterYCommonFinal, \
                                                      measRasterZCommonFinal, \
                                                      truthRasterXCommonFinal, \
                                                      truthRasterYCommonFinal, \
                                                      truthRasterZCommonFinal)
                # Make plots
                # if(makePlots):
                    
                #     # Make contour plot
                #     print('   Making Plots...')
                #     plotContour(resultsCrossTrackShift, resultsAlongTrackShift, resultsMAE, 
                #                 atlMeasuredData, atlCorrections, rasterResolution, showPlots, 
                #                 outFilePath, i)
                
                #     # Make ZY scatter plots
                #     plotZY(measRasterYCommonFinal, measRasterZCommonFinal, truthRasterYCommonFinal, 
                #            truthRasterZCommonFinal, zErrorCommonFinal, segmentErrorY_plot, 
                #            segmentErrorZ_plot, atlMeasuredData, atlCorrections, rasterResolution, 
                #            useMeasSigConf, filterData, showPlots, outFilePath, i)
                    
                #     # Make ZT scatter plots
                #     plotZT(measRasterTCommonFinal, measRasterZCommonFinal, zErrorCommonFinal, 
                #            atlMeasuredData, atlCorrections, rasterResolution, 
                #            useMeasSigConf, filterData, showPlots, outFilePath, i)
        
                # endIf 
    
                i += 1
                
            # endFor
            
        else:
                
            # Hard-coded offsets case
            print('   Using preset values...')
            
            # Set cross-track/along-track offsets
            correctionsCrossTrack = crossTrackBounds
            correctionsAlongTrack = alongTrackBounds
            
            # Get CT/AT offsets in Easting/Northing plane
            R_mat = rotationData.R_mat
            xRotPt = rotationData.xRotPt
            yRotPt = rotationData.yRotPt
            correctionsEasting, correctionsNorthing, _, _, _ = getCoordRotRev(correctionsCrossTrack, correctionsAlongTrack, R_mat, xRotPt, yRotPt)
            correctionsEasting -= xRotPt
            correctionsNorthing -= yRotPt
            
            # Get vertical offset
            if(offsets.useVerticalShift):
                correctionsVertical = np.array([float(offsets.verticalShift)])
            else:
                correctionsVertical = np.array([0.0])
            # endIf
            
            # Rasterize MEASURED data
            rasterResolution = rasterResolutions[-1]
            gridMethod = 'Mean'
            fillValue = -999
            print('   Gridding ICESat-2 Data at %s m Resolution using %s Values...' % (rasterResolution, gridMethod))
            
            # Add in corrections
            measXRot += correctionsCrossTrack
            measYRot += correctionsAlongTrack
            measZ += correctionsVertical
            
            # Raster MEASURED data
            measRasterRot = getRaster(measXRot, measYRot, measZ, rasterResolution, gridMethod, fillValue, measT)
            
            # Reduce TRUTH data to width of cross-track offsets
    #        measWidth = (np.ceil( ( measXRot.max() - measXRot.min() ) / 2.0 )).astype(int)
    #        minTruthBound = (crossTrackBounds - measWidth).astype(int) - 5
    #        maxTruthBound = (crossTrackBounds + measWidth).astype(int) + 5
    #        truthInOffsetsInds = (truthXRot >= minTruthBound) & (truthXRot <= maxTruthBound)
    #        truthXRotReduced = truthXRot[truthInOffsetsInds]
    #        truthYRotReduced = truthYRot[truthInOffsetsInds]
    #        truthZReduced = truthZ[truthInOffsetsInds]
            truthXRotReduced = truthXRot
            truthYRotReduced = truthYRot
            truthZReduced = truthZ
                
            # Rasterize TRUTH data
            print('   Gridding Reference Data at %s m Resolution using %s Values...' % (rasterResolution, gridMethod))
            truthRasterRot = getRaster(truthXRotReduced, truthYRotReduced, truthZReduced, rasterResolution, gridMethod, fillValue)
            
            # Remove NaN data (-999) in MEASURED raster
            measIndsToKeep = measRasterRot.grid != -999
            measRasterXRot = measRasterRot.x[measIndsToKeep]
            measRasterYRot = measRasterRot.y[measIndsToKeep]
            measRasterZ    = measRasterRot.grid[measIndsToKeep]
            # measRasterT    = measRasterRot.t[measIndsToKeep]
            
            # Get final MEASURED raster data
            measRasterXRotFinal= measRasterXRot
            measRasterYRotFinal = measRasterYRot
            measRasterZFinal = measRasterZ
            # measRasterTFinal = measRasterT
                
            # Remove NaN data (-999) in TRUTH raster
            truthIndsToKeep = truthRasterRot.grid != -999
            truthRasterXRot = truthRasterRot.x[truthIndsToKeep]
            truthRasterYRot = truthRasterRot.y[truthIndsToKeep]
            truthRasterZ    = truthRasterRot.grid[truthIndsToKeep]
                
            # Get common MEASURED and TRUTH raster indices 
            xyMeasRasterRotFinal = np.concatenate((np.c_[measRasterXRotFinal],np.c_[measRasterYRotFinal]),axis = 1)
            xyTruthRasterRot = np.concatenate((np.c_[truthRasterXRot],np.c_[truthRasterYRot]),axis = 1)
            commonValsFinal, xyCommonMeasIndsFinal, xyCommonTruthIndsFinal = getIntersection2d(xyMeasRasterRotFinal, xyTruthRasterRot, assume_unique=True)
            
            # Get common MEASURED raster points
            measRasterXRotCommonFinal = measRasterXRotFinal[xyCommonMeasIndsFinal]
            measRasterYRotCommonFinal = measRasterYRotFinal[xyCommonMeasIndsFinal]
            measRasterZCommonFinal = measRasterZFinal[xyCommonMeasIndsFinal]
            # measRasterTCommonFinal = measRasterTFinal[xyCommonMeasIndsFinal]
            
            # Get common TRUTH raster points
            truthRasterXRotCommonFinal = truthRasterXRot[xyCommonTruthIndsFinal]
            truthRasterYRotCommonFinal = truthRasterYRot[xyCommonTruthIndsFinal]
            truthRasterZCommonFinal = truthRasterZ[xyCommonTruthIndsFinal]
            
            # Get Z Error raster (TRUTH - MEASURED)
            zErrorCommonFinal = measRasterZCommonFinal - truthRasterZCommonFinal
            
            # Get 100 m segment Z Mean Error data
            errorResolution = 100
            gridMethod = 'Mean'
            fillValue = np.nan
            try:
                segmentError = getRaster(measRasterXRotCommonFinal, measRasterYRotCommonFinal, zErrorCommonFinal, errorResolution, gridMethod, fillValue)
                segmentErrorX = np.mean(segmentError.x, axis = 1)
                segmentErrorY = np.mean(segmentError.y, axis = 1)
                segmentErrorZ = np.nanmean(segmentError.grid, axis = 1)
            except:
                segmentErrorX = np.NaN
                segmentErrorY = np.NaN
                segmentErrorZ = np.NaN
            # endIf 
            
            # Format Y,Z values for plotting
            segmentErrorXYZ = np.column_stack([segmentErrorX, segmentErrorY, segmentErrorZ])
            segmentErrorXYZsorted = segmentErrorXYZ[segmentErrorXYZ[:,1].argsort(),]
            segmentErrorLo = segmentErrorXYZsorted[:,1] - errorResolution*0.5
            segmentErrorHi = segmentErrorXYZsorted[:,1] + errorResolution*0.5
            segmentErrorXin = (np.column_stack([segmentErrorXYZsorted[:,0], segmentErrorXYZsorted[:,0]])).flatten()
            segmentErrorYin = (np.column_stack([segmentErrorLo, segmentErrorHi])).flatten()
            # segmentErrorZin = (np.column_stack([segmentErrorXYZsorted[:,2], segmentErrorXYZsorted[:,2]])).flatten()
    
            # Convert to easting/northing plane
            _, segmentErrorY_plot, _, _, _ = getCoordRotRev(segmentErrorXin, segmentErrorYin, R_mat, xRotPt, yRotPt)        
            # segmentErrorZ_plot = segmentErrorZin
    
            # Convert common measured and truth data back to Easting/Northing plane
            try:
                measRasterXCommonFinal, measRasterYCommonFinal, _, _, _ = getCoordRotRev(measRasterXRotCommonFinal, measRasterYRotCommonFinal, R_mat, xRotPt, yRotPt)
                truthRasterXCommonFinal, truthRasterYCommonFinal, _, _, _ = getCoordRotRev(truthRasterXRotCommonFinal, truthRasterYRotCommonFinal, R_mat, xRotPt, yRotPt)
            except:
                measRasterXCommonFinal = np.array([])
                measRasterYCommonFinal = np.array([])
                truthRasterXCommonFinal = np.array([])
                truthRasterYCommonFinal = np.array([])
            # endIf 
    
            MAE = np.nanmean(abs(zErrorCommonFinal))
            RMSE = np.sqrt(np.nanmean(zErrorCommonFinal**2))
            ME = np.nanmean(zErrorCommonFinal)
                            
            # Set error values
            correctionsMAE = np.array([MAE])
            correctionsRMSE = np.array([RMSE])
            correctionsME = np.array([ME])
        
            # Get shifted data
            shiftedEasting = np.array(atl03.df.easting) + correctionsEasting
            shiftedNorthing = np.array(atl03.df.northing) + correctionsNorthing
            shiftedCrossTrack = np.array(atl03.df.crosstrack) + correctionsCrossTrack
            shiftedAlongTrack = np.array(atl03.df.alongtrack) + correctionsAlongTrack
            shiftedVertical = np.array(atl03.df.h_ph) + correctionsVertical
            # shiftedVerticalMsl = atlMeasuredData.zMsl + correctionsVertical
            
            # Get lat/lon
            # shiftedLat, shiftedLon = getUTM2LatLon(shiftedEasting, 
            #                                        shiftedNorthing, 
            #                                        atlMeasuredData.zone, 
            #                                        atlMeasuredData.hemi)
            
            # Generate atlCorrections struct
            atlCorrections = atlCorrectionsStruct(shiftedEasting, shiftedNorthing, \
                                                  shiftedCrossTrack, shiftedAlongTrack, \
                                                  [], [], \
                                                  shiftedVertical, [], \
                                                  np.array(atl03.df.time), [], \
                                                  np.array(atl03.df.classification), \
                                                  np.array(atl03.df.signal_conf_ph), \
                                                  [], [], \
                                                  correctionsCrossTrack, correctionsAlongTrack, \
                                                  correctionsEasting, correctionsNorthing, \
                                                  correctionsVertical, \
                                                  correctionsMAE, correctionsRMSE, correctionsME, \
                                                  measRasterXCommonFinal, \
                                                  measRasterYCommonFinal, \
                                                  measRasterZCommonFinal, \
                                                  truthRasterXCommonFinal, \
                                                  truthRasterYCommonFinal, \
                                                  truthRasterZCommonFinal)

            # atlCorrections = atlCorrectionsStruct(shiftedEasting, shiftedNorthing, \
            #                                       shiftedCrossTrack, shiftedAlongTrack, \
            #                                       shiftedLat, shiftedLon, \
            #                                       shiftedVertical, shiftedVerticalMsl, \
            #                                       np.array(atl03.df.time), atlMeasuredData.deltaTime, \
            #                                       np.array(atl03.df.classification), \
            #                                       np.array(atl03.df.signal_conf_ph), \
            #                                       atlMeasuredData.zone, atlMeasuredData.hemi, \
            #                                       correctionsCrossTrack, correctionsAlongTrack, \
            #                                       correctionsEasting, correctionsNorthing, \
            #                                       correctionsVertical, \
            #                                       correctionsMAE, correctionsRMSE, correctionsME, \
            #                                       measRasterXCommonFinal, \
            #                                       measRasterYCommonFinal, \
            #                                       measRasterZCommonFinal, \
            #                                       truthRasterXCommonFinal, \
            #                                       truthRasterYCommonFinal, \
            #                                       truthRasterZCommonFinal)
                            
            # Make plots
            # if(makePlots):
                
            #     # Make contour plot
            #     print('   Making Plots...')
            
            #     # Plot counter
            #     i = 0
            
            #     # Make ZY scatter plots
            #     plotZY(measRasterYCommonFinal, measRasterZCommonFinal, truthRasterYCommonFinal,
            #            truthRasterZCommonFinal, zErrorCommonFinal, segmentErrorY_plot, 
            #            segmentErrorZ_plot, atlMeasuredData, atlCorrections, rasterResolution, 
            #            useMeasSigConf, filterData, showPlots, outFilePath, i)
                
            #     # Make ZT scatter plots
            #     plotZT(measRasterTCommonFinal, measRasterZCommonFinal, zErrorCommonFinal, 
            #            atlMeasuredData, atlCorrections, rasterResolution, 
            #            useMeasSigConf, filterData, showPlots, outFilePath, i)
    
            # endIf 
                
        # endIf
    
    else:
        
        # Either no measured data or truth data after filtering, so send empty data out
        # atlCorrections = []
        atlCorrections = atlCorrectionsStruct([], [], \
                                                              [], [], \
                                                              [], [], \
                                                              [], [], \
                                                              [], [], \
                                                              [], \
                                                              [], \
                                                              [], [], \
                                                              [], [], \
                                                              [], [], \
                                                              [], \
                                                              [], [], [], \
                                                              [], \
                                                              [], \
                                                              [], \
                                                              [], \
                                                              [], \
                                                              [])
        
        if(not measDataExists):
            print('   WARNING: No ICESat-2 data in selected filter options: %s' %filterData)
        # endIf
        
        if(not truthDataExists):
            print('   WARNING: No Reference data in selected filter options: %s' %filterData)
        # endIf
        
    # endIf
    
    # End timer
    timeEnd = runTime.time()
    timeElapsedTotal = timeEnd - timeStart
    timeElapsedMin = np.floor(timeElapsedTotal / 60)
    timeElapsedSec = timeElapsedTotal % 60
    
    # Print completion message
    print('')
    print('   Module Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
    print('\n')
    
    # Return output data
    return atlCorrections


# Unit test
if __name__ == "__main__":
    
    print('GET GEOLOCATION ERROR (UNIT TEST):\n')
    
#     ##### Start Inputs for getAtlMeasuredSwath

#     # Path to ATL03 Input File
#     # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20181126114738_08990103_001_01.h5' # FINLAND
#     # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20181030110205_04860106_001_01.h5' # SONONMA
#     # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20190101195003_00670202_001_01_sreq_2696.h5' # SONOMA
#     # atl03FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL03_20190228170214_09510202_001_02_sreq_2696.h5' # SONOMA
#     # atl03FilePath = '//bigtex/laserpewpew/data/release/001/ATL03_r001/ATL03_20190426213703_04370308_001_01.h5' # Brazil    
    
#     # Path to ATL08 Input File
#     # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20181126114738_08990103_952_01.h5' # FINLAND
#     # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20181030110205_04860106_952_01.h5' # SONOMA
#     # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20190101195003_00670202_952_01.h5' # SONOMA
#     # atl08FilePath = '//lidar-server/lidar/USERS/eric/benjelly_atl08/ATL08_20190228170214_09510202_952_02.h5' # SONOMA
#     # atl08FilePath = False    
    
#     # Path to Output Directory
#     # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_finland_20181126_python' # FINLAND
#     # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_sonoma_20181030_python' # SONOMA
#     # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_sonoma_20190101_python' # SONOMA
#     # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/r001_sonoma_20190228_python' # SONOMA
#     # outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test' # TEST
    
# #    atl03FilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/ATL03_20181030110205_04860106_001_01_sub_218.h5' # WSMR
# #    atl08FilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/ATL08_20181030110205_04860106_001_01_sub_218.h5' # WSMR
# #    outFilePath = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test'
    
#     atl03FilePath = 'Z:/data/release/002/ATL03_r002/ATL03_20190928175636_00280506_002_01_sreq_3181.h5' # WSMR
#     atl08FilePath = []
#     outFilePath = 'C:/Users/malonzo/GLAM/delete'
    
# #    atl03FilePath = 'z:/data/release/002/ATL03_r002/ATL03_20181018202814_03090102_002_01.h5' # N Carolina 20181018 r002
# #    atl08FilePath = []
# #    outFilePath = 'N:/USERS/mike/iceSat2/atl03_validation/r002_ncarolina_20181018_test'
    
#     # Ground track number to analyze
#     gtNum = 'gt2r'
    
#     # User options
#     trimInfo = 'auto'   # OPTIONS: ('none', 'auto', or 'manual')
#                         # none - does not trim data at all
#                         # auto - trims ATL03 track to extent of truth bounding region
#                         # manual - trims ATL03 track by latitude or time
#                             # Example: 'manual,lat,38,39'
#                             # Only uses data between latitude 38 and 39 deg
#                             # Example: 'manual,time,3,4'
#                             # Only uses data between time 3 and 4 seconds
                                        
#     createAtl03LasFile = True     # Option to create output measured ATL03 .las file
#     createAtl03KmlFile = True     # Option to create output measured ATL03 .kml file
#     createAtl03CsvFile = False    # Option to create output measured ATL03 .csv file
#     createAtl08KmlFile = False    # Option to create output measured ATL08 .kml file
#     createAtl08CsvFile = False    # Option to create output measured ATL08 .csv file
        
#     ##### End Inputs for getAtlMeasuredSwath
    
    
#     ##### Start Inputs for getAtlTruthSwath
    
#     buffer = 50                 # Distance in cross-track (meters) around ATL03 track to look for truth data 
#     useExistingTruth = True     # Option to use existing truth data if it exists
#     # truthSwathDir = '//lidar-server/lidar/USERS/mike/iceSat2/truth_data/'  # Path to existing truth data
#     # truthSwathDir = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/sonoma_ATL03_20181030110205_04860106_001_01_gt2r_TRUTH_50L50Rm_buffer_classified_short_dtm_utm.tif'
#     truthSwathDir = '//lidar-server/lidar/USERS/mike/iceSat2/truth_data/wsmr_ATL03_20190928175636_00280506_R002_01_sreq_3002_gt2r_TRUTH_50L50Rm_buffer.las'
# #    truthSwathDir = 'N:/USERS/mike/iceSat2/atl03_validation/r002_ncarolina_20181018/ATL03_20181018202814_03090102_002_01_gt1r_REFERENCE_50L50Rm_buffer.las'
#     # truthSwathDir = '//lidar-server/lidar/USERS/mike/iceSat2/truth_data/wsmr_ATL03_20191031162419_05310506_R002_01_sub_344_gt2r_TRUTH_150L150Rm_buffer.las'
#     # truthSwathDir = '//lidar-server/lidar/USERS/mike/iceSat2/atl03_validation/test/sonoma_ATL03_20181030110205_04860106_001_01_gt2r_TRUTH_50L50Rm_buffer_classified_short.las'
#     truthFileType = '.las'
#     createTruthFile = False      # Option to create output truth .las file
    
#     ##### End Inputs for getAtlTruthSwath
    
    
#     ##### Start Inputs for getMeasurementError
    
#     offsetsCrossTrackBounds = np.array([-50,50])      # Cross-track limits to search for geolocation error
#     offsetsAlongTrackBounds = np.array([-50,50])      # Along-track limits to search for geolocation error
#     offsetsRasterResolutions = np.array([8, 4, 2, 1])  # Multi-resolutional step-down raster resolutions (in meters)
#     refHeightType = 'HAE'              # 'HAE' or 'MSL'
#     offsetsUseVerticalShift = False    # Option to use a vertical shift
#     offsetsVerticalShift = 0           # Vertical shift to use if above set to True (in meters)
#     useMeasSigConf = True             # Use measured signal confidence (or use ground truth)
#                                       # Meas Classes (0 = Unclass, 1 = Ground, 2 = Low Veg, 3 = High Veg), Texpert Truth Classes (0 = Unclass, 2 = Ground, 4 = Veg, 6 = Building)
#     filterData = [3,4]               # Signal Confidence (0, 1, 2, 3, 4)
#     offsets = offsetsStruct(offsetsCrossTrackBounds, offsetsAlongTrackBounds, offsetsRasterResolutions, offsetsUseVerticalShift, offsetsVerticalShift )
#     createMeasCorrFile = True     # Option to create ouput measured corrected .las file
#     makePlots = True              # Option to make output plots
#     showPlots = False             # Option to show output plot windows
      
#     ##### End Inputs for getMeasurementError  
    
    
#     ##### CODE BELOW -- DO NOT EDIT ###############################################
    
#     timeStart = runTime.time()
    
#     # Call getAtlMeasuredSwath
#     print('RUNNING getAtlMeasuredSwath...\n')
#     atl03Data, atl08Data, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, gtNum, trimInfo, createAtl03LasFile, createAtl03KmlFile, createAtl08KmlFile, createAtl03CsvFile, createAtl08CsvFile)
          
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
#     atlTruthData = getAtlTruthSwath(atl03Data, rotationData, 
#                                     truthHeaderDF, truthFilePaths,
#                                     buffer, outFilePath, createTruthFile, 
#                                     truthFileType, useExistingTruth, logFileID=False)
    
#     # Call getMeasurementError
#     print('RUNNING getMeasurementError...\n')
#     atlCorrections = getMeasurementError(atl03Data, atlTruthData, refHeightType, 
#                                           rotationData, outFilePath, useMeasSigConf, 
#                                           filterData, offsets, createMeasCorrFile, 
#                                           makePlots, showPlots, logFileID=False)

#     # End timer
#     timeEnd = runTime.time()
#     timeElapsedTotal = timeEnd - timeStart
#     timeElapsedMin = np.floor(timeElapsedTotal / 60)
#     timeElapsedSec = timeElapsedTotal % 60
        
#     # Print completion message
#     print('   Script Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
#     print('\n')
    
# endIf
    
    
