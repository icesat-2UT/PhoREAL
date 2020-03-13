# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 13:14:46 2019

This module will take in ATL03 MEASURED and TRUTH data and find the 
geolocation error in the MEASURED data based on TRUTH data.

@author: malonzo
"""

# Import modules
import os
import numpy as np
import time as runTime
import sys
import warnings
from getAtlMeasuredSwath_auto import getAtlMeasuredSwath
from getAtlTruthSwath_auto import getAtlTruthSwath
from icesatIO import writeLas
from icesatUtils import (ismember, getRaster, getIntersection, getCoordRotRev)
from icesatPlot import (plotContour, plotZY, plotZT)


class offsetsStruct:
    
    # Define class with designated fields
    def __init__(self, crossTrackBounds, alongTrackBounds, rasterResolutions, useVerticalShift, verticalShift):
        
        self.crossTrackBounds = crossTrackBounds
        self.alongTrackBounds = alongTrackBounds
        self.rasterResolutions = rasterResolutions
        self.useVerticalShift = useVerticalShift
        self.verticalShift = verticalShift


class atlMeasuredDataReducedStruct:
    
    # Define class with designated fields
    def __init__(self, easting, northing, z, crossTrack, alongTrack, time, classification, signalConf):
        
        self.easting = easting
        self.northing = northing
        self.z = z
        self.crossTrack = crossTrack
        self.alongTrack = alongTrack
        self.time = time
        self.classification = classification
        self.signalConf = signalConf


class atlCorrectionsStruct:
    
    # Define class with designated fields
    def __init__(self, correctionsCrossTrack, correctionsAlongTrack, correctionsEasting, correctionsNorthing, correctionsVertical, correctionsMAE, correctionsRMSE, correctionsME):
        
        self.crossTrack = correctionsCrossTrack
        self.alongTrack = correctionsAlongTrack
        self.easting = correctionsEasting
        self.northing = correctionsNorthing
        self.z = correctionsVertical
        self.mae = correctionsMAE
        self.rmse = correctionsRMSE
        self.me = correctionsME
        
        
def getMeasurementError(atlMeasuredData, atlTruthData, rotationData, outFilePath, useMeasSigConf, filterData, offsets, createMeasCorrFile, makePlots, showPlots):
    
    # Start timer
    timeStart = runTime.time()
    
    # Get and print ground track number
    gtNum = atlMeasuredData.gtNum
    print('   Ground Track Number: %s' % gtNum)

    # Turn off runtime warnings (from np.nanmean)
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    # EndIf
    
    # Get array of offsets
    crossTrackBounds = offsets.crossTrackBounds
    alongTrackBounds = offsets.alongTrackBounds
    rasterResolutions = offsets.rasterResolutions
    
    # Get reduced MEASURED data (limited to bounds of TRUTH data)
    measInTruthInds = (atlMeasuredData.alongTrack >= atlTruthData.alongTrack.min()) & (atlMeasuredData.alongTrack <= atlTruthData.alongTrack.max())
    atlMeasuredDataReducedX = atlMeasuredData.easting[measInTruthInds]
    atlMeasuredDataReducedY = atlMeasuredData.northing[measInTruthInds]
    atlMeasuredDataReducedZ = atlMeasuredData.z[measInTruthInds]
    atlMeasuredDataReducedXrot = atlMeasuredData.crossTrack[measInTruthInds]
    atlMeasuredDataReducedYrot = atlMeasuredData.alongTrack[measInTruthInds]
    atlMeasuredDataReducedTime = atlMeasuredData.time[measInTruthInds]
    atlMeasuredDataReducedClass = atlMeasuredData.classification[measInTruthInds]
    atlMeasuredDataReducedSignal = atlMeasuredData.signalConf[measInTruthInds]
    
    # Get MEASURED reduced class
    atlMeasuredDataReduced = atlMeasuredDataReducedStruct(atlMeasuredDataReducedX, \
                                                          atlMeasuredDataReducedY, \
                                                          atlMeasuredDataReducedZ, \
                                                          atlMeasuredDataReducedXrot, \
                                                          atlMeasuredDataReducedYrot, \
                                                          atlMeasuredDataReducedTime, \
                                                          atlMeasuredDataReducedClass, \
                                                          atlMeasuredDataReducedSignal)
    
    # Filter Ground data from ATL08 if available
    if(useMeasSigConf):
    
        print('   Using Signal Confidence for Measurement Offset Computation')
        
        # Filter MEASURED data points based on signal confidence
        measGroundInds = ismember(atlMeasuredDataReduced.signalConf,filterData)[0]
        measXRot = atlMeasuredDataReduced.crossTrack[measGroundInds]
        measYRot = atlMeasuredDataReduced.alongTrack[measGroundInds]
        measZ = atlMeasuredDataReduced.z[measGroundInds]
        measT = atlMeasuredDataReduced.time[measGroundInds]
        
        # Do not filter ground TRUTH data points
        truthXRot = atlTruthData.crossTrack
        truthYRot = atlTruthData.alongTrack
        truthZ = atlTruthData.z
    
    else:
    
        print('   Using Ground Truth for Measurement Offset Computation')
        
        # Filter ground MEASURED data points
        icesatGroundIndex = 1
        measGroundInds = ismember(atlMeasuredDataReduced.classification,icesatGroundIndex)[0]
        measXRot = atlMeasuredDataReduced.crossTrack[measGroundInds]
        measYRot = atlMeasuredDataReduced.alongTrack[measGroundInds]
        measZ = atlMeasuredDataReduced.z[measGroundInds]
        measT = atlMeasuredDataReduced.time[measGroundInds]
        
        # Filter ground TRUTH data points
        truthGroundInds = ismember(atlTruthData.classification,filterData)[0]
        truthXRot = atlTruthData.crossTrack[truthGroundInds]
        truthYRot = atlTruthData.alongTrack[truthGroundInds]
        truthZ = atlTruthData.z[truthGroundInds]
        
    # endIf
        
    if((len(crossTrackBounds) > 1) and (len(alongTrackBounds) > 1)):
            
        # Initialize parameters
        correctionsCrossTrack = 0
        correctionsAlongTrack = 0

        crossTrackEdgeError = 0
        alongTrackEdgeError = 0
        numErr = 0
        crossTrackBoundLimits = [-48, 48]
        alongTrackBoundLimits = [-1000, 1000]

        def limBounds(x, y, lim):
            # x = max([x, lim[0]])
            # y = min([lim[1], y])
            return x, y
        
        # Get array offsets
        crossTrackMinBound = (np.fix(crossTrackBounds[0]/rasterResolutions[0])*rasterResolutions[0]).astype(int)
        crossTrackMaxBound = (np.fix((crossTrackBounds[-1])/rasterResolutions[0])*rasterResolutions[0]).astype(int)
        crossTrackOffsetsOrig = np.arange(crossTrackMinBound, crossTrackMaxBound + 1, rasterResolutions[0])
        alongTrackMinBound = (np.fix(alongTrackBounds[0]/rasterResolutions[0])*rasterResolutions[0]).astype(int)
        alongTrackMaxBound = (np.fix((alongTrackBounds[-1])/rasterResolutions[0])*rasterResolutions[0]).astype(int)
        alongTrackOffsetsOrig = np.arange(alongTrackMinBound, alongTrackMaxBound + 1, rasterResolutions[0])
        
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
                
                # Set cross-track/along-track offsets
                # crossTrackOffsets = crossTrackOffsetsOrig
                # alongTrackOffsets = alongTrackOffsetsOrig

                crossTrackBoundStart = crossTrackBounds[0] + crossTrackEdgeError
                crossTrackBoundEnd = crossTrackBounds[-1] + crossTrackEdgeError
                crossTrackBoundStart, crossTrackBoundEnd = limBounds(crossTrackBoundStart, crossTrackBoundEnd, crossTrackBoundLimits)

                crossTrackMinBound = (np.fix(crossTrackBoundStart/rasterResolutions[0])*rasterResolutions[0]).astype(int)
                crossTrackMaxBound = (np.fix(crossTrackBoundEnd/rasterResolutions[0])*rasterResolutions[0]).astype(int)
                crossTrackOffsets = np.arange(crossTrackMinBound, crossTrackMaxBound + 1, rasterResolutions[0])

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
            print('   Cross Track Offsets = [%s, %s]' %(crossTrackOffsets[0], crossTrackOffsets[-1]))
            print('   Along Track Offsets = [%s, %s]' %(alongTrackOffsets[0], alongTrackOffsets[-1]))
            
            # Make Offsets Int's and flip along-track offsets array so larger numbers at top
            crossTrackOffsets = crossTrackOffsets.astype(int)
            alongTrackOffsets = (np.flipud(alongTrackOffsets)).astype(int)
            
            # Reduce TRUTH data to width of cross-track offsets
            measWidth = (np.ceil( ( measXRot.max() - measXRot.min() ) / 2.0 )).astype(int)
    #        minTruthBound = (crossTrackOffsets.min() - measWidth).astype(int)
    #        maxTruthBound = (crossTrackOffsets.max() + measWidth).astype(int)
            minTruthBound = (crossTrackOffsets.min() - measWidth).astype(int) - 5
            maxTruthBound = (crossTrackOffsets.max() + measWidth).astype(int) + 5
            truthInOffsetsInds = (truthXRot >= minTruthBound) & (truthXRot <= maxTruthBound)
            truthXRotReduced = truthXRot[truthInOffsetsInds]
            truthYRotReduced = truthYRot[truthInOffsetsInds]
            truthZReduced = truthZ[truthInOffsetsInds]
            
            # Rasterize MEASURED data
            gridMethod = 'Mean'
            fillValue = -999
            print('   Gridding Measured Data at %s m Resolution using %s Values...' % (rasterResolution, gridMethod))
            measRasterRot = getRaster(measXRot, measYRot, measZ, rasterResolution, gridMethod, fillValue, measT)
            
            # Rasterize TRUTH data
            print('   Gridding Truth Data at %s m Resolution using %s Values...' % (rasterResolution, gridMethod))
            truthRasterRot = getRaster(truthXRotReduced, truthYRotReduced, truthZReduced, rasterResolution, gridMethod, fillValue)
            
            # Find offsets with minimum MAE
            print('   Finding Offsets with Minimum MAE...')
            
            # Remove NaN data (-999) in MEASURED raster
            measIndsToKeep = measRasterRot.grid != -999
            measRasterXRot = measRasterRot.x[measIndsToKeep]
            measRasterYRot = measRasterRot.y[measIndsToKeep]
            measRasterZ    = measRasterRot.grid[measIndsToKeep]
            measRasterT    = measRasterRot.t[measIndsToKeep]
            
            # Remove NaN data (-999) in TRUTH raster
            truthIndsToKeep = truthRasterRot.grid != -999
            truthRasterXRot = truthRasterRot.x[truthIndsToKeep]
            truthRasterYRot = truthRasterRot.y[truthIndsToKeep]
            truthRasterZ    = truthRasterRot.grid[truthIndsToKeep]

            
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
                    
                    # Get new MEASURED line at new CT/AT offsets
                    measRasterXRotShifted = measRasterXRot + crossTrackShift
                    measRasterYRotShifted = measRasterYRot + alongTrackShift
                    
                    # Get common TRUTH and MEASURED raster points
                    xyMeasRasterRotShifted = np.concatenate((np.c_[measRasterXRotShifted],np.c_[measRasterYRotShifted]),axis = 1)
                    xyTruthRasterRot = np.concatenate((np.c_[truthRasterXRot],np.c_[truthRasterYRot]),axis = 1)
                    commonVals, xyCommonMeasInds, xyCommonTruthInds = getIntersection(xyMeasRasterRotShifted, xyTruthRasterRot)
            
                    # Find percent common points (Require 20% min to process data)
                    pctCommon = len(xyCommonMeasInds)/len(xyMeasRasterRotShifted)
            
                    # Process data if there are common TRUTH and MEASURED points
                    if(commonVals.any() and pctCommon>=0.2):
                        
                        # Get common TRUTH and MEASURED Z points
                        measRasterZCommon = measRasterZ[xyCommonMeasInds]
                        truthRasterZCommon = truthRasterZ[xyCommonTruthInds]
                        
                        # Find Z error
                        zError = truthRasterZCommon - measRasterZCommon
                        
                        # Determine vertical shift
                        if(offsets.useVerticalShift == 1):
                            
                            # Use defined vertical shift
                            verticalShift = offsets.verticalShift
                            
                        else:
                            
                            # MAE is tied to median error (RMSE to mean error)
                            verticalShift = np.median(zError)
                            
                        # Endif
                        
                        # Apply vertical shift
                        measRasterZShiftedCommon = measRasterZCommon + verticalShift
                        
                        # Get new Z error
                        zError = truthRasterZCommon - measRasterZShiftedCommon
                        
                        # Get Mean Absolute Error and RMSE
                        MAE = np.mean(abs(zError))
                        RMSE = np.sqrt(np.mean(zError**2))
                        ME = np.mean(zError)
                        
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
                    
                # endIf
            # endIf
        
    #        # CHANGE AFTER DEMO!!!
    #        # Set first/last rows/columns to NaN on first iteration to avoid choosing edge cases           
    #        if(i==0):
    #            resultsMAE[0,:] = np.nan
    #            resultsMAE[-1,:] = np.nan
    #            resultsMAE[:,0] = np.nan
    #            resultsMAE[:,-1] = np.nan
    #        # endIf
            
            # Find row, col for min MAE
            minRow = np.where(resultsMAE == np.nanmin(resultsMAE))[0]
            minCol = np.where(resultsMAE == np.nanmin(resultsMAE))[1]
            if minRow == [] or minCol == []:
                print('Error: Gridding all NaN, ensure truth/meas files are correct')
                return []

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
            correctionsCrossTrack = resultsCrossTrackShift[minRow,minCol]
            correctionsAlongTrack = resultsAlongTrackShift[minRow,minCol]
            correctionsVertical = resultsVerticalShift[minRow,minCol]
            correctionsMAE = resultsMAE[minRow,minCol]
            correctionsRMSE = resultsRMSE[minRow,minCol]
            correctionsME = resultsME[minRow,minCol]
        
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
            
            # Store data in class structure
            atlCorrections = atlCorrectionsStruct(correctionsCrossTrack, correctionsAlongTrack, \
                                                  correctionsEasting, correctionsNorthing, \
                                                  correctionsVertical, \
                                                  correctionsMAE, correctionsRMSE, correctionsME)
            
            # Get final MEASURED raster data
            measRasterXRotFinal= measRasterXRot + atlCorrections.crossTrack
            measRasterYRotFinal = measRasterYRot + atlCorrections.alongTrack
            measRasterZFinal = measRasterZ + atlCorrections.z
            measRasterTFinal = measRasterT
    
            # Get common MEASURED and TRUTH raster indices 
            xyMeasRasterRotFinal = np.concatenate((np.c_[measRasterXRotFinal],np.c_[measRasterYRotFinal]),axis = 1)
            commonValsFinal, xyCommonMeasIndsFinal, xyCommonTruthIndsFinal = getIntersection(xyMeasRasterRotFinal, xyTruthRasterRot)
            
            # Get common MEASURED raster points
            measRasterXRotCommonFinal = measRasterXRotFinal[xyCommonMeasIndsFinal]
            measRasterYRotCommonFinal = measRasterYRotFinal[xyCommonMeasIndsFinal]
            measRasterZCommonFinal = measRasterZFinal[xyCommonMeasIndsFinal]
            measRasterTCommonFinal = measRasterTFinal[xyCommonMeasIndsFinal]
            
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
            segmentError = getRaster(measRasterXRotCommonFinal, measRasterYRotCommonFinal, zErrorCommonFinal, errorResolution, gridMethod, fillValue)
            segmentErrorX = np.mean(segmentError.x, axis = 1)
            segmentErrorY = np.mean(segmentError.y, axis = 1)
            segmentErrorZ = np.nanmean(segmentError.grid, axis = 1)
            
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
            _, measRasterYCommonFinal, _, _, _ = getCoordRotRev(measRasterXRotCommonFinal, measRasterYRotCommonFinal, R_mat, xRotPt, yRotPt)
            _, truthRasterYCommonFinal, _, _, _ = getCoordRotRev(truthRasterXRotCommonFinal, truthRasterYRotCommonFinal, R_mat, xRotPt, yRotPt)
    
            # Make plots
            if(makePlots):
                
                # Make contour plot
                print('   Making Plots...')
                plotContour(resultsCrossTrackShift, resultsAlongTrackShift, resultsMAE, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, i)
            
                # Make ZY scatter plots
                plotZY(measRasterYCommonFinal, measRasterZCommonFinal, truthRasterYCommonFinal, truthRasterZCommonFinal, zErrorCommonFinal, segmentErrorY_plot, segmentErrorZ_plot, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, i)
                
                # Make ZT scatter plots
                plotZT(measRasterTCommonFinal, measRasterZCommonFinal, zErrorCommonFinal, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, i)
    
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
        print('   Gridding Measured Data at %s m Resolution using %s Values...' % (rasterResolution, gridMethod))
        
        # Add in corrections
        measXRot += correctionsCrossTrack
        measYRot += correctionsAlongTrack
        measZ += correctionsVertical
        
        # Raster MEASURED data
        measRasterRot = getRaster(measXRot, measYRot, measZ, rasterResolution, gridMethod, fillValue, measT)
        
        # Reduce TRUTH data to width of cross-track offsets
        measWidth = (np.ceil( ( measXRot.max() - measXRot.min() ) / 2.0 )).astype(int)
        minTruthBound = (crossTrackBounds - measWidth).astype(int) - 5
        maxTruthBound = (crossTrackBounds + measWidth).astype(int) + 5
        truthInOffsetsInds = (truthXRot >= minTruthBound) & (truthXRot <= maxTruthBound)
        truthXRotReduced = truthXRot[truthInOffsetsInds]
        truthYRotReduced = truthYRot[truthInOffsetsInds]
        truthZReduced = truthZ[truthInOffsetsInds]
            
        # Rasterize TRUTH data
        print('   Gridding Truth Data at %s m Resolution using %s Values...' % (rasterResolution, gridMethod))
        truthRasterRot = getRaster(truthXRotReduced, truthYRotReduced, truthZReduced, rasterResolution, gridMethod, fillValue)
        
        # Remove NaN data (-999) in MEASURED raster
        measIndsToKeep = measRasterRot.grid != -999
        measRasterXRot = measRasterRot.x[measIndsToKeep]
        measRasterYRot = measRasterRot.y[measIndsToKeep]
        measRasterZ    = measRasterRot.grid[measIndsToKeep]
        measRasterT    = measRasterRot.t[measIndsToKeep]
        
        # Get final MEASURED raster data
        measRasterXRotFinal= measRasterXRot
        measRasterYRotFinal = measRasterYRot
        measRasterZFinal = measRasterZ
        measRasterTFinal = measRasterT
            
        # Remove NaN data (-999) in TRUTH raster
        truthIndsToKeep = truthRasterRot.grid != -999
        truthRasterXRot = truthRasterRot.x[truthIndsToKeep]
        truthRasterYRot = truthRasterRot.y[truthIndsToKeep]
        truthRasterZ    = truthRasterRot.grid[truthIndsToKeep]
            
        # Get common MEASURED and TRUTH raster indices 
        xyMeasRasterRotFinal = np.concatenate((np.c_[measRasterXRotFinal],np.c_[measRasterYRotFinal]),axis = 1)
        xyTruthRasterRot = np.concatenate((np.c_[truthRasterXRot],np.c_[truthRasterYRot]),axis = 1)
        commonValsFinal, xyCommonMeasIndsFinal, xyCommonTruthIndsFinal = getIntersection(xyMeasRasterRotFinal, xyTruthRasterRot)
        
        # Get common MEASURED raster points
        measRasterXRotCommonFinal = measRasterXRotFinal[xyCommonMeasIndsFinal]
        measRasterYRotCommonFinal = measRasterYRotFinal[xyCommonMeasIndsFinal]
        measRasterZCommonFinal = measRasterZFinal[xyCommonMeasIndsFinal]
        measRasterTCommonFinal = measRasterTFinal[xyCommonMeasIndsFinal]
        
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
        segmentError = getRaster(measRasterXRotCommonFinal, measRasterYRotCommonFinal, zErrorCommonFinal, errorResolution, gridMethod, fillValue)
        segmentErrorX = np.mean(segmentError.x, axis = 1)
        segmentErrorY = np.mean(segmentError.y, axis = 1)
        segmentErrorZ = np.nanmean(segmentError.grid, axis = 1)
        
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
        _, measRasterYCommonFinal, _, _, _ = getCoordRotRev(measRasterXRotCommonFinal, measRasterYRotCommonFinal, R_mat, xRotPt, yRotPt)
        _, truthRasterYCommonFinal, _, _, _ = getCoordRotRev(truthRasterXRotCommonFinal, truthRasterYRotCommonFinal, R_mat, xRotPt, yRotPt)


        MAE = np.mean(abs(zErrorCommonFinal))
        RMSE = np.sqrt(np.mean(zErrorCommonFinal**2))
        ME = np.mean(zErrorCommonFinal)
                        
        # Set error values
        correctionsMAE = np.array([MAE])
        correctionsRMSE = np.array([RMSE])
        correctionsME = np.array([ME])
        
        # Generate atlCorrections struct
        atlCorrections = atlCorrectionsStruct(correctionsCrossTrack, correctionsAlongTrack, \
                                              correctionsEasting, correctionsNorthing, \
                                              correctionsVertical, \
                                              correctionsMAE, correctionsRMSE, correctionsME)
                        
        # Make plots
        if(makePlots):
            
            # Make contour plot
            print('   Making Plots...')
        
            # Plot counter
            i = 0
        
            # Make ZY scatter plots
            plotZY(measRasterYCommonFinal, measRasterZCommonFinal, truthRasterYCommonFinal, truthRasterZCommonFinal, zErrorCommonFinal, segmentErrorY_plot, segmentErrorZ_plot, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, i)
            
            # Make ZT scatter plots
            plotZT(measRasterTCommonFinal, measRasterZCommonFinal, zErrorCommonFinal, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, i)

        # endIf 
            
    # endIf
    
    # Create output file
    if(createMeasCorrFile): 
        
        # Get measured corrected data
        atlMeasCorrX = atlMeasuredData.easting + atlCorrections.easting
        atlMeasCorrY = atlMeasuredData.northing + atlCorrections.northing
        atlMeasCorrZ = atlMeasuredData.z + atlCorrections.z
        
        # Get E,N,U text for output file name
        xDirTxt = 'e' # east
        yDirTxt = 'n' # north
        zDirTxt = 'u' # up
        
        # Get directional text strings
        if(atlCorrections.easting < 0):
            xDirTxt = 'w' # west
        # endIf
        
        if(atlCorrections.northing < 0):
            yDirTxt = 's' # south
        # endIf
        
        if(atlCorrections.z < 0):
            zDirTxt = 'd' # down
        # endIf
        
        # Get measured corrected string text
        enuTxt = xDirTxt + '{0:0.1f}'.format(abs(atlCorrections.easting[0])) + '_' + \
                 yDirTxt + '{0:0.1f}'.format(abs(atlCorrections.northing[0])) + '_' + \
                 zDirTxt + '{0:0.1f}'.format(abs(atlCorrections.z[0]))
                 
        # Get output path name
        print('')
        print('   Writing measured corrected .las file...')
        outName = atlMeasuredData.atl03FileName + '_' + atlMeasuredData.gtNum + '_SHIFTED_' + enuTxt + '.las'
        outPath = os.path.normpath(outFilePath + '/' + outName)
        
        # If output directory does not exist, create it
        if(not os.path.exists(os.path.normpath(outFilePath))):
            os.mkdir(os.path.normpath(outFilePath))
        # endIf
        
        # Write .las file
        writeLas(np.ravel(atlMeasCorrX),np.ravel(atlMeasCorrY),np.ravel(atlMeasCorrZ),'utm',outPath,np.ravel(atlMeasuredData.classification),np.ravel(atlMeasuredData.intensity),np.ravel(atlMeasuredData.signalConf),atlMeasuredData.hemi,atlMeasuredData.zone)
    
    # endIf
    
    # End timer
    timeEnd = runTime.time()
    timeElapsedTotal = timeEnd - timeStart
    timeElapsedMin = np.floor(timeElapsedTotal / 60)
    timeElapsedSec = timeElapsedTotal % 60
    
    # Print completion message
    print('   Module Completed in %d min %d sec.' % (timeElapsedMin, timeElapsedSec))
    print('\n')
    
    # Return output data
    return atlCorrections


# Unit test
if __name__ == "__main__":
    
    print('Test')
