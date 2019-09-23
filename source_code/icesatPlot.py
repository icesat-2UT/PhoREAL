# -*- coding: utf-8 -*-
"""
This script that provides basic plotting functionality for PhoREAL

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
import matplotlib.pyplot as plt
import os
import numpy as np


# Function to make contour plot
def plotContour(resultsCrossTrackShift, resultsAlongTrackShift, resultsMAE, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, counter):
    
    # Make contour plot
    plt.figure()
    plt.contour(resultsCrossTrackShift, resultsAlongTrackShift, resultsMAE, colors = 'black')
    plt.contourf(resultsCrossTrackShift, resultsAlongTrackShift, resultsMAE, cmap = 'viridis')
    legendStr = 'Min MAE = ' +  '{:0.3f}'.format(atlCorrections.mae[0]) + ' m'
    plt.plot(atlCorrections.crossTrack, atlCorrections.alongTrack, 'ro', markersize = 7, markerfacecolor = 'r', label = legendStr)
    plt.axis('equal')
    plt.axis('tight')
    plt.grid(b = True, which = 'major', axis = 'both')
    plt.xlabel('Cross-Track Offset (m)')
    plt.ylabel('Along-Track Offset (m)')
    title = atlMeasuredData.atl03FileName + ' (' + atlMeasuredData.gtNum + '): ' + atlMeasuredData.trackDirection + '\n' + \
            'Observed Data Correction (' + str(rasterResolution) + ' m Raster)\n' + \
            'Easting, Northing, Vertical: ' + '{:0.1f}'.format(atlCorrections.easting[0]) + ', ' + '{:0.1f}'.format(atlCorrections.northing[0]) + ', ' + '{:0.1f}'.format(atlCorrections.z[0]) + ' m \n' + \
            'Cross-Track, Along-Track: ' + '{:0.0f}'.format(atlCorrections.crossTrack[0]) + ', ' + '{:0.0f}'.format(atlCorrections.alongTrack[0]) + ' m'
    plt.title(title, fontsize = 10, fontweight = 'bold')
    plt.legend(loc = 'upper left')
    cbar = plt.colorbar();
    cbar.set_label('Vertical Mean Absolute Error (m)')
    
    # Save plot
    if(not os.path.exists(os.path.normpath(outFilePath))):
        os.mkdir(os.path.normpath(outFilePath))
    # EndIf
    plotNum = counter + 1
    outName = atlMeasuredData.atl03FileName + '_' + atlMeasuredData.gtNum + '_figContour_' + '{:0.0f}'.format(plotNum) + '.png'
    outPath = os.path.normpath(outFilePath + '/' + outName)
    plt.savefig(outPath)
    
    # Show plot
    if(showPlots):
        plt.show()
    # EndIf
    
    
# Function to make Z vs Y plot
def plotZY(measRasterYCommonFinal, measRasterZCommonFinal, truthRasterYCommonFinal, truthRasterZCommonFinal, zErrorCommonFinal, segmentErrorY, segmentErrorZ, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, counter):
    
    # Define plot colors
    myBlue = (0/255, 114/255, 178/255)
    myOrange = (213/255, 94/255, 0/255)
    myYellow = (230/255, 159/255, 0/255)
    
    # Make Z vs Y subplot 1
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.scatter(truthRasterYCommonFinal/1000, truthRasterZCommonFinal, color = myBlue, label = 'Truth Mean Raster')
    ax1.scatter(measRasterYCommonFinal/1000, measRasterZCommonFinal, color = myOrange, label = 'ATL03 Mean Raster')
    ax1.axis('tight')
    ax1.grid(b = True, which = 'major', axis = 'both')
    ax1.set_ylabel('Z (m)')
    ax1.legend(loc = 'upper left')
    titleStr = 'Truth and Corrected Measured Position Data (' + str(rasterResolution) + ' m Raster)'
    ax1.set_title(titleStr)
    
    # Make Z Error vs Y subplot 2
    ax2.scatter(truthRasterYCommonFinal/1000, zErrorCommonFinal, color = myBlue, label = 'Z Error')
    ax2.plot(segmentErrorY/1000, segmentErrorZ, color = myYellow, label = '100 m Segment Mean Error')
    ax2.axis('tight')
    ax2.grid(b = True, which = 'major', axis = 'both')
    ax2.set_ylabel('Z Error (m)')
    ax2.set_xlabel('UTM Northing (km)')
    ax2.legend(loc = 'upper left')
    titleStr = 'MAE = ' + '{:0.2f}'.format(atlCorrections.mae[0]) + ' m, RMSE = ' + '{:0.2f}'.format(atlCorrections.rmse[0]) + ' m, Mean Error = ' + '{:0.2f}'.format(atlCorrections.me[0]) + ' m'
    ax2.set_title(titleStr)
    
    # Save plot
    if(not os.path.exists(os.path.normpath(outFilePath))):
        os.mkdir(os.path.normpath(outFilePath))
    # EndIf
    plotNum = counter + 1
    outName = atlMeasuredData.atl03FileName + '_' + atlMeasuredData.gtNum + '_figZY_' + '{:0.0f}'.format(plotNum) + '.png'
    outPath = os.path.normpath(outFilePath + '/' + outName)
    plt.savefig(outPath)
    
    # Show plot
    if(showPlots):
        plt.show()
    # EndIf
    
    
# Function to make Z vs T plot
def plotZT(measRasterTCommonFinal, measRasterZCommonFinal, zErrorCommonFinal, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, counter):
    
    # Define plot colors
    myBlue = (0/255, 114/255, 178/255)
    myOrange = (213/255, 94/255, 0/255)
    
    # Make Z vs Y subplot 1
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.scatter(measRasterTCommonFinal, measRasterZCommonFinal, color = myOrange)
    ax1.axis('tight')
    ax1.grid(b = True, which = 'major', axis = 'both')
    ax1.set_ylabel('Z (m)')
    titleStr = 'Measured Time Data (' + str(rasterResolution) + ' m Raster)'
    ax1.set_title(titleStr)
    
    # Make Z Error vs Y subplot 2
    ax2.scatter(measRasterTCommonFinal, zErrorCommonFinal, color = myBlue)
    ax2.axis('tight')
    ax2.grid(b = True, which = 'major', axis = 'both')
    ax2.set_ylabel('Z Error (m)')
    ax2.set_xlabel('Time (sec)')
    titleStr = 'MAE = ' + '{:0.2f}'.format(atlCorrections.mae[0]) + ' m, RMSE = ' + '{:0.2f}'.format(atlCorrections.rmse[0]) + ' m, Mean Error = ' + '{:0.2f}'.format(atlCorrections.me[0]) + ' m'
    ax2.set_title(titleStr)
    
    # Save plot
    if(not os.path.exists(os.path.normpath(outFilePath))):
        os.mkdir(os.path.normpath(outFilePath))
    # EndIf
    plotNum = counter + 1
    outName = atlMeasuredData.atl03FileName + '_' + atlMeasuredData.gtNum + '_figZT_' + '{:0.0f}'.format(plotNum) + '.png'
    outPath = os.path.normpath(outFilePath + '/' + outName)
    plt.savefig(outPath)
    
    # Show plot
    if(showPlots):
        plt.show()
    # EndIf
    
# Function to plot any x,y data from getAtlMeasuredSwath_auto.py
def getPlot(xData, yData, xLabel, yLabel, title, filterType = [], filterData = [], filterNum = []):
    
    # Define plot colors
    myColors = np.array([[0/255, 114/255, 178/255],[230/255, 159/255, 0/255],[240/255, 228/255, 66/255],[0/255, 158/255, 115/255],[213/255, 94/255, 0/255]])
    
    # Color Blind-friendly colors: https://www.nature.com/articles/nmeth.1618
    # Blue
    # Orange
    # Yellow
    # Green
    # Red

    # Open figure window
    fig1 = plt.figure()
    
    if(np.size(filterData) > 0):
        
        # Loop through filter numbers
        for i in range(0,np.size(filterNum)):
            
            # Filter data
            matchVal = filterNum[i]
            matchingInds = filterData == matchVal
            xDataPlot = xData[matchingInds]
            yDataPlot = yData[matchingInds]
            
            if('class' in filterType.lower()):
                if(matchVal==0):
                    matchStr = 'Unclassified'
                elif(matchVal==1):
                    matchStr = 'Ground'
                elif(matchVal==2):
                    matchStr = 'Low Veg'
                elif(matchVal==3):
                    matchStr = 'High Veg'
                # endIf
            else:
                matchStr = matchVal
            # endIf
            
            # Plot filtered data
            plt.scatter(xDataPlot, yDataPlot, color = myColors[matchVal,:], label = matchStr, s=0.5)
        
        # endFor
        
        plt.legend(loc = 'upper left')
    
    else:
        
        # Plot all data
        plt.scatter(xData, yData, color = myColors[0,:], s=0.5)
    
    # endIf
    
    # Figure properties
    plt.axis('tight')
    plt.grid(b = True, which = 'major', axis = 'both')
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(title)
    fig1.show()