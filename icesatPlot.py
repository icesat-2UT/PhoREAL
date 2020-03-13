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
from icesatIO import writeArrayToCSV
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import glob
import sys
import warnings
from matplotlib.backend_tools import ToolToggleBase
plt.rcParams['toolbar'] = 'toolmanager'


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
    plt.tight_layout()
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
    
    # Close figure
    plt.close()
    
    
# Function to make Z vs Y plot
def plotZY(measRasterYCommonFinal, measRasterZCommonFinal, truthRasterYCommonFinal, truthRasterZCommonFinal, zErrorCommonFinal, segmentErrorY, segmentErrorZ, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, counter):
    
    # Define plot colors
    myBlue = (0/255, 114/255, 178/255)
    myOrange = (213/255, 94/255, 0/255)
    myYellow = (230/255, 159/255, 0/255)
    
    # Make Z vs Y subplot 1
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.scatter(measRasterYCommonFinal/1000, measRasterZCommonFinal, color = myOrange, s=0.7, label = 'Shifted ICESat-2 Mean Raster', zorder=1)
    ax1.scatter(truthRasterYCommonFinal/1000, truthRasterZCommonFinal, color = myBlue, s=0.7, label = 'Reference Mean Raster', zorder=0)
    ax1.axis('tight')
    ax1.grid(b = True, which = 'major', axis = 'both')
    ax1.set_ylabel('Z (m)')
    ax1.legend(loc = 'upper left')
    titleStr = 'Shifted ICESat-2 and Reference Data (' + str(rasterResolution) + ' m Raster)'
    ax1.set_title(titleStr)
    
    # Make Z Error vs Y subplot 2
    ax2.scatter(truthRasterYCommonFinal/1000, zErrorCommonFinal, color = myBlue, s=0.7, label = 'Z Error')
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
    
    # Close figure
    plt.close()
    
    
# Function to make Z vs T plot
def plotZT(measRasterTCommonFinal, measRasterZCommonFinal, zErrorCommonFinal, atlMeasuredData, atlCorrections, rasterResolution, showPlots, outFilePath, counter):
    
    # Define plot colors
    myBlue = (0/255, 114/255, 178/255)
    myOrange = (213/255, 94/255, 0/255)
    
    # Make Z vs Y subplot 1
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.scatter(measRasterTCommonFinal, measRasterZCommonFinal, color = myOrange, s=0.7)
    ax1.axis('tight')
    ax1.grid(b = True, which = 'major', axis = 'both')
    ax1.set_ylabel('Z (m)')
    titleStr = 'IceSat-2 Time Data (' + str(rasterResolution) + ' m Raster)'
    ax1.set_title(titleStr)
    
    # Make Z Error vs Y subplot 2
    ax2.scatter(measRasterTCommonFinal, zErrorCommonFinal, color = myBlue, s=0.7)
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
    
    # Close figure
    plt.close()
    
    
# Lasso Selector Tool    
class SelectFromCollection(object):
    
    def __init__(self, ax, collection, alpha_other=0.3):
        

#        xys = []
#        Npts = []
#        fc = []
#        ec = []
#        ind = []
#        
#        for i in range(0,len(collection)):
#            
#            xysSingle = collection[i].get_offsets().data
#            NptsSingle = len(xysSingle)
#            fcSingle = np.tile(collection[i].get_facecolors(),(NptsSingle,1))
#            ecSingle = np.tile(collection[i].get_edgecolors(),(NptsSingle,1))
#            indSingle = []
#            
#            xys.append(xysSingle)
#            Npts.append(NptsSingle)
#            fc.append(fcSingle)
#            ec.append(ecSingle)
#            ind.append(indSingle)
#            
#        # endFor
#        
#        self.xys = xys
#        self.Npts = Npts
#        self.fc = fc
#        self.ec = ec
#        self.ind = ind
        
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.ax = ax
        self.alpha_other = alpha_other
        self.xys = collection.get_offsets().data
        self.Npts = len(self.xys)
        self.baseColor = collection.get_facecolors()
        self.fc = []
        self.fc = np.tile(collection.get_facecolors(),(self.Npts,1))
        self.ind = []
        self.lasso = LassoSelector(ax, onselect=self.onselect)       
        self.active = True
        

    def onselect(self, verts):
                
        # Get vertices of lasso
        path = Path(verts)
        
#        for i in range(0,len(self.collection)): 
#            self.fc[i] = self.ec[i]
#            self.ind[i] = np.nonzero(path.contains_points(self.xys[i]))[0]
#            self.fc[i][self.ind[i],0] = 1
#            self.fc[i][self.ind[i],1] = 0
#            self.fc[i][self.ind[i],2] = 0
#            self.collection[i].set_edgecolors(self.fc[i][:,:])
#        # endFor

        self.ind = []
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        fc = self.fc
        fc[:,:] = self.baseColor
        fc[self.ind,0] = 1
        fc[self.ind,1] = 0
        fc[self.ind,2] = 0
        self.collection.set_edgecolors(fc)
        self.canvas.blit(self.ax.bbox)
        self.active = True
        
        
    def disconnect(self):
        self.lasso.disconnect_events()
        
#        for i in range(0,len(self.collection)): 
#            self.collection[i].set_edgecolors(self.ec[i][:,:])
#        # endFor

        self.ind = []
        self.collection.set_edgecolors(self.baseColor)        
        self.canvas.blit(self.ax.bbox)
        self.active = False
        
        
    def revert(self):
        
#        for i in range(0,len(self.collection)): 
#            self.collection[i].set_edgecolors(self.ec[i][:,:])
#        # endFor
        
        self.ind = []
        self.collection.set_edgecolors(self.baseColor)
        self.canvas.blit(self.ax.bbox)
        self.active = True
        
 
def getPlotPts(self, fig, ax, pts, outPath, origTitle, atl03Data, onState):  
    
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    # endif
    
    # Call point selector    
    selector = SelectFromCollection(ax, pts)
        
    if(onState):
    
        # Set directions for user in figure title
        ax.set_title('Left-Click to Select Points, "Enter" to Accept, "Esc" to Exit')
    
        # Keyboard Button callback
        def accept(event):
            
    #        # Get selected points
    #        selected_time = []
    #        selected_lat = []
    #        selected_lon = []
    #        selected_easting = []
    #        selected_northing = []
    #        selected_crossTrack = []
    #        selected_alongTrack = []
    #        selected_z = []
    #        selected_classification = []
    #        selected_signalConf = []
    #        
    #        for i in range(0,len(selector.ind)):
    #            
    #            selected_time = np.append(selected_time,atl03Data.time[selector.ind[i]])
    #            selected_lat = np.append(selected_lat,atl03Data.lat[selector.ind[i]])
    #            selected_lon = np.append(selected_lon,atl03Data.lon[selector.ind[i]])
    #            selected_easting = np.append(selected_easting,atl03Data.easting[selector.ind[i]])
    #            selected_northing = np.append(selected_northing,atl03Data.northing[selector.ind[i]])
    #            selected_crossTrack = np.append(selected_crossTrack,atl03Data.crossTrack[selector.ind[i]])
    #            selected_alongTrack = np.append(selected_alongTrack,atl03Data.alongTrack[selector.ind[i]])
    #            selected_z = np.append(selected_z,atl03Data.z[selector.ind[i]])
    #            selected_classification = np.append(selected_classification,atl03Data.classification[selector.ind[i]])
    #            selected_signalConf = np.append(selected_signalConf,atl03Data.signalConf[selector.ind[i]])
    #            
    #        # endFor
                
            # 'Enter' button callback
            if(event.key == 'enter'):
                
                # Get selected points
                selected_time = atl03Data.time[selector.ind]
                selected_lat = atl03Data.lat[selector.ind]
                selected_lon = atl03Data.lon[selector.ind]
                selected_easting = atl03Data.easting[selector.ind]
                selected_northing = atl03Data.northing[selector.ind]
                selected_crossTrack = atl03Data.crossTrack[selector.ind]
                selected_alongTrack = atl03Data.alongTrack[selector.ind]
                selected_z = atl03Data.z[selector.ind]
                selected_classification = atl03Data.classification[selector.ind]
                selected_signalConf = atl03Data.signalConf[selector.ind]
                
                # Get total number of selected points
                numPts = len(selected_time)
    
                if(numPts>0):
        
                    # Update directions in title
                    ax.set_title('%s Points Saved\nLeft-Click to Select Again, "Enter" to Accept, "Esc" to Exit' % numPts)
    
                    # Get output file name               
                    outName = 'selectedPoints*.csv'
                    
                    # Get .csv headers
                    if(atl03Data.zone=='3413' or atl03Data.zone=='3976'):
                    
                        namelist = ['Time (sec)', 'Latitude (deg)', 'Longitude (deg)', \
                                'Polar Stereo X (m)', 'Polar Stereo Y (m)', \
                                'Cross-Track (m)', 'Along-Track (m)', \
                                'Height (m)', \
                                'Classification', 'Signal Confidence']
                    else:
                    
                        namelist = ['Time (sec)', 'Latitude (deg)', 'Longitude (deg)', \
                                'Easting (m)', 'Northing (m)', \
                                'Cross-Track (m)', 'Along-Track (m)', \
                                'Height (m)', \
                                'Classification', 'Signal Confidence']
                    # endIf
                    
                    # Get .csv data
                    datalist = [selected_time, selected_lat, selected_lon, \
                                selected_easting, selected_northing, \
                                selected_crossTrack, selected_alongTrack,\
                                selected_z, \
                                selected_classification, selected_signalConf] 
                
                    # Get output file name (append _N.txt for file number)
                    outFilePath = os.path.normpath(outPath + '/' + outName)
                    matchingOutFiles = glob.glob(outFilePath)
                    if(matchingOutFiles):
                        matchingNums = np.array([])
                        for numFile in range(0,len(matchingOutFiles)):
                            matchingFile = os.path.basename(matchingOutFiles[numFile])
                            matchingNum = matchingFile.split('_')[2]
                            matchingNums = np.append(matchingNums,matchingNum)
                        # endFor
                        matchingNumMax = max(matchingNums)
                        matchingNumNew = str(int(matchingNumMax) + 1)
                        outFile = os.path.normpath(outPath + '/selectedPoints_file_' + matchingNumNew + '_pts_' + str(numPts) + '.csv')
                    else:
                        outFile = os.path.normpath(outPath + '/selectedPoints_file_1_pts_' + str(numPts) + '.csv')
                    # endIf
    
                    # Write output data to .csv file
                    writeArrayToCSV(outFile, namelist, datalist)
                    
                    # Revert to original color scheme
                    selector.revert()
    
                else:
                    
                    # Update directions in title
                    ax.set_title('No Points Saved\nLeft-Click to Select Again, "Enter" to Accept, "Esc" to Exit')
    
                # endIf
            # endIf
            
            # 'Escape' button callback
            if(event.key == 'escape'):
                selector.disconnect()
                ax.set_title(origTitle)
                
                self.disable()
                
            # endIf
            
            # Update figure
            fig.canvas.draw_idle()
            
        # endDef
            
        # Mouse-Button Release Callback
        def release(event):
            
            if(selector.active):
            
    #            # Get number of selected points
    #            numPts = 0
    #            for i in range(0,len(selector.ind)):
    #                numPts += len(selector.ind[i])
    #            # endFor
                            
                # Get number of selected points
                numPts = len(selector.ind)
                
                # Update figure/title
                ax.set_title('%s Points Selected\nLeft-Click to Select Again, "Enter" to Accept, "Esc" to Exit' % numPts)
                fig.canvas.draw_idle()
                        
            # endIf
        # endDef
            
            
        # Create callback when keyboard buttons are pressed
        fig.canvas.mpl_connect('key_press_event', accept)
        
        # Create callback when mouse-button is released
        fig.canvas.mpl_connect('button_release_event', release)
    
    else:
    
        selector.disconnect()
        ax.set_title(origTitle)
        fig.canvas.draw_idle()
        
    # endIf
            
# endDef
    

class getPoints(ToolToggleBase):
        
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    # endif
    
    default_toggled = False
    
    def __init__(self, *args, figNum, ax, pts, outPath, origTitle, atl03Data, **kwargs):
        
        # Initialize base variables
        self.figNum = figNum
        self.ax = ax
        self.pts = pts
        self.outPath = outPath
        self.origTitle = origTitle
        self.atl03Data = atl03Data
        
        # Call super-initialized variables
        super().__init__(*args, **kwargs)
        
        
    def enable(self, *args, **kwargs):
        
        # Get inputs for getPlotPts function
        figNum = self.figNum
        ax = self.ax
        pts = self.pts
        outPath = self.outPath 
        origTitle = self.origTitle
        atl03Data = self.atl03Data
        onState = True
        
        # Call getPlotPts function
        getPlotPts(self, figNum, ax, pts, outPath, origTitle, atl03Data, onState)
        
        
    def disable(self, *args, **kwargs):
        
        # Get inputs for getPlotPts function
        figNum = self.figNum
        ax = self.ax
        pts = self.pts
        outPath = self.outPath 
        origTitle = self.origTitle
        atl03Data = self.atl03Data
        onState = False
        
        # Call getPlotPts function
        getPlotPts(self, figNum, ax, pts, outPath, origTitle, atl03Data, onState)
    
    
# Function to plot any x,y data from getAtlMeasuredSwath_auto.py
def getPlot(xData, yData, xLabel, yLabel, title, outPath, origTitle, atl03Data, filterType = [], filterData = [], filterNum = []):
    
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    # endif
        
    # Open figure window
    fig1 = plt.figure()
    
    if(np.size(filterData) > 0):
        
        # Loop through filter numbers
        for i in range(0,np.size(filterNum)):
            
            # Filter data
            matchVal = filterNum[i]
            matchingInds = filterData == matchVal
            if(len(matchingInds>0)):
                
                xDataPlot = xData[matchingInds]
                yDataPlot = yData[matchingInds]
                            
                if('class' in filterType.lower()):
                    if(matchVal < 1):
                        matchStr = 'ATL03 Unclassified'
                        myColor = [194/255, 197/255, 204/255]
                    elif(matchVal==1):
                        matchStr = 'ATL03 Ground'
                        myColor = [210/255, 184/255, 38/255]
                    elif(matchVal==2):
                        matchStr = 'ATL03 Canopy'
                        myColor = [69/255, 129/255, 26/255]
                    elif(matchVal==3):
                        matchStr = 'ATL03 Top of Canopy'
                        myColor = [133/255, 243/255, 52/255]
                    # endIf
                else:
                    if(matchVal==0):
                        myColor = [194/255, 197/255, 204/255]
                    elif(matchVal==1):
                        myColor = [0, 0, 0]
                    elif(matchVal==2):
                        myColor = [69/255, 129/255, 26/255]
                    elif(matchVal==3):
                        myColor = [1, 0, 0]
                    elif(matchVal==4):
                        myColor = [0, 0, 1]
                    # endIf
                    
                    matchStr = 'ATL03 Sig Conf = ' + str(matchVal)
                    
                # endIf
                            
                # Plot filtered data
                plt.scatter(xDataPlot, yDataPlot, color=myColor, label=matchStr, s=0.7)
            
            # endIf
        
        # endFor
        
        plt.legend(loc = 'upper left')
    
    else:
        
        # Plot all data
        plt.scatter(xData, yData, color=[194/255, 197/255, 204/255], s=0.7, label='ATL03 Data')
        plt.legend(loc = 'upper left')
        
    # endIf
    
    # Figure properties
    plt.axis('tight')
    plt.grid(b = True, which = 'major', axis = 'both')
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(title)
    
    # Get figure handles
    figNum = plt.gcf()
    ax = plt.gca()
    pts = ax.collections[-1]
    
    # Route 'Select Points' button to matplotlib figure window
    fig1.canvas.manager.toolmanager.add_tool('Select Points', getPoints, figNum = figNum, ax = ax, pts = pts, outPath = outPath, origTitle = origTitle, atl03Data = atl03Data)
    fig1.canvas.manager.toolbar.add_tool('Select Points', 'navigation', 3)
    
    # Show plot
    fig1.show()
    
    
# Function to plot ATL08 data from getAtlMeasuredSwath_auto.py
def getPlot_atl08(xData, yData, xLabel, yLabel, title, yName):
    
#    # Define plot colors
#    myColors = np.array([[0/255, 114/255, 178/255],[0/255, 158/255, 115/255],[213/255, 94/255, 0/255]])

    # Open figure window
    figNum = plt.gcf().number
    fig1 = plt.figure(figNum)
    
    if('canopy' in yName.lower()):
        matchStr = 'ATL08 Max Canopy'
        myColor = [69/255, 129/255, 26/255]
    elif('bestfit' in yName.lower()):
        matchStr = 'ATL08 Terrain Best Fit'
        myColor = [210/255, 184/255, 38/255]
    elif('median' in yName.lower()): 
        matchStr = 'ATL08 Terrain Median'
        myColor = [148/255, 0/255, 211/255]
    # endIf

    # Plot all data
    plt.plot(xData, yData, 'o', markerfacecolor=myColor, markeredgecolor='k', label=matchStr)
    
    # Figure properties
    plt.legend(loc = 'upper left')
    plt.axis('tight')
    plt.grid(b = True, which = 'major', axis = 'both')
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(title)
    plt.draw()
    fig1.show()


# Function to plot Truth data from getAtlTruthSwath_auto.py
def getPlot_truth(xData, yData, xLabel, yLabel, title, yName):
    
    # Define plot colors
    myColor = [0.3, 0.3, 0.3]

    # Open figure window
    figNum = plt.gcf().number
    fig1 = plt.figure(figNum)

    # Plot all data
    plt.scatter(xData, yData, color=myColor, label=yName, s=1.0, zorder=0)
    
    # Figure properties
    plt.legend(loc = 'upper left')
    plt.axis('tight')
    plt.grid(b = True, which = 'major', axis = 'both')
    plt.draw()
    fig1.show()
    
    
# Function to plot Corrected Measured data from getMeasurementError_auto.py
def getPlot_measCorr(xData, yData, xLabel, yLabel, title, yName):
    
    # Define plot colors
    myColor = [1.0, 0.65, 0.0]

    # Open figure window
    figNum = plt.gcf().number
    fig1 = plt.figure(figNum)

    # Plot all data
    plt.scatter(xData, yData, color=myColor, label=yName, s=0.7)
    
    # Figure properties
    plt.legend(loc = 'upper left')
    plt.axis('tight')
    plt.grid(b = True, which = 'major', axis = 'both')
    plt.draw()
    fig1.show()