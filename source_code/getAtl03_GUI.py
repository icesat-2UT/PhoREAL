# -*- coding: utf-8 -*-
"""
This script loads the PhoREAL GUI which executes getAtlMeasuredSwath_auto.py

Copyright 2019 Applied Research Laboratories, University of Texas at Austin

This package is free software; the copyright holder gives unlimited
permission to copy and/or distribute, with or without modification, as
long as this notice is preserved.

Authors:
    Mike Alonzo
    Eric Guenther
    
Date: September 20, 2019

"""

import tkinter as tk
from tkinter import filedialog
from tkinter.ttk import Combobox
from tkinter.ttk import Progressbar
import numpy as np
from getAtlMeasuredSwath_auto import getAtlMeasuredSwath
from icesatPlot import getPlot
from gui_logo import images


# Print opening message
print('\n')
print('***********************************************')
print('ATL03 GUI is opening, this may take a moment...')
print('***********************************************')
print('\n')


window = tk.Tk()

# GUI title
window.title('PhoREAL v1.0 - Applied Research Labs (The University of Texas at Austin)')

# GUI size
window.geometry('540x910')

# Set icon image for GUI
window.wm_iconbitmap(images)


### ATL03 INPUT FILE

# ATL03 File Text Label
lbl = tk.Label(window, text='ATL03 File:', font=('Arial Bold', 16), anchor = 'w', justify='left')
lbl.place(x=10, y=15)

# ATL03 File Text Entry Box
atl03_textBox = tk.Entry(window, width=70)
atl03_textBox.place(x=10, y=50)

# Browse ATL03 File Button Callback
def browseAtl03():
    currentData = atl03_textBox.get()
    atl03File = filedialog.askopenfilename(title = 'Select ATL03 .h5 file')
    atl03_textBox.delete(0,len(currentData))
    atl03_textBox.insert(0,atl03File) 
# endDef
    
# ATL03 File Browse Button
atl03BrowseButton = tk.Button(window, text='Browse', font=('Arial Bold', 12), command=browseAtl03) 
atl03BrowseButton.place(x=450, y=40)


### ATL08 INPUT FILE

# ATL08 File Text Label
lbl = tk.Label(window, text='ATL08 File (Optional):', font=('Arial Bold', 16), anchor = 'w', justify='left')
lbl.place(x=10, y=85)

# ATL08 File Text Entry Box
atl08_textBox = tk.Entry(window, width=70) 
atl08_textBox.place(x=10, y=120)

# Browse ATL08 File Button Callback
def browseAtl08():
    currentData = atl08_textBox.get()
    atl08File = filedialog.askopenfilename(title = 'Select ATL08 .h5 file')
    atl08_textBox.delete(0,len(currentData))
    atl08_textBox.insert(0,atl08File) 
    # endIf
# endDef
    
# ATL03 File Browse Button
atl08BrowseButton = tk.Button(window, text='Browse', font=('Arial Bold', 12), command=browseAtl08) 
atl08BrowseButton.place(x=450, y=110)


### Output Path

# Output Path Text Label
lbl = tk.Label(window, text='Output Directory:', font=('Arial Bold', 16), anchor = 'w', justify='left')
lbl.place(x=10, y=155)

# Output Path Text Entry
outPath_textBox = tk.Entry(window, width=70) 
outPath_textBox.place(x=10, y=190)

# Browse Output File Button Callback
def browseOutput():
    currentData = outPath_textBox.get()
    outPath = filedialog.askdirectory(title = 'Select Output File Path')
    outPath_textBox.delete(0,len(currentData))
    outPath_textBox.insert(0,outPath) 
# endDef
    
# ATL03 File Browse Button
outPathBrowseButton = tk.Button(window, text='Browse', font=('Arial Bold', 12), command=browseOutput) 
outPathBrowseButton.place(x=450, y=180)


### GT Number Text Label
lbl = tk.Label(window, text='Ground Track Numbers:', font=('Arial Bold', 16), anchor = 'w', justify='left')
lbl.place(x=10, y=230)

gtNumsAll = np.array(['GT1R','GT1L','GT2R','GT2L','GT3R','GT3L'])

# GT Number Checkboxes
gtNum1rChkState = tk.BooleanVar()
gtNum1rChkState.set(True)
gtNum1r_checkBox = tk.Checkbutton(window, text = 'GT1R', font=('Arial Bold', 14), var = gtNum1rChkState) 
gtNum1r_checkBox.place(x=30, y=260)

gtNum2rChkState = tk.BooleanVar()
gtNum2rChkState.set(True)
gtNum2r_checkBox = tk.Checkbutton(window, text = 'GT2R', font=('Arial Bold', 14), var = gtNum2rChkState) 
gtNum2r_checkBox.place(x=30, y=290)

gtNum3rChkState = tk.BooleanVar()
gtNum3rChkState.set(True)
gtNum3r_checkBox = tk.Checkbutton(window, text = 'GT3R', font=('Arial Bold', 14), var = gtNum3rChkState) 
gtNum3r_checkBox.place(x=30, y=320)

gtNum1lChkState = tk.BooleanVar()
gtNum1lChkState.set(False)
gtNum1l_checkBox = tk.Checkbutton(window, text = 'GT1L', font=('Arial Bold', 14), var = gtNum1lChkState) 
gtNum1l_checkBox.place(x=150, y=260)

gtNum2lChkState = tk.BooleanVar()
gtNum2lChkState.set(False)
gtNum2l_checkBox = tk.Checkbutton(window, text = 'GT2L', font=('Arial Bold', 14), var = gtNum2lChkState) 
gtNum2l_checkBox.place(x=150, y=290)

gtNum3lChkState = tk.BooleanVar()
gtNum3lChkState.set(False)
gtNum3l_checkBox = tk.Checkbutton(window, text = 'GT3L', font=('Arial Bold', 14), var = gtNum3lChkState) 
gtNum3l_checkBox.place(x=150, y=320)


### Check None Button Callback
def checkNone():
    if(trimNoneModeChkState):
        trimAutoModeChkState.set(False)
        trimManualModeChkState.set(False)
        latModeChkState.set(False)
        latMode_checkBox.config(state = 'disabled')
        latMin_textBox.config(state = 'disabled')
        latMax_textBox.config(state = 'disabled')
        timeModeChkState.set(False)
        timeMode_checkBox.config(state = 'disabled')
        timeMin_textBox.config(state = 'disabled')
        timeMax_textBox.config(state = 'disabled')
    # endIf
# endDef
    
### Check Auto Button Callback
def checkAuto():
    if(trimAutoModeChkState):
        trimNoneModeChkState.set(False)
        trimManualModeChkState.set(False)
        latModeChkState.set(False)
        latMode_checkBox.config(state = 'disabled')
        latMin_textBox.config(state = 'disabled')
        latMax_textBox.config(state = 'disabled')
        timeModeChkState.set(False)
        timeMode_checkBox.config(state = 'disabled')
        timeMin_textBox.config(state = 'disabled')
        timeMax_textBox.config(state = 'disabled')
    # endIf
# endDef
    
### Check Manual Button Callback
def checkManual():
    if(trimManualModeChkState):
        trimNoneModeChkState.set(False)
        trimAutoModeChkState.set(False)
        latModeChkState.set(False)
        latMode_checkBox.config(state = 'normal')
        latMin_textBox.config(state = 'normal')
        latMax_textBox.config(state = 'normal')
        timeModeChkState.set(False)
        timeMode_checkBox.config(state = 'normal')
        timeMin_textBox.config(state = 'normal')
        timeMax_textBox.config(state = 'normal')
    # endIf
# endDef
        

### Trim Info Text Label
lbl = tk.Label(window, text='Trim Info Options:', font=('Arial Bold', 16))
lbl.place(x=10, y=370)

# Trim Info Checkboxes
trimNoneModeChkState = tk.BooleanVar()
trimNoneModeChkState.set(True)
trimNoneMode_checkBox = tk.Checkbutton(window, text = 'None', font=('Arial Bold', 14), var = trimNoneModeChkState, command = checkNone) 
trimNoneMode_checkBox.place(x=30, y=400)

trimAutoModeChkState = tk.BooleanVar()
trimAutoModeChkState.set(False)
trimAutoMode_checkBox = tk.Checkbutton(window, text = 'Auto', font=('Arial Bold', 14), var = trimAutoModeChkState, command = checkAuto) 
trimAutoMode_checkBox.place(x=30, y=430)

trimManualModeChkState = tk.BooleanVar()
trimManualModeChkState.set(False)
trimManualMode_checkBox = tk.Checkbutton(window, text = 'Manual', font=('Arial Bold', 14), var = trimManualModeChkState, command = checkManual) 
trimManualMode_checkBox.place(x=30, y=460)


### Check Latitude Button Callback
def checkLat():
    if(latModeChkState):
        timeModeChkState.set(False)
    # endIf
# endDef
    
### Check Time Button Callback
def checkTime():
    if(timeModeChkState):
        latModeChkState.set(False)
    # endIf
# endDef
    
    
### Lat Mode Checkbox
latModeChkState = tk.BooleanVar()
latModeChkState.set(False)
latMode_checkBox = tk.Checkbutton(window, text = 'Latitude', font=('Arial Bold', 14), var = latModeChkState, command = checkLat, state = 'disabled') 
latMode_checkBox.place(x=60, y=490)

### Lat Min Text Label
lbl = tk.Label(window, text='Min:', font=('Arial Bold', 14))
lbl.place(x=200, y=495)

# Lat Min Text Entry
latMin_textBox = tk.Entry(window, width=10, state = 'disabled') 
latMin_textBox.place(x=250, y=500)

### Lat Max Text Label
lbl = tk.Label(window, text='Max:', font=('Arial Bold', 14))
lbl.place(x=340, y=495)

# Lat Max Text Entry
latMax_textBox = tk.Entry(window, width=10, state = 'disabled') 
latMax_textBox.place(x=400, y=500)

### Lat Degrees Text Label
lbl = tk.Label(window, text='Degrees', font=('Arial Bold', 10))
lbl.place(x=460, y=500)

### Time Mode Checkbox
timeModeChkState = tk.BooleanVar()
timeModeChkState.set(False)
timeMode_checkBox = tk.Checkbutton(window, text = 'Time', font=('Arial Bold', 14), var = timeModeChkState, command = checkTime, state = 'disabled') 
timeMode_checkBox.place(x=60, y=520)

### Time Min Text Label
lbl = tk.Label(window, text='Min:', font=('Arial Bold', 16))
lbl.place(x=200, y=525)

# Time Min Text Entry
timeMin_textBox = tk.Entry(window, width=10, state = 'disabled') 
timeMin_textBox.place(x=250, y=530)

### Time Max Text Label
lbl = tk.Label(window, text='Max:', font=('Arial Bold', 16))
lbl.place(x=340, y=525)

# Time Max Text Entry
timeMax_textBox = tk.Entry(window, width=10, state = 'disabled') 
timeMax_textBox.place(x=400, y=530)

### Time Seconds Text Label
lbl = tk.Label(window, text='Seconds', font=('Arial Bold', 10))
lbl.place(x=460, y=530)


### Create .las Text Entry
createLasChkState = tk.BooleanVar()
createLasChkState.set(True)
createLas_checkBox = tk.Checkbutton(window, text = 'Create .las File', font=('Arial Bold', 16), var = createLasChkState) 
createLas_checkBox.place(x=10, y=565)


### Create .kml Text Entry
createKmlChkState = tk.BooleanVar()
createKmlChkState.set(True)
createKml_checkBox = tk.Checkbutton(window, text = 'Create .kml File', font=('Arial Bold', 16), var = createKmlChkState) 
createKml_checkBox.place(x=10, y=595)


### Create .csv Text Entry
createCsvChkState = tk.BooleanVar()
createCsvChkState.set(True)
createCsv_checkBox = tk.Checkbutton(window, text = 'Create .csv File', font=('Arial Bold', 16), var = createCsvChkState) 
createCsv_checkBox.place(x=10, y=625)


### Make status bar
lbl = tk.Label(window, text='Progress:', font=('Arial Bold', 10))
lbl.place(x=250, y=635)
statusBar = Progressbar(window, length=190)
statusBar['value'] = 0
statusBar.place(x=320, y=635)
    

# Run Button Callback
def runAtl03():
    
    # Update status bar
    statusBar['value'] = 0
    window.update()
        
    # Make atlMeasuredData a global variable
    global atlMeasuredData
        
    # Initialize variables
    atlMeasuredData = []
    atlMeasuredDataSingle = []
          
    # Try code
    try:
        
        # Get inputs
        atl03FilePath = atl03_textBox.get().strip()
        atl08FilePath = atl08_textBox.get().strip()
        outFilePath = outPath_textBox.get().strip()
        gtNumsTF = [gtNum1rChkState.get(), gtNum1lChkState.get(), gtNum2rChkState.get(), gtNum2lChkState.get(), gtNum3rChkState.get(), gtNum3lChkState.get()]    
        gtNums = gtNumsAll[gtNumsTF]    
    
        if(len(atl08FilePath)>0):
            atl08FileExists = True
        else:
            atl08FileExists = False
        # endIf
        
        # Get trim inputs
        if(trimNoneModeChkState.get()):
            trimInfo = 'none'
        elif(trimAutoModeChkState.get()):
            trimInfo = 'auto'
        elif(trimManualModeChkState.get()):
            trimMode = 'manual'
            if(latModeChkState.get()):
                trimType = 'lat'
                trimMin = latMin_textBox.get()
                trimMax = latMax_textBox.get()
            elif(timeModeChkState.get()):
                trimType = 'time'
                trimMin = timeMin_textBox.get()
                trimMax = timeMax_textBox.get()
            trimInfo = trimMode + ',' + trimType + ',' + trimMin + ',' + trimMax
            
        # Get output file options
        createLasFile = createLasChkState.get()
        createKmlFile = createKmlChkState.get()
        createCsvFile = createCsvChkState.get()
            
        # Update status bar
        statusBar['value'] = 25
        window.update()
        
        # Loop through all gt nums
        headerData = False
        rotationData = False
        for i in range(0,len(gtNums)):
        
            # Get GT Nums
            gtNum = gtNums[i].lower()
                
            # Run getAtlMeasuredSwath function
            atlMeasuredDataSingle, headerData, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, gtNum, trimInfo, createLasFile, createKmlFile, createCsvFile)
                        
            # Append measured data for each run
            atlMeasuredData.append(atlMeasuredDataSingle)
                        
            # Update status bar
            progressNum = ( (i+1) / len(gtNums) ) * 75
            statusBar['value'] = 25 + progressNum
            window.update()
                        
        # endFor
                
        if(atlMeasuredData):
            
            # Set Ground Tracks to plot
            gtNumsTuple = tuple(gtNums)
            gtNumPlotBox['values'] = gtNumsTuple
            gtNumPlotBox.current(0)
            
            # Set X Vals to plot
            xValsTuple = ('Time (sec)', 'Latitude (deg)', 'Longitude (deg)', \
                            'UTM Easting (m)', 'UTM Northing (m)', \
                            'Cross-Track (m)', 'Along-Track (m)', \
                            'Height (m)', \
                            'Classification', 'Signal Confidence')
            xValsBox['values'] = xValsTuple
            xValsBox.current(0)
            
            # Set Y Vals to plot
            yValsTuple = ('Time (sec)', 'Latitude (deg)', 'Longitude (deg)', \
                            'UTM Easting (m)', 'UTM Northing (m)', \
                            'Cross-Track (m)', 'Along-Track (m)', \
                            'Height (m)', \
                            'Classification', 'Signal Confidence')
            yValsBox['values'] = yValsTuple
            yValsBox.current(7)
            
            # Set X Label
            xAxisVal = xValsBox.get()
            currentData = xlabel_textBox.get()
            xlabel_textBox.delete(0,len(currentData))
            xlabel_textBox.insert(0,xAxisVal) 
            
            # Set Y label
            yAxisVal = yValsBox.get()
            currentData = ylabel_textBox.get()
            ylabel_textBox.delete(0,len(currentData))
            ylabel_textBox.insert(0,yAxisVal) 
            
            # Set Title label
            gtNumVal = gtNumPlotBox.get()
            currentData = title_textBox.get()
            title_textBox.delete(0,len(currentData))
            title_textBox.insert(0,gtNumVal) 
            
            # Set Vals to filter on
            if(atl08FileExists):
                filterTuple = ('  ','Classification', 'Signal Confidence')
            else:
                filterTuple = ('  ', 'Signal Confidence')
            # endIf
            filterBox['values'] = filterTuple
            filterBox.current(0)
            
            # Set Filter Number Checkboxes
            filter0_checkBox.place_forget()
            filter1_checkBox.place_forget()
            filter2_checkBox.place_forget()  
            filter3_checkBox.place_forget()        
            filter4_checkBox.place_forget()
            filter0ChkState.set(False)
            filter1ChkState.set(False)
            filter2ChkState.set(False)
            filter3ChkState.set(False)
            filter4ChkState.set(False)
            
        # endIf
        
        # Print completion message
        print('RUN COMPLETE.')
    
    except:
        print('Could not get ATL03 data. Please check inputs.')
    #endTry
        
# endDef
    
# Run Button
btn = tk.Button(window, text='RUN', font=('Arial Bold', 20), width = 15, command=runAtl03) 
btn.place(x=250, y=575)


### Plotting Section

# Plot text
lbl = tk.Label(window, text='Plotting Options:', font=('Arial Bold', 16), anchor = 'w', justify='left')
lbl.place(x=10, y=675)

# X Axis Callback
def xAxisCallback(event):
    xAxisVal = xValsBox.get()
    currentData = xlabel_textBox.get()
    xlabel_textBox.delete(0,len(currentData))
    xlabel_textBox.insert(0,xAxisVal)   
# endDef
    
# Y Axis Callback
def yAxisCallback(event):
    yAxisVal = yValsBox.get()
    currentData = ylabel_textBox.get()
    ylabel_textBox.delete(0,len(currentData))
    ylabel_textBox.insert(0,yAxisVal)   
# endDef
    
# GT Num Plot Callback
def gtNumPlotCallback(event):
    gtNumVal = gtNumPlotBox.get()
    currentData = title_textBox.get()
    title_textBox.delete(0,len(currentData))
    title_textBox.insert(0,gtNumVal)   
# endDef
  
    
# Filter Number Checkboxes
filter0ChkState = tk.BooleanVar()
filter0ChkState.set(False)
filter0_checkBox = tk.Checkbutton(window, text = '0', font=('Arial Bold', 12), var = filter0ChkState) 

filter1ChkState = tk.BooleanVar()
filter1ChkState.set(False)
filter1_checkBox = tk.Checkbutton(window, text = '1', font=('Arial Bold', 12), var = filter1ChkState) 

filter2ChkState = tk.BooleanVar()
filter2ChkState.set(False)
filter2_checkBox = tk.Checkbutton(window, text = '2', font=('Arial Bold', 12), var = filter2ChkState) 

filter3ChkState = tk.BooleanVar()
filter3ChkState.set(False)
filter3_checkBox = tk.Checkbutton(window, text = '3', font=('Arial Bold', 12), var = filter3ChkState) 

filter4ChkState = tk.BooleanVar()
filter4ChkState.set(False)
filter4_checkBox = tk.Checkbutton(window, text = '4', font=('Arial Bold', 12), var = filter4ChkState) 
     
# Filter Choice Callback
def filterChoiceCallback(event):
    
    # Get  filter value
    filterChoice = filterBox.get()
    
    if('  ' in filterChoice.lower()):
        
        # Filter Number Checkboxes
        filter0_checkBox.place_forget()
        filter1_checkBox.place_forget()
        filter2_checkBox.place_forget()  
        filter3_checkBox.place_forget()        
        filter4_checkBox.place_forget()
        
    elif('class' in filterChoice.lower()):
        
        # Filter Number Checkboxes
        filter0_checkBox.place(x=30, y=870) 
        filter0_checkBox.config(text = 'Unclassified')
        filter0ChkState.set(False)
        filter1_checkBox.place(x=160, y=870)
        filter1_checkBox.config(text = 'Ground')
        filter1ChkState.set(False)
        filter2_checkBox.place(x=250, y=870)
        filter2_checkBox.config(text = 'Low Veg') 
        filter2ChkState.set(False)
        filter3_checkBox.place(x=350, y=870)  
        filter3_checkBox.config(text = 'High Veg') 
        filter3ChkState.set(False)
        filter4_checkBox.place_forget()
        filter4ChkState.set(False)
        
    elif('signal' in filterChoice.lower()):
        
        # Filter Number Checkboxes
        filter0_checkBox.place(x=30, y=870) 
        filter0_checkBox.config(text = '0')
        filter0ChkState.set(False)
        filter1_checkBox.place(x=80, y=870)  
        filter1_checkBox.config(text = '1')
        filter1ChkState.set(False)
        filter2_checkBox.place(x=130, y=870)
        filter2_checkBox.config(text = '2') 
        filter2ChkState.set(False)        
        filter3_checkBox.place(x=180, y=870)
        filter3_checkBox.config(text = '3')
        filter3ChkState.set(False)
        filter4_checkBox.place(x=230, y=870)
        filter4_checkBox.config(text = '4')
        filter4ChkState.set(False)
            
    # endIf
# endDef
    
# GT Num Plot Combo Box
lbl = tk.Label(window, text='Ground Track:', font=('Arial Bold', 14), anchor = 'w', justify='left')
lbl.place(x=30, y=715)
gtNumPlotBox = Combobox(window, width=10)
gtNumPlotBox.place(x= 180, y = 720)
gtNumPlotBox.bind("<<ComboboxSelected>>", gtNumPlotCallback)

# X Axis Combo Box
lbl = tk.Label(window, text='X Axis:', font=('Arial Bold', 14), anchor = 'w', justify='left')
lbl.place(x=30, y=745)
xValsBox = Combobox(window)
xValsBox.place(x= 110, y = 750)
xValsBox.bind("<<ComboboxSelected>>", xAxisCallback)

# X Label Entry Box
lbl = tk.Label(window, text='X Label:', font=('Arial Bold', 14), anchor = 'w', justify='left')
lbl.place(x=280, y=745)
xlabel_textBox = tk.Entry(window, width=20)
xlabel_textBox.place(x=370, y=750)

# Y Axis Combo Box
lbl = tk.Label(window, text='Y Axis:', font=('Arial Bold', 14), anchor = 'w', justify='left')
lbl.place(x=30, y=775)
yValsBox = Combobox(window)
yValsBox.place(x= 110, y = 780)
yValsBox.bind("<<ComboboxSelected>>", yAxisCallback)

# Y Label Entry Box
lbl = tk.Label(window, text='Y Label:', font=('Arial Bold', 14), anchor = 'w', justify='left')
lbl.place(x=280, y=775)
ylabel_textBox = tk.Entry(window, width=20)
ylabel_textBox.place(x=370, y=780)

# Title Entry Box
lbl = tk.Label(window, text='Title:', font=('Arial Bold', 14), anchor = 'w', justify='left')
lbl.place(x=30, y=805)
title_textBox = tk.Entry(window, width=32)
title_textBox.place(x=110, y=810)

# Filter On Combo Box
lbl = tk.Label(window, text='Filter On:', font=('Arial Bold', 14), anchor = 'w', justify='left')
lbl.place(x=30, y=840)
filterBox = Combobox(window, width=25)
filterBox.place(x= 130, y = 845)
filterBox.bind("<<ComboboxSelected>>", filterChoiceCallback)

# Run Button Callback
plotList = ('time', 'lat', 'lon', 'easting', 'northing', \
            'crossTrack', 'alongTrack', 'z', 'classification', 'signalConf')
filterState = np.array([0,1,2,3,4])
def plotAtl03():
    
    # Try plot code
    try:
        
        # Get
        gtNumToPlot = gtNumPlotBox.current()
        
        # Get x,y combo box number selections
        xVarNum = xValsBox.current()
        yVarNum = yValsBox.current()
        
        # Get x,y combo bxx text selections
        xData = eval('atlMeasuredData[' + str(gtNumToPlot) + '].' + plotList[xVarNum])
        yData = eval('atlMeasuredData[' + str(gtNumToPlot) + '].' + plotList[yVarNum])
        
        # Get labels
        xLabel = xlabel_textBox.get()
        yLabel = ylabel_textBox.get()
        title = title_textBox.get()
        
        # Get Filter data type (classification or signal confidence) and filter numbers
        filterChoice = filterBox.get()
        if('  ' in filterChoice.lower()):
            filterType = []
            filterData = []
            filterNum = []
        elif('class' in filterChoice.lower()):
            filterType = filterChoice
            filterData = eval('atlMeasuredData[' + str(gtNumToPlot) + '].classification')
            filterTF = [filter0ChkState.get(), filter1ChkState.get(), filter2ChkState.get(), filter3ChkState.get(), filter4ChkState.get()]    
            filterNum = filterState[filterTF]
        elif('signal' in filterChoice.lower()):
            filterType = filterChoice
            filterData = eval('atlMeasuredData[' + str(gtNumToPlot) + '].signalConf')
            filterTF = [filter0ChkState.get(), filter1ChkState.get(), filter2ChkState.get(), filter3ChkState.get(), filter4ChkState.get()]    
            filterNum = filterState[filterTF]
        # endIf
        
        # Call getPlot function
        getPlot(xData, yData, xLabel, yLabel, title, filterType, filterData, filterNum)
    
    except:
        
        print('Cannot plot data. Please check inputs')
        
    # endTry
# endDef
    
# Plot Button
btn = tk.Button(window, text='PLOT', font=('Arial Bold', 20), width = 10, command=plotAtl03) 
btn.place(x=335, y=810)

# Open GUI window
window.mainloop()