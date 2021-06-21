# -*- coding: utf-8 -*-
"""
This script loads the PhoREAL GUI which:
    - Reads ICESat-2 ATL03 and ATL08 .h5 files
    - Reads Reference data in .las, .laz, and .tiff formats
    - Finds geolocation offsets in the ICESat-2 data with respect to the Reference data
    - Generates comparison and error plots
    - Creates an output.log file of all output from the console window 
    - Allows users to make their own output plots

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
from tkinter import ttk
from tkinter import messagebox
from tkinter.filedialog import asksaveasfile 

import numpy as np
import os
import sys
import time as runTime
from datetime import datetime
import webbrowser
import pandas as pd
from scipy import interpolate
import pickle as pkl

from getAtlMeasuredSwath_auto import getAtlMeasuredSwath
from getAtlTruthSwath_auto import getAtlTruthSwath
from getMeasurementError_auto import getMeasurementError, offsetsStruct

from icesatPlot import (getPlot, getPlot_atl08, getPlot_truth, 
                        getPlot_measCorr, getPklPlot, addStatsToPlot)
from icesatIO import (createHTMLChart, writeLog, write_mat, getTruthFilePaths,
                      getTruthHeaders, swbeamToGT, beamNumToGT)
from icesatUtils import (superFilter, getNameParts, interp_vals)
from icesatBin import get_bin_df

from gui_logo import images
from gui_addins import (viewerBlank_html, viewerBlankOnline_html)


# Print opening message
print('\n')
print('*************************************************')
print('PhoREAL GUI is running...')
print('*************************************************')
print('\n')


window = tk.Tk()

# GUI title
phoRealVersion = 'v3.27'
window.title('PhoREAL %s - Applied Research Labs (The University of Texas at Austin)' %phoRealVersion)

# GUI size
window.geometry('1140x500')
window.resizable(width=False, height=False)

# Set icon image for GUI
if os.name == 'nt':
    window.wm_iconbitmap(images)
# endIf

# Create top level tabs
tabControl = ttk.Notebook(window)  
       
# Create tab1
tab1 = ttk.Frame(tabControl)            
tabControl.add(tab1, text='Home')      
tabControl.pack(expand=1, fill='both') 

# Create tab2
tab2 = ttk.Frame(tabControl)            
tabControl.add(tab2, text='Stats')      
tabControl.pack(expand=1, fill='both') 

# Create tab3 
tab3 = ttk.Frame(tabControl)            
tabControl.add(tab3, text='Plot')      
tabControl.pack(expand=1, fill='both') 

# Create tab4 
tab4 = ttk.Frame(tabControl)            
tabControl.add(tab4, text='Help')      
tabControl.pack(expand=1, fill='both')

# Create tab5 
tab5 = ttk.Frame(tabControl)            
tabControl.add(tab5, text='About')      
tabControl.pack(expand=1, fill='both')

# Set tab font sytle
s = ttk.Style()
s.configure('TNotebook.Tab', font=('Arial','10','bold'))


###############################################################################
#
# TAB 1: HOME
#
###############################################################################

# Define global variables
global cwd, atl03FileExists, atl03FilePaths, atl08FileExists, atl08FilePaths, \
       outPathExists, truthDataExists, atl03DF_all, truthFilePaths, truthFileType, statsNumList
       
# Initialize parameters
cwd = os.path.normpath(os.getcwd())
atl03FileExists = False
atl03FilePaths = []
atl08FileExists = False
atl08FilePaths = []
outPathExists = True
truthDataExists = False
atl03DF_all = []
truthFilePaths = []
truthFileType = '.las'
statsNumList = np.array([], dtype=int)

### Panel label
measLabelframe = tk.LabelFrame(tab1, width=545, height=450, text='Get ICESat-2 Data Input', font=('Arial Bold', 14))
measLabelframe.place(x=15, y=10)

# Check ATL08 Button Callback
def checkATL03():
    
    global atl03FileExists, atl03FilePaths
    
    atl03FileInputRaw = atl03_textBox.get()
    atl03FileInputSplit = atl03FileInputRaw.split(',')
    atl03FileInput = [x.strip() for x in atl03FileInputSplit]
    atl03FilePaths = getTruthFilePaths(atl03FileInput, 'ATL03*.h5', logFileID=False)
    
    if(atl03FilePaths):
        atl03FileExists = True
    else:
        atl03FileExists = False
    # endIf
    
# endDef
    
# Function to determine the ATL03 file input type
def selectAtl03FileType(event):
    
    global atl03FilePaths
    
    atl03_textBox.delete(0,'end')
    atl03_textBox.insert(0,'')
    atl03FilePaths = []
    
# endDef
    
# Browse ATL03 File Button Callback
def browseAtl03():
    
    global atl03FileExists, atl03FilePaths
    
    # Get GUI input parameters
    atl03FileTypeValue = atl03FileTypeBox.get()
    
    # Determine whether to select a file or a directory
    if('file' in atl03FileTypeValue.lower()):
        # Option for .h5 file(s)
        helpStr = 'Select ATL03 .h5 File(s)'
        extStr = 'ATL03*.h5'
        atl03Input = filedialog.askopenfilenames(title = helpStr, filetypes = [(extStr, extStr)])
        inpIsDir = False
    else:
        # Option for .h5 directory
        helpStr = 'Select Directory with ATL03 .h5 File(s)'
        atl03Input = os.path.normpath(filedialog.askdirectory(title = helpStr))
        inpIsDir = True
    # endIf
    
    # Get ATL03 file paths from user input
    if(atl03Input!='.' and atl03Input!=''):
        atl03FilePaths = getTruthFilePaths(atl03Input, 'ATL03*.h5', logFileID=False)
        atl03_textBox.delete(0,'end')
        if(inpIsDir):
            atl03_textBox.insert(0,atl03Input)
        else:
            atl03_textBox.insert(0,(', ').join(atl03FilePaths))
        # endIf
        atl03FileExists = True
    else:
        atl03FileExists = False
    # endIf
# endDef
    
# ATL03 File Type Combo Box
atl03FileTypeBox = ttk.Combobox(measLabelframe, width=12)
atl03FileTypeBox.place(x=150, y=11)
atl03FileTypes = ['.h5 File(s)','.h5 Directory']
atl03FileTypeBox['values'] = atl03FileTypes
atl03FileTypeBox.bind("<<ComboboxSelected>>", selectAtl03FileType)
atl03FileTypeBox.current(0)

# ATL03 File Entry Box
lbl = tk.Label(measLabelframe, text='ATL03 File Type:', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=10)
atl03_textBox = tk.Entry(measLabelframe, width=70)
atl03_textBox.place(x=10, y=40)
#atl03_textBox.bind('<Return>', (lambda event: checkATL03()))
#atl03_textBox.bind('<Tab>', (lambda event: checkATL03()))

# ATL03 File Browse Button
atl03BrowseButton = tk.Button(measLabelframe, text='Browse', font=('Arial Bold', 12), command=browseAtl03) 
atl03BrowseButton.place(x=450, y=30)

# Check ATL08 Button Callback
def checkATL08():
    
    global atl08FileExists, atl08FilePaths
    
    atl08FileInputRaw = atl08_textBox.get()
    atl08FileInputSplit = atl08FileInputRaw.split(',')
    atl08FileInput = [x.strip() for x in atl08FileInputSplit]
    atl08FilePaths = getTruthFilePaths(atl08FileInput, 'ATL08*.h5', logFileID=False)
    
    if(atl08FilePaths):
        atl08FileExists = True
    else:
        atl08FileExists = False
    # endIf   
# endDef
    
# Function to determine the ATL03 file input type
def selectAtl08FileType(event):
    
    global atl08FilePaths
    
    atl08_textBox.delete(0,'end')
    atl08_textBox.insert(0,'')
    atl08FilePaths = []
    
# endDef
    
# Browse ATL08 File Button Callback
def browseAtl08():
    
    global atl08FileExists, atl08FilePaths
    
    # Get GUI input parameters
    currentData = atl08_textBox.get()
    atl08FileTypeValue = atl08FileTypeBox.get()
    useMeasErrorSection = useMeasErrorSectionChkState.get()
    
    # Determine if single ATL03 file or not
    singleAtl03File = False
    if(len(atl03FilePaths)==1):
        if(atl03FilePaths[0][-3:]=='.h5'):
            singleAtl03File = True
        # endIf
    # endIf
    
    # Get matching ATL08 file name path if single ATL03 file
    try:
        if(atl03FileExists and singleAtl03File):
            
            # Get ATL03 file and name parts
            atl03Path = atl03_textBox.get().strip()
            atl03File_w_ext = os.path.basename(atl03Path)
            atl03NameParts = getNameParts(atl03File_w_ext)
            
            # Build ATL08 .h5 string matching title of ATL03 file
            filterText = 'ATL08_' + atl03NameParts.year + atl03NameParts.month + \
            atl03NameParts.day + atl03NameParts.hour + atl03NameParts.minute + \
            atl03NameParts.second + '_' + atl03NameParts.trackNum + \
            atl03NameParts.unknown + '_' + atl03NameParts.releaseNum + '*.h5'
        else:
            
            # String for all ATL08 .h5 files
            filterText = 'ATL08_*.h5' 
        # endIf

    except:
        
        # String for all ATL08 .h5 files
        filterText = 'ATL08_*.h5'
    # endTry
        
    # Determine whether to select a file or directory
    if('file' in atl08FileTypeValue.lower()):
        # Option for .h5 file(s)
        helpStr = 'Select ATL08 .h5 File(s)'
        atl08Input = filedialog.askopenfilenames(title = helpStr, filetypes = [(filterText, filterText)])
        inpIsDir = False
    else:
        # Option for .h5 directory
        helpStr = 'Select Directory with ATL08 .h5 File(s)'
        atl08Input = os.path.normpath(filedialog.askdirectory(title = helpStr))
        inpIsDir = True
    # endIf
    
    # Get ATL08 file paths from user input
    if(atl08Input!='.' and atl08Input!=''):
        atl08FilePaths = getTruthFilePaths(atl08Input, 'ATL08*.h5', logFileID=False)
        atl08_textBox.delete(0,len(currentData))
        if(inpIsDir):
            atl08_textBox.insert(0,atl08Input)
        else:
            atl08_textBox.insert(0,(', ').join(atl08FilePaths))
        # endIf
        atl08FileExists = True
        if(useMeasErrorSection):
            useGroundIndex_checkBox.config(state = 'normal')
        # endIf
    else:
        atl08FileExists = False
        useGroundIndex_checkBox.config(state = 'disabled')
    # endIf
# endDef
    
# ATL08 File Type Combo Box
atl08FileTypeBox = ttk.Combobox(measLabelframe, width=12)
atl08FileTypeBox.place(x=150, y=71)
atl08FileTypes = ['.h5 File(s)','.h5 Directory']
atl08FileTypeBox['values'] = atl08FileTypes
atl08FileTypeBox.bind("<<ComboboxSelected>>", selectAtl08FileType)
atl08FileTypeBox.current(0)

# ATL08 File Entry Box
lbl = tk.Label(measLabelframe, text='ATL08 File Type:', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=70)
atl08_textBox = tk.Entry(measLabelframe, width=70) 
atl08_textBox.place(x=10, y=100)
lbl = tk.Label(measLabelframe, text='(Optional)', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=255, y=70)
#atl08_textBox.bind('<Return>', (lambda event: checkATL08()))
#atl08_textBox.bind('<Tab>', (lambda event: checkATL08()))
    
# ATL08 File Browse Button
atl08BrowseButton = tk.Button(measLabelframe, text='Browse', font=('Arial Bold', 12), command=browseAtl08) 
atl08BrowseButton.place(x=450, y=90)

# Browse Output File Button Callback
def browseOutput():
    
    currentData = outPath_textBox.get()
    outPath = os.path.normpath(filedialog.askdirectory(title = 'Select Output File Path'))
    if(outPath != '.'):
        outPath_textBox.delete(0,len(currentData))
        outPath_textBox.insert(0,outPath) 
    # endIf
    
# endDef

# Output Path Entry Box
lbl = tk.Label(measLabelframe, text='Output Directory:', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=130)
outPath_textBox = tk.Entry(measLabelframe, width=70) 
outPath_textBox.place(x=10, y=160)
#outPath_textBox.insert(0,cwd)

# Output Directory Browse Button
outPathBrowseButton = tk.Button(measLabelframe, text='Browse', font=('Arial Bold', 12), command=browseOutput) 
outPathBrowseButton.place(x=450, y=150)


### GT Number Text Label
lbl = tk.Label(measLabelframe, text='Ground Track Numbers:', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=190)

# Function to display GT types based on user choice
def gtTypeCallback(event):
    
    # Get GT selection type
    gtType = gtTypeBox.get()
    
    # Set GT selection text
    if(gtType == 'GT Right/Left'):
        gtNum1r_checkBox.config(text='GT1R')
        gtNum2r_checkBox.config(text='GT2R')
        gtNum3r_checkBox.config(text='GT3R')
        gtNum1l_checkBox.config(text='GT1L')
        gtNum2l_checkBox.config(text='GT2L')
        gtNum3l_checkBox.config(text='GT3L')
    elif(gtType == 'GT Strong/Weak'):
        gtNum1r_checkBox.config(text='GT1S')
        gtNum2r_checkBox.config(text='GT2S')
        gtNum3r_checkBox.config(text='GT3S')
        gtNum1l_checkBox.config(text='GT1W')
        gtNum2l_checkBox.config(text='GT2W')
        gtNum3l_checkBox.config(text='GT3W')
    elif(gtType == 'Beam Number'):
        gtNum1r_checkBox.config(text='# 1')
        gtNum2r_checkBox.config(text='# 3')
        gtNum3r_checkBox.config(text='# 5')
        gtNum1l_checkBox.config(text='# 2')
        gtNum2l_checkBox.config(text='# 4')
        gtNum3l_checkBox.config(text='# 6')
    # endIf
    
# endDef
    
# Function to get GT nums based on GT selection type and user input
def getGtNums(atl03FilePath):
    
    # Get GT selection type
    gtType = gtTypeBox.get()
    
    # Get GT user choices
    gtNumsTF = [gtNum1rChkState.get(), gtNum1lChkState.get(), gtNum2rChkState.get(), gtNum2lChkState.get(), gtNum3rChkState.get(), gtNum3lChkState.get()]    

    # Determine GT based on GT selection type and user input
    if(gtType == 'GT Right/Left'):
        
        gtNumsAll = np.array(['GT1R','GT1L','GT2R','GT2L','GT3R','GT3L'])
        gtNums = gtNumsAll[gtNumsTF]
        
    elif(gtType == 'GT Strong/Weak'):
        
        groupNumAll = np.array(['1','1','2','2','3','3'])
        groupNum = groupNumAll[gtNumsTF]
        swAll = np.array(['strong','weak','strong','weak','strong','weak'])
        sw = swAll[gtNumsTF]
        gtNums = swbeamToGT(atl03FilePath, groupNum, sw)
        
    elif(gtType == 'Beam Number'):
        
        beamNumAll = np.array(['1','2','3','4','5','6'])
        beamNum = beamNumAll[gtNumsTF]
        gtNums = beamNumToGT(atl03FilePath, beamNum)
        
    # endIf
    
    return gtNums    

# endDef

# GT Selection Type
lbl = tk.Label(measLabelframe, text='Select By:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=220)
gtTypeBox = ttk.Combobox(measLabelframe, width=16)
gtTypeBox.place(x= 90, y = 222)
gtTypes = ['GT Right/Left','GT Strong/Weak','Beam Number']
gtTypeBox['values'] = gtTypes
gtTypeBox.current(0)
gtTypeBox.bind("<<ComboboxSelected>>", gtTypeCallback)

### GT Number Checkboxes
gtNum1rChkState = tk.BooleanVar()
gtNum1rChkState.set(True)
gtNum1r_checkBox = tk.Checkbutton(measLabelframe, text = 'GT1R', font=('Arial', 12), var = gtNum1rChkState) 
gtNum1r_checkBox.place(x=250, y=190)

gtNum2rChkState = tk.BooleanVar()
gtNum2rChkState.set(True)
gtNum2r_checkBox = tk.Checkbutton(measLabelframe, text = 'GT2R', font=('Arial', 12), var = gtNum2rChkState) 
gtNum2r_checkBox.place(x=350, y=190)

gtNum3rChkState = tk.BooleanVar()
gtNum3rChkState.set(True)
gtNum3r_checkBox = tk.Checkbutton(measLabelframe, text = 'GT3R', font=('Arial', 12), var = gtNum3rChkState) 
gtNum3r_checkBox.place(x=450, y=190)

gtNum1lChkState = tk.BooleanVar()
gtNum1lChkState.set(False)
gtNum1l_checkBox = tk.Checkbutton(measLabelframe, text = 'GT1L', font=('Arial', 12), var = gtNum1lChkState) 
gtNum1l_checkBox.place(x=250, y=220)

gtNum2lChkState = tk.BooleanVar()
gtNum2lChkState.set(False)
gtNum2l_checkBox = tk.Checkbutton(measLabelframe, text = 'GT2L', font=('Arial', 12), var = gtNum2lChkState) 
gtNum2l_checkBox.place(x=350, y=220)

gtNum3lChkState = tk.BooleanVar()
gtNum3lChkState.set(False)
gtNum3l_checkBox = tk.Checkbutton(measLabelframe, text = 'GT3L', font=('Arial', 12), var = gtNum3lChkState) 
gtNum3l_checkBox.place(x=450, y=220)


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
lbl = tk.Label(measLabelframe, text='Trim ICESat-2 Data Options:', font=('Arial Bold', 12))
lbl.place(x=10, y=260)

### Trim Info Checkboxes
trimNoneModeChkState = tk.BooleanVar()
trimNoneModeChkState.set(True)
trimNoneMode_checkBox = tk.Checkbutton(measLabelframe, text = 'None', font=('Arial', 12), var = trimNoneModeChkState, command = checkNone) 
trimNoneMode_checkBox.place(x=10, y=290)

trimAutoModeChkState = tk.BooleanVar()
trimAutoModeChkState.set(False)
trimAutoMode_checkBox = tk.Checkbutton(measLabelframe, text = 'Auto', font=('Arial', 12), var = trimAutoModeChkState, command = checkAuto) 
trimAutoMode_checkBox.place(x=10, y=320)

trimManualModeChkState = tk.BooleanVar()
trimManualModeChkState.set(False)
trimManualMode_checkBox = tk.Checkbutton(measLabelframe, text = 'Manual', font=('Arial', 12), var = trimManualModeChkState, command = checkManual) 
trimManualMode_checkBox.place(x=100, y=290)

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
latMode_checkBox = tk.Checkbutton(measLabelframe, text = 'Latitude', font=('Arial', 12), var = latModeChkState, command = checkLat, state = 'disabled') 
latMode_checkBox.place(x=200, y=290)

### Lat Min Text Label
lbl = tk.Label(measLabelframe, text='Min', font=('Arial', 10))
lbl.place(x=310, y=270)

# Lat Min Text Entry
latMin_textBox = tk.Entry(measLabelframe, width=7, state = 'disabled') 
latMin_textBox.place(x=305, y=295)

### Lat Max Text Label
lbl = tk.Label(measLabelframe, text='Max', font=('Arial', 10))
lbl.place(x=385, y=270)

# Lat Max Text Entry
latMax_textBox = tk.Entry(measLabelframe, width=7, state = 'disabled') 
latMax_textBox.place(x=380, y=295)

### Lat Degrees Text Label
lbl = tk.Label(measLabelframe, text='Degrees', font=('Arial', 10))
lbl.place(x=440, y=295)

### Time Mode Checkbox
timeModeChkState = tk.BooleanVar()
timeModeChkState.set(False)
timeMode_checkBox = tk.Checkbutton(measLabelframe, text = 'Time', font=('Arial', 12), var = timeModeChkState, command = checkTime, state = 'disabled') 
timeMode_checkBox.place(x=200, y=320)

# Time Min Text Entry
timeMin_textBox = tk.Entry(measLabelframe, width=7, state = 'disabled') 
timeMin_textBox.place(x=305, y=325)

# Time Max Text Entry
timeMax_textBox = tk.Entry(measLabelframe, width=7, state = 'disabled') 
timeMax_textBox.place(x=380, y=325)

### Time Seconds Text Label
lbl = tk.Label(measLabelframe, text='Seconds', font=('Arial', 10))
lbl.place(x=440, y=325)

### Create Output Files Text Label
lbl = tk.Label(measLabelframe, text='Create Output Files:', font=('Arial Bold', 12))
lbl.place(x=10, y=360)

### Create ATL03 .las Text Entry
createLasChkState = tk.BooleanVar()
createLasChkState.set(True)
createLas_checkBox = tk.Checkbutton(measLabelframe, text = 'ATL03 .las File', font=('Arial', 12), var = createLasChkState) 
createLas_checkBox.place(x=10, y=390)

### Create ATL03 .kml Text Entry
createKmlChkState = tk.BooleanVar()
createKmlChkState.set(True)
createKml_checkBox = tk.Checkbutton(measLabelframe, text = 'ATL03 .kml File', font=('Arial', 12), var = createKmlChkState) 
createKml_checkBox.place(x=190, y=360)

### Create ATL03 .csv Text Entry
createCsvChkState = tk.BooleanVar()
createCsvChkState.set(False)
createCsv_checkBox = tk.Checkbutton(measLabelframe, text = 'ATL03 .csv File', font=('Arial', 12), var = createCsvChkState) 
createCsv_checkBox.place(x=190, y=390)

### Create ATL08 .kml Text Entry
createATL03PklChkState = tk.BooleanVar()
createATL03PklChkState.set(False)
createATL03Pkl_checkBox = tk.Checkbutton(measLabelframe, text = 'ATL03 .pkl File', font=('Arial', 12), var = createATL03PklChkState) 
createATL03Pkl_checkBox.place(x=370, y=360)

### Create ATL08 .csv Text Entry
createATL08CsvChkState = tk.BooleanVar()
createATL08CsvChkState.set(False)
createATL08Csv_checkBox = tk.Checkbutton(measLabelframe, text = 'ATL08 .csv File', font=('Arial', 12), var = createATL08CsvChkState) 
createATL08Csv_checkBox.place(x=370, y=390)
    

###############################################################################
#
# TRUTH SECTION
#
###############################################################################

### TRUTH DATA Panel label
truthLabelframe = tk.LabelFrame(tab1, width=545, height=120, text='Get Reference Data Input', font=('Arial Bold', 14))
truthLabelframe.place(x=580, y=10)

# Use reference section callback
def useTruthSectionCallback():
    useTruthSection = useTruthSectionChkState.get()
    if(useTruthSection):
        for child in truthLabelframe.winfo_children():
            child.configure(state='normal')
        # endFor
        checkUseExistingTruth()
    else:
        for child in truthLabelframe.winfo_children():
            child.configure(state='disable')
        # endFor    
        useMeasErrorSectionChkState.set(False)
        useMeasErrorSectionCallback()
    # endIf
# endDef

# Browse Truth Data Directory File Button Callback
def browseTruthDir():
    
    global truthDataExists, truthFilePaths, truthFileType
    
    # Get GUI input parameters
    currentData = truthDataDir_textBox.get()
    truthFileTypeValue = truthFileTypeBox.get()
    inpIsDir = False
    
    # Set listbox options depending on user input
    if('file' in truthFileTypeValue.lower()):
        # Option for reference file(s)
        if('(s)' in truthFileTypeValue.lower()):
            # Option for las/tif reference files
            if('las' in truthFileTypeValue.lower()):
                helpStr = 'Select Reference .las File(s)'
                extStr = '*.las'
                truthFileType = '.las'
            else:
                helpStr = 'Select Reference .tif File(s)'
                extStr = '*.tif'
                truthFileType = '.tif'
            # endIf
            truthData = filedialog.askopenfilenames(title = helpStr, filetypes = [(extStr, extStr)])
            truthDataExists = True
        else:
            # Option for las/tif buffer reference file 
            if('las' in truthFileTypeValue.lower()):
                helpStr = 'Select Reference Buffer .las File'
                extStr = '*buffer.las'
                truthFileType = '*buffer.las'
            else:
                helpStr = 'Select Reference Buffer .tif File'
                extStr = '*.tif'
                truthFileType = '.tif'
            # endIf
            truthData = os.path.normpath(filedialog.askopenfilename(title = helpStr, filetypes = [(extStr, extStr)]))
            truthDataExists = os.path.exists(truthData)
        # endIf
    else:
        # Reference directory option
        if('las' in truthFileTypeValue.lower()):
            helpStr = 'Select Reference Directory with .las File(s)'
            truthFileType = '.las'
        else:
            helpStr = 'Select Reference Directory with .tif File(s)'
            truthFileType = '.tif'
        # endIf
        truthData = os.path.normpath(filedialog.askdirectory(title = helpStr))
        truthDataExists = os.path.isdir(truthData)
        inpIsDir = True
    # endIf
    
    if(truthDataExists and truthData!='.' and truthData!=''):
        truthFilePaths = getTruthFilePaths(truthData, truthFileType, logFileID=False)
        truthDataDir_textBox.delete(0,len(currentData))
        if(inpIsDir):
            truthDataDir_textBox.insert(0,truthData)
        else:
            truthDataDir_textBox.insert(0,(', ').join(truthFilePaths))
        # endIf
    else:
        truthDataExists = False
    # endIf
# endDef
    
# Use Existing Truth button callback
def checkUseExistingTruth():
    
    global truthFilePaths
    
    # Set values in listbox
    useTruth = useExistingTruthChkState.get()
    if(useTruth):
        truthBuffer_textBox.config(state = 'disabled')
        truthFileTypeBox['values'] = ['.las File']
    else:
        truthBuffer_textBox.config(state = 'normal')
        truthFileTypeBox['values'] = ['.las File(s)','.las Directory','.tif File(s)','.tif Directory']
    # endIf
    truthFileTypeBox.current(0)
    
    # Clear reference entry box
    truthDataDir_textBox.delete(0,'end')
    truthDataDir_textBox.insert(0,'')
    truthFilePaths = []
# endDef
    
# Browse Truth Data Directory File Button Callback
def checkTruthDir():    
    
    global truthDataExists, truthFilePaths
    
    truthFileInputRaw = truthDataDir_textBox.get()
    truthFileInputSplit = truthFileInputRaw.split(',')
    truthFileInput = [x.strip() for x in truthFileInputSplit]
    truthFilePaths = getTruthFilePaths(truthFileInput, truthFileType, logFileID=False)
    
    if(truthFilePaths):
        truthDataExists = True
    else:
        truthDataExists = False
    # endIf
    
# endDef
    
# Function to clear listbox contents when selected
def selectTruthFileType(event):
    
    global truthFilePaths, truthFileType
    
    truthDataDir_textBox.delete(0,'end')
    truthDataDir_textBox.insert(0,'')
    truthFilePaths = []
    
    truthFileTypeValue = truthFileTypeBox.get()
    if('las' in truthFileTypeValue.lower()):
        truthFileType = '.las'
    else:
        truthFileType = '.tif'
    # endIf
    
# endDef

### Activate Truth section
useTruthSectionChkState = tk.BooleanVar()
useTruthSectionChkState.set(False)
useTruthSection_checkBox = tk.Checkbutton(tab1, text = '', font=('Arial', 12), var = useTruthSectionChkState, command = useTruthSectionCallback) 
useTruthSection_checkBox.place(x=1110, y=10)

### Use Existing Truth Check Box
useExistingTruthChkState = tk.BooleanVar()
useExistingTruthChkState.set(False)
useExistingTruth_checkBox = tk.Checkbutton(truthLabelframe, text = 'Use Existing Buffer', font=('Arial', 12), var = useExistingTruthChkState, command = checkUseExistingTruth) 
useExistingTruth_checkBox.place(x=10, y=10)

### Truth Buffer Size Entry Box
lbl = tk.Label(truthLabelframe, text='Buffer Size (m):', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=180, y=12)
truthBuffer_textBox = tk.Entry(truthLabelframe, width=5)
truthBuffer_textBox.place(x=300, y=16)
truthBuffer_textBox.insert(0,'50')
truthBuffer_textBox.config(state = 'disabled')

### Create Output Truth File Check Box
createTruthFileChkState = tk.BooleanVar()
createTruthFileChkState.set(True)
createTruthFile_checkBox = tk.Checkbutton(truthLabelframe, text = 'Save Reference File', font=('Arial', 12), var = createTruthFileChkState) 
createTruthFile_checkBox.place(x=350, y=10)

### Truth file Type listbox
lbl = tk.Label(truthLabelframe, text='Type:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=12, y=53)
truthFileTypeBox = ttk.Combobox(truthLabelframe, width=12)
truthFileTypeBox.place(x= 60, y = 55)
truthFileTypes = ['.las File(s)','.las Directory','.tif File(s)','.tif Directory']
truthFileTypeBox['values'] = truthFileTypes
truthFileTypeBox.bind("<<ComboboxSelected>>", selectTruthFileType)
    
### Truth Data Directory Entry Box
truthDataDir_textBox = tk.Entry(truthLabelframe, width=43)
truthDataDir_textBox.place(x=170, y=55)
#truthDataDir_textBox.bind('<Return>', (lambda event: checkTruthDir()))
#truthDataDir_textBox.bind('<Tab>', (lambda event: checkTruthDir()))

### Truth Data Dir Browse Button
truthDataDirBrowseButton = tk.Button(truthLabelframe, text='Browse', font=('Arial Bold', 12), command=browseTruthDir) 
truthDataDirBrowseButton.place(x=450, y=50)

# Set Truth section to disabled (default mode)
for child in truthLabelframe.winfo_children():
    child.configure(state='disable')
# endFor
    

###############################################################################
#
# MEASUREMENT ERROR SECTION
#
###############################################################################

### MEASUREMENT ERROR Panel label
measErrorLabelframe = tk.LabelFrame(tab1, width=545, height=250, text='Find ICESat-2 Offsets Relative to Reference Data', font=('Arial Bold', 14))
measErrorLabelframe.place(x=580, y=140)

# Use Measurement Error section callback
def useMeasErrorSectionCallback():
    useMeasErrorSection = useMeasErrorSectionChkState.get()
    if(useMeasErrorSection):
        for child in measErrorLabelframe.winfo_children():
            child.configure(state='normal')
        # endFor
        useVertShiftCallback()
        useMeasSigConfCallback()
        useGroundIndexConfCallback()
        if(not atl08FileExists):
            useGroundIndex_checkBox.config(state = 'disabled')
        # endIf
        
        # Enable truth section
        if(not(useTruthSectionChkState.get())):
            useTruthSectionChkState.set(True)
            useTruthSectionCallback()
        # endIf
    else:
        for child in measErrorLabelframe.winfo_children():
            child.configure(state='disable')
        # endFor
        useMeasErrorSectionChkState.set(False)
    # endIf
# endDef
    
### Activate Measurement Error section
useMeasErrorSectionChkState = tk.BooleanVar()
useMeasErrorSectionChkState.set(False)
useMeasErrorSection_checkBox = tk.Checkbutton(tab1, text = '', font=('Arial', 12), var = useMeasErrorSectionChkState, command = useMeasErrorSectionCallback) 
useMeasErrorSection_checkBox.place(x=1110, y=140)

### Cross-Track Bounds Entry Box
lbl = tk.Label(measErrorLabelframe, text='Cross-Track Bounds (m):', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=10)
crossTrackBounds_textBox = tk.Entry(measErrorLabelframe, width=10)
crossTrackBounds_textBox.place(x=195, y=15)
crossTrackBounds_textBox.insert(0,'-50, 50')

### Along-Track Bounds Entry Box
lbl = tk.Label(measErrorLabelframe, text='Along-Track Bounds (m):', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=270, y=10)
alongTrackBounds_textBox = tk.Entry(measErrorLabelframe, width=10)
alongTrackBounds_textBox.place(x=455, y=15)
alongTrackBounds_textBox.insert(0,'-50, 50')

### Multi-Resolutional Stepdown Entry Box
lbl = tk.Label(measErrorLabelframe, text='Grid Resolution(s) (m):', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=80)
multiresStepdown_textBox = tk.Entry(measErrorLabelframe, width=9)
multiresStepdown_textBox.place(x=175, y=85)
multiresStepdown_textBox.insert(0,'8, 4, 2, 1')

### Reference Height Selection
lbl = tk.Label(measErrorLabelframe, text='Use ICESat-2 Heights:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=245, y=80)
refHeightTypeBox = ttk.Combobox(measErrorLabelframe, width=17)
refHeightTypeBox.place(x= 410, y = 82)
refHeightTypes = ['Ellipsoidal (HAE)','Orthometric (MSL)']
refHeightTypeBox['values'] = refHeightTypes
refHeightTypeBox.current(0)

# Use vertical shift checkbox callback
def useVertShiftCallback():
    useVertShift = useFixedVertShiftChkState.get()
    if(useVertShift):
        verticalShift_textBox.config(state = 'normal')
    else:
        verticalShift_textBox.config(state = 'disabled')
    # endIf
# endDef
    
### Use Fixed Vertical Shift Check Box
useFixedVertShiftChkState = tk.BooleanVar()
useFixedVertShiftChkState.set(False)
useFixedVertShift_checkBox = tk.Checkbutton(measErrorLabelframe, text = 'Use Fixed Vertical Shift', font=('Arial', 12), var=useFixedVertShiftChkState, command=useVertShiftCallback) 
useFixedVertShift_checkBox.place(x=10, y=43)

### Vertical Shift Amount Entry Box
lbl = tk.Label(measErrorLabelframe, text='Vertical Shift (m):', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=240, y=45)
verticalShift_textBox = tk.Entry(measErrorLabelframe, width=10)
verticalShift_textBox.place(x=370, y=50)
verticalShift_textBox.insert(0,'0')
verticalShift_textBox.config(state = 'disabled')

# Use measured signal confidence checkbox callback
def useMeasSigConfCallback():
    useMeasSigConf = useMeasSigConfChkState.get()
    if(useMeasSigConf):
        measSigConfIndex_textBox.config(state = 'normal')
        truthGroundIndex_textBox.config(state = 'disabled')
        useGroundIndexChkState.set(False)
    else:
        measSigConfIndex_textBox.config(state = 'disabled')
        truthGroundIndex_textBox.config(state = 'normal')
        useGroundIndexChkState.set(True)
    # endIf
# endDef
    
### Use Measured Signal Confidence Entry Box
useMeasSigConfChkState = tk.BooleanVar()
useMeasSigConfChkState.set(True)
useMeasSigConf_checkBox = tk.Checkbutton(measErrorLabelframe, text = 'Use ICESat-2 Signal Confidence Value(s):', font=('Arial', 12), var=useMeasSigConfChkState, command=useMeasSigConfCallback) 
useMeasSigConf_checkBox.place(x=10, y=115)

### Measured Signal Confidence Entry Box
measSigConfIndex_textBox = tk.Entry(measErrorLabelframe, width=8)
measSigConfIndex_textBox.place(x=340, y=120)
measSigConfIndex_textBox.insert(0,'3, 4')

# Use measured signal confidence checkbox callback
def useGroundIndexConfCallback():
    useGroundIndex = useGroundIndexChkState.get()
    if(useGroundIndex):
        truthGroundIndex_textBox.config(state = 'normal')
        measSigConfIndex_textBox.config(state = 'disabled')
        useMeasSigConfChkState.set(False)
    else:
        truthGroundIndex_textBox.config(state = 'disabled')
        measSigConfIndex_textBox.config(state = 'normal')
        useMeasSigConfChkState.set(True)
    # endIf
# endDef
    
### Use Truth Ground Index Entry Box
useGroundIndexChkState = tk.BooleanVar()
useGroundIndexChkState.set(False)
useGroundIndex_checkBox = tk.Checkbutton(measErrorLabelframe, text = 'Use Reference Ground Index:', font=('Arial', 12), var=useGroundIndexChkState, command=useGroundIndexConfCallback) 
useGroundIndex_checkBox.place(x=10, y=150)

### Truth Ground Index Entry Box
truthGroundIndex_textBox = tk.Entry(measErrorLabelframe, width=5)
truthGroundIndex_textBox.place(x=245, y=155)
truthGroundIndex_textBox.insert(0,'2')
truthGroundIndex_textBox.config(state = 'disabled')
lbl = tk.Label(measErrorLabelframe, text='(Requires ATL08 File)', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=285, y=150)

### Create Measured Corrected File Check Box
createMeasCorrFileChkState = tk.BooleanVar()
createMeasCorrFileChkState.set(True)
createMeasCorrFile_checkBox = tk.Checkbutton(measErrorLabelframe, text = 'Save Shifted ICESat-2 File', font=('Arial', 12), var=createMeasCorrFileChkState) 
createMeasCorrFile_checkBox.place(x=10, y=185)

### Make Plots Check Box
makePlotsChkState = tk.BooleanVar()
makePlotsChkState.set(True)
makePlots_checkBox = tk.Checkbutton(measErrorLabelframe, text = 'Make Output Plots', font=('Arial', 12), var=makePlotsChkState) 
makePlots_checkBox.place(x=280, y=185)

# Set measurement correction section to disabled (default mode)
for child in measErrorLabelframe.winfo_children():
    child.configure(state='disabled')
# endFor
    
# Function to refresh stats
def refreshStats():
    
    global statsDF
    
    # Update status bar
    statsStatusBar['value'] = 0
    statuslbl.config(text='')
    
    # Update status bar
    statsStatusBarRH['value'] = 0
    statuslblRH.config(text='   Status:')

    # Clear Add Stats listbox in Stats Section
    addStatsTuple = ('')
    addStatsBox['values'] = addStatsTuple
    addStatsBox.delete(0,'end')
    
    # Clear stats
    statsDF = []
    
    # Update window
    window.update()
                
# end Def
    
# Function to load ATL03 plot info after a run
def loadAtl03_info():
    
    # Set file name being analyzed
    fileName_textBox.delete(0,'end')
    fileName_textBox.insert(0,atl03Data[0].atl03FileName)
    
    # Set file name being analyzed (in Stats Tab section)
    fileName_textBoxRH.delete(0,'end')
    fileName_textBoxRH.insert(0,atl03Data[0].atl03FileName)
    
    # Set Ground Tracks to plot
    gtNumsText = [i + ' (Beam #' + j + ', ' + k + ' beam)' for i, j, k in zip(gtNumsGood, beamNumsGood, beamStrengthGood)]
    gtNumsTuple = tuple(gtNumsText)
    gtNumPlotBox['values'] = gtNumsTuple
    gtNumPlotBox.current(0)
    
    # Set Ground Tracks to plot (in Stats Tab section)
    gtNumPlotBoxRH['values'] = gtNumsTuple
    gtNumPlotBoxRH.current(0)
    
    if(atl03Data[0].zone=='3413' or atl03Data[0].zone=='3976'):
        
        plotVarsTuple = ('Time (sec)', 'Delta Time (sec)', \
                         'Latitude (deg)', 'Longitude (deg)', \
                         'Polar Stereo X (m)', 'Polar Stereo Y (m)', \
                         'Cross-Track (m)', 'Along-Track (m)', \
                         'Height (m HAE)', 'Height (m MSL)', \
                         'Classification', 'Signal Confidence')
        
    else:
        
        plotVarsTuple = ('Time (sec)', 'Delta Time (sec)', \
                         'Latitude (deg)', 'Longitude (deg)', \
                         'UTM Easting (m)', 'UTM Northing (m)', \
                         'Cross-Track (m)', 'Along-Track (m)', \
                         'Height (m HAE)', 'Height (m MSL)', \
                         'Classification', 'Signal Confidence')
        
    # endIf
        
    # Set X Vals to plot
    xValsBox['values'] = plotVarsTuple
    xValsBox.current(0)
    
    # Set Y Vals to plot
    yValsBox['values'] = plotVarsTuple
    yValsBox.current(8)
    
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
    
    # Set Vals to filter on
    if(atl03Data[0].dataIsMapped):
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
    
    # Refresh stats
    refreshStats()
    
    # Update window
    window.update()
    
# endDef
           
# Function to load ATL08 plot info after a run
def loadAtl08_info():

    # Set Y Vals to plot
    yValsTuple_atl08 = ('Max Canopy (m HAE)', 'Terrain Best Fit (m HAE)', 'Terrain Median (m HAE)',
                        'Max Canopy (m MSL)', 'Terrain Best Fit (m MSL)', 'Terrain Median (m MSL)')
    yValsBox_atl08['values'] = yValsTuple_atl08
    yValsBox_atl08.current(0)
    
    # Set Segment By in Stats sections
    segmentByTuple = ('Time (sec)','Latitude (deg)','UTM Northing (m)','Along-Track (m)')
    segmentByBox['values'] = segmentByTuple
    segmentByBox.current(0)
    
    # Set increment units in Stats section
    currentData = incrementBox.get()
    incrementBox.delete(0,len(currentData))
    incrementBox.insert(0,'1')
    incrementText.config(text = 'sec')

    # Refresh stats
    refreshStats()
    
    # Update window
    window.update()
                
# endDef
    
# Function to check and match up ATL03/ATL08 file paths
def getMatchingAtlFiles(atl03FilePaths, atl08FilePaths):
    
    # Get ATL08 formatted string to compare against
    atl08FilePathsFormatted = np.empty(np.shape(atl08FilePaths), dtype=object)
    for i in range(0,len(atl08FilePathsFormatted)):
        
        # Get ATL03 file
        atl08FilePath = atl08FilePaths[i]
        
        # Get ATL03 file parts
        atl08File_w_ext = os.path.basename(atl08FilePath)
        atl08NameParts = getNameParts(atl08File_w_ext)
        
        # Build ATL08 formatted string
        atl08FilePathsFormatted[i] = 'ATL03_' + atl08NameParts.year + atl08NameParts.month + \
        atl08NameParts.day + atl08NameParts.hour + atl08NameParts.minute + \
        atl08NameParts.second + '_' + atl08NameParts.trackNum + \
        atl08NameParts.unknown
                
    # endFor
    
    # Loop through each ATL03 file and find matching ATL08 file
    atl03FilePathsAll = np.empty(np.shape(atl03FilePaths), dtype=object)
    atl08FilePathsAll = np.empty(np.shape(atl03FilePaths), dtype=object)
    for i in range(0,len(atl03FilePaths)):
        
        # Get ATL03 file
        atl03FilePath = atl03FilePaths[i]
        
        # Get ATL03 file parts
        atl03File_w_ext = os.path.basename(atl03FilePath)
        atl03NameParts = getNameParts(atl03File_w_ext)
        
        # Build ATL03 formatted string
        atl03FilePathsFormatted = 'ATL03_' + atl03NameParts.year + atl03NameParts.month + \
        atl03NameParts.day + atl03NameParts.hour + atl03NameParts.minute + \
        atl03NameParts.second + '_' + atl03NameParts.trackNum + \
        atl03NameParts.unknown
        
        # Compare ATL03 formatted string to ATL08 formatted string
        atl03FilePathsAll[i] = atl03FilePaths[i]
        if(np.isin(atl03FilePathsFormatted,atl08FilePathsFormatted)):
            ind = np.where(atl08FilePathsFormatted==atl03FilePathsFormatted)[0][0]
            atl08FilePathsAll[i] = atl08FilePaths[ind]
        else:
            atl08FilePathsAll[i] = ''
        # endIf
    # endFor
    
    return atl03FilePathsAll, atl08FilePathsAll
    
# endDef

# Run Button Callback
def runAtl03():
    
    # Update status bar
    statusBar['value'] = 0
    window.update()
        
    # Make atlMeasuredData a global variable
    global atl03Data, atl08Data, atlTruthDataFiltered, \
    truthDataExists, atlCorrections, atl03FileExists, atl08FileExists, \
    outFilePath, outPathExists, atl03DF_all, \
    gtNumsGood, beamNumsGood, beamStrengthGood
    
    # Check ATL03/ATL08 files
    checkATL03()
    checkATL08()
    atl03FilePathsAll, atl08FilePathsAll = getMatchingAtlFiles(atl03FilePaths,atl08FilePaths)
    
    # Check output file path, make output directory if it doesn't exist
    outFilePath = outPath_textBox.get().strip()
    if('' == outFilePath or '.' == outFilePath):
        outPathExists = False
    else:
        outPathExists = os.path.isdir(outFilePath)
        if(not outPathExists):
            os.mkdir(outFilePath)
            outPathExists = True
        # endIf
    # endIf
    
    # Get reference section status
    useTruthSection = useTruthSectionChkState.get()
    
    # Check Truth Directory
    truthOkToRun = True
    if(useTruthSection):
        checkTruthDir()
        if(not truthDataExists):
            truthOkToRun = False
        # endIf
    # endIf
          
    if(atl03FileExists and outPathExists and truthOkToRun):
    
        # Disable run button
        RunButton.config(state=tk.DISABLED)
        
        # Open .log file for writing
        logFileName = 'temp.txt'
        logFilePath = os.path.normpath(outFilePath + '/' + logFileName)
        if(os.path.exists(logFilePath)):
            os.remove(logFilePath)
        # endIf
        try:
            logFileID = open(logFilePath, 'w')
        except:
            logFileID = False
        # endTry
        
        # Try code
        try:
            
            # Start timer
            timeStart = runTime.time()
        
            # Get current date/time
            currentDateTime = (datetime.now()).strftime('%m/%d/%Y %H:%M:%S')
            writeLog('-------------------------------------', logFileID)
            writeLog('Run Initiated On: %s' %currentDateTime, logFileID)
            writeLog('-------------------------------------', logFileID)
            writeLog('', logFileID)
        
            # Loop through all ATL03/ATL08 files
            totalFiles = len(atl03FilePathsAll)
            for numFile in range(0,totalFiles):
                
                # Initialize variables
                atl03Data = []
                atl03DataSingle = []
                atl03DF = []
                atl03DF_all = []
                atl08Data = []
                atl08DataSingle = []
                atlTruthDataSingle = []
                atlTruthDataFiltered = []
                atlTruthDataFilteredSingle = []
                atlCorrections = []
                atlCorrectionsSingle = []
                gtNumsGood = []
                beamNumsGood = []
                beamStrengthGood = []
                
                # Get ATL03/ATL08 file paths
                atl03FilePath = atl03FilePathsAll[numFile] 
                atl08FilePath = atl08FilePathsAll[numFile]
                
                # Print updates
                writeLog('File #%d of %d' %(numFile+1,totalFiles),logFileID)
                writeLog('File Name: %s' %(os.path.basename(atl03FilePath)),logFileID)
                writeLog('', logFileID)         
                         
                # Get GT Nums
                gtNums = getGtNums(atl03FilePath)
            
                # Get Truth Data inputs
                if(useTruthSection):
                    useExistingTruth = useExistingTruthChkState.get()
                    # truthSwathDir = truthDataDir_textBox.get().strip()
                    bufferInput = truthBuffer_textBox.get().strip()
                    if('' == bufferInput):
                        buffer = 0
                    else:
                        buffer = int(bufferInput)
                    # endIf
                    createTruthFile = createTruthFileChkState.get()
                # endIf
                
                # Get Corrected Measured inputs
                useMeasErrorSection = useMeasErrorSectionChkState.get()
                if(useMeasErrorSection):
                    offsetsCrossTrackBounds = eval('[' + crossTrackBounds_textBox.get().strip() + ']')
                    offsetsAlongTrackBounds = eval('[' + alongTrackBounds_textBox.get().strip() + ']')
                    offsetsRasterResolutions = eval('[' + multiresStepdown_textBox.get().strip() + ']')
                    refHeightType = refHeightTypeBox.get()
                    offsetsUseVerticalShift = useFixedVertShiftChkState.get()
                    offsetsVerticalShift = float(verticalShift_textBox.get().strip())
                    offsets = offsetsStruct(offsetsCrossTrackBounds, offsetsAlongTrackBounds, offsetsRasterResolutions, offsetsUseVerticalShift, offsetsVerticalShift)
                    useMeasSigConf = useMeasSigConfChkState.get()
                    if(useMeasSigConf):
                        filterData = eval('[' + measSigConfIndex_textBox.get().strip() + ']')
                    else:
                        filterData = int(truthGroundIndex_textBox.get().strip())
                    # endIf
                    createMeasCorrFile = createMeasCorrFileChkState.get()
                    makePlots = makePlotsChkState.get()
                    showPlots = False
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
                else:
                    trimInfo = 'none'
                    trimNoneModeChkState.set(True)
                # endIf
                    
                # Get output file options
                createLasFile = createLasChkState.get()
                createKmlFile = createKmlChkState.get()
                createCsvFile = createCsvChkState.get()
                createATL08KmlFile = False
                createATL03PklFile = createATL03PklChkState.get()
                createATL08CsvFile = createATL08CsvChkState.get()
                    
                # Read truth data once before looping over other GT nums
                if(useTruthSection and truthDataExists and not(useExistingTruth)):
                    
                    # Get truth file header info
                    writeLog('Reading Reference Header Data...\n', logFileID)
                    truthHeaderDF = getTruthHeaders(truthFilePaths, truthFileType, logFileID)
                
                elif(useTruthSection and truthDataExists and useExistingTruth):
                    
                    truthHeaderDF = False
                            
                # endIf
                    
                # Loop through all GT nums
                rotationData = False
                for i in range(0,len(gtNums)):
                    
                    # Update status bar
                    statusBar['value'] = 0
                    window.update()
                
                    # Get GT Nums
                    gtNum = gtNums[i].lower()
                        
                    # Get ICESat-2 data
                    writeLog('Reading ICESat-2 Data...\n', logFileID)
                    atl03DataSingle, atl08DataSingle, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, 
                                                                                         gtNum, trimInfo, createLasFile, createKmlFile, 
                                                                                         createATL08KmlFile, createCsvFile, createATL08CsvFile, 
                                                                                         logFileID)
                    
                    # Update status bar
                    progressNum = 30
                    statusBar['value'] = progressNum
                    window.update()
                    
                    if(bool(atl03DataSingle)):
                        
                        # Create ATL03 array for dataframe
                        atl03Array = np.column_stack([atl03DataSingle.time,
                                                      atl03DataSingle.deltaTime,
                                                      atl03DataSingle.lat,
                                                      atl03DataSingle.lon,
                                                      atl03DataSingle.easting,
                                                      atl03DataSingle.northing,
                                                      atl03DataSingle.crossTrack,
                                                      atl03DataSingle.alongTrack,
                                                      atl03DataSingle.z,
                                                      atl03DataSingle.zMsl,
                                                      atl03DataSingle.classification])
            
                        # Set ATL03 column names for dataframe
                        colNames = ['Time (sec)','Delta Time (sec)',
                                    'Latitude (deg)','Longitude (deg)',
                                    'UTM Easting (m)','UTM Northing (m)',
                                    'Cross-Track (m)','Along-Track (m)',
                                    'Height (m HAE)','Height (m MSL)',
                                    'classification']
                        
                        # Create ATL03 dataframe
                        atl03DF = pd.DataFrame(data=atl03Array, columns=colNames)
                        
                        # Append measured data for each run
                        atl03Data.append(atl03DataSingle)
                        atl03DF_all.append(atl03DF)
                        gtNumsGood.append(gtNum)
                        beamNumsGood.append(atl03DataSingle.beamNum)
                        beamStrengthGood.append(atl03DataSingle.beamStrength)
                        if(atl08FileExists):
                            atl08Data.append(atl08DataSingle)
                        # endIf
                        
                        # Write output stats file if selected
                        if(createStatsFileChkState.get() and len(statsNumList)>0 and atl03Data[i].dataIsMapped):
                            try:
                                writeLog('   Writing output stats .csv file...\n', logFileID)
                                outputCsvBaseName = atl03DataSingle.atl03FileName + '_' + atl03DataSingle.gtNum + '_stats.csv'
                                outputCsvName = os.path.normpath(outFilePath + '/' + outputCsvBaseName)
                                computeStatsRH(outputCsvName, verbose=False)
                            except:
                                pass
                            # endTry
                        # endIf
                    
                    # endIf
    
                    if(bool(atl03DataSingle) and useTruthSection and truthDataExists):
                            
                        writeLog('Reading Reference Data...\n', logFileID)
                            
                        # Call getAtlTruthSwath
                        atlTruthDataSingle = getAtlTruthSwath(atl03DataSingle, rotationData, truthHeaderDF, 
                                                              truthFilePaths, buffer, outFilePath, createTruthFile, 
                                                              truthFileType, useExistingTruth, logFileID)
    
                        # Run superfilter on data
                        writeLog('   Filtering Reference Data...', logFileID)
                        atlTruthDataFilteredSingle, _ = superFilter(atl03DataSingle, atlTruthDataSingle, xBuf = 1, classCode = [])
                        writeLog('\n', logFileID)
                
                        # Interpolate to get reference time and deltaTime
                        atlTruthDataFilteredSingle.time = interp_vals(np.ravel(atl03DataSingle.alongTrack),
                                                                      np.ravel(atl03DataSingle.time), 
                                                                      np.ravel(atlTruthDataFilteredSingle.alongTrack))
                        
                        atlTruthDataFilteredSingle.deltaTime = interp_vals(np.ravel(atl03DataSingle.alongTrack),
                                                                           np.ravel(atl03DataSingle.deltaTime),
                                                                           np.ravel(atlTruthDataFilteredSingle.alongTrack))
                        
                        # Append filtered reference data for each run
                        atlTruthDataFiltered.append(atlTruthDataFilteredSingle)
                        
                        # Update status bar
                        progressNum = 60
                        statusBar['value'] = progressNum
                        window.update()
                    
                    # endIf
                    
                    if(bool(atl03DataSingle) and bool(atlTruthDataSingle) and useMeasErrorSection):
                        
                        # Get ICESat-2 geolocation offsets
                        writeLog('Finding ICESat-2 Offsets...\n', logFileID)
                        atlTruthEmpty = (len(atlTruthDataFilteredSingle.easting)==0) and (len(atlTruthDataFilteredSingle.northing)==0)
                        if(not atlTruthEmpty):
                            atlCorrectionsSingle = getMeasurementError(atl03DataSingle, atlTruthDataSingle, refHeightType, 
                                                                       rotationData, outFilePath, useMeasSigConf, filterData, 
                                                                       offsets, createMeasCorrFile, makePlots, showPlots, logFileID)
                        else:
                            writeLog('   WARNING: No offsets to find since no truth data.\n', logFileID)
                        # endIf
                        
                        # Append corrected measured data for each run
                        atlCorrections.append(atlCorrectionsSingle)
    
                        # Update status bar
                        progressNum = 90
                        statusBar['value'] = progressNum
                        window.update()
                    
                    # endIf
                    
                    # Update status bar
                    progressNum = 100
                    statusBar['value'] = progressNum
                    window.update()
                               
                    if(atl03DataSingle):
                        
                        # Store data into .mat file
                        matFileName = atl03DataSingle.atl03FileName + '_' + atl03DataSingle.gtNum + '.mat'
                        matFilePath = os.path.normpath(outFilePath + '/' + matFileName)
                        mat_atl03FileName = atl03DataSingle.atl03FileName
                        mat_gtNum = atl03DataSingle.gtNum
                        mat_trackDirection = atl03DataSingle.trackDirection
                    
                    # endIf
                    
                    # Prepare output .mat file data
                    if(atlCorrectionsSingle):
                        
                        # Collect .mat file data
                        mat_eastingCorrection = np.array([np.round(atlCorrectionsSingle.easting[0],1)])
                        mat_northingCorrection = np.array([np.round(atlCorrectionsSingle.northing[0],1)])
                        mat_verticalCorrection = np.array([np.round(atlCorrectionsSingle.z[0],1)])
                        mat_crossTrackCorrection = np.array([np.round(atlCorrectionsSingle.crossTrack[0],0)])
                        mat_alongTrackCorrection = np.array([np.round(atlCorrectionsSingle.alongTrack[0],0)])        
                        mat_mae = np.array([np.round(atlCorrectionsSingle.mae[0],2)])
                        mat_rmse = np.array([np.round(atlCorrectionsSingle.rmse[0],2)])
                        mat_me = np.array([np.round(atlCorrectionsSingle.me[0],2)])                      
                        mat_measX_raster = np.c_[atlCorrectionsSingle.measX_raster]
                        mat_measY_raster = np.c_[atlCorrectionsSingle.measY_raster]
                        mat_truthX_raster = np.c_[atlCorrectionsSingle.truthX_raster]
                        mat_truthY_raster = np.c_[atlCorrectionsSingle.truthY_raster]
            
                        # Concatenate .mat file data
                        matData = [mat_atl03FileName, mat_gtNum, mat_trackDirection, \
                                   mat_eastingCorrection, mat_northingCorrection, mat_verticalCorrection, \
                                   mat_crossTrackCorrection, mat_alongTrackCorrection, \
                                   mat_mae, mat_rmse, mat_me, \
                                   mat_measX_raster, mat_measY_raster, \
                                   mat_truthX_raster, mat_truthY_raster]
                    
                        # Make .mat file variable names
                        matNames = ['atl03FileName','gtNum','trackDirection', \
                                    'eastingCorrection','northingCorrection','verticalCorrection', \
                                    'crossTrackCorrection','alongTrackCorrection', \
                                    'meanAbsError','rmse','meanError', \
                                    'icesat2_corrected_utmn_raster','icesat2_corrected_z_raster', \
                                    'reference_utmn_raster','reference_z_raster']
                        
                        # Write output .mat file
                        try:
                            writeLog('Writing output .mat file...\n', logFileID)
                            write_mat(matFilePath, matData, matNames)
                        except:
                            pass
                        # endTry
                        
                    # endIf
                # endFor
                
                # Save output pickle file with all data
                if(createATL03PklFile):
                    try:
                        writeLog('   Writing ATL03 .pkl file...\n', logFileID)
                        pklOutName = atl03Data[0].atl03FileName + '_data.pkl'
                        pklOutPath = os.path.normpath(outFilePath + '/' + pklOutName)
                        saveVars = [atl03Data, atl03DF_all, atl08Data, atlTruthDataFiltered,
                                    atlCorrections, gtNumsGood, beamNumsGood, beamStrengthGood]
                        with open(pklOutPath, 'wb') as f:
                            pkl.dump(saveVars, f)
                        # endWith
                    except:
                        writeLog('   WARNING: Could not write .pkl file.\n', logFileID)
                        pass
                    # endTry
                # endIf
                
                if(len(gtNums)==0):
                    
                    # Reset status bar and print warning
                    statusBar['value'] = 0
                    writeLog('\nWARNING: No Ground Track Selected.\n', logFileID)
                    
                # endIf
                        
                # Load ATL03 plot info
                if(atl03Data):
                    loadAtl03_info()
                    if(atl03Data[0].dataIsMapped):
                        loadAtl08_info()
                    # endIf
                # endIf
                
                # End timer
                timeEnd = runTime.time()
                timeElapsedTotal = timeEnd - timeStart
                timeElapsedMin = np.floor(timeElapsedTotal / 60)
                timeElapsedSec = timeElapsedTotal % 60
                
                # Print completion message
                writeLog('RUN COMPLETE (Total Run Time = %d min %d sec).\n' % (timeElapsedMin, timeElapsedSec), logFileID)
            
            # endFor
            
        except:
            
            # Enable run button
            RunButton.config(state=tk.NORMAL)
        
            # Get error message
            errorMessage = sys.exc_info()[0]
            
            # Print error message
            writeLog('ERROR:', logFileID)
            writeLog('%s\n\n' %errorMessage, logFileID)
            writeLog('Could not process data. Please check inputs.', logFileID)
                        
        #endTry
        
        # Get current date/time
        currentDateTime = (datetime.now()).strftime('%m/%d/%Y %H:%M:%S')
        writeLog('', logFileID)
        writeLog('-------------------------------------', logFileID)
        writeLog('Run Completed On: %s' %currentDateTime, logFileID)
        writeLog('-------------------------------------', logFileID)
        writeLog('', logFileID)
          
        # Close .log file
        logFileID.close()
        
        # Rename .log file
        logNameNew = os.path.splitext(os.path.basename(atl03FilePath))[0] + '_log.txt'
        logFilePathNew = os.path.normpath(outFilePath + '\\' + logNameNew)
        if(os.path.exists(logFilePathNew)):
            os.remove(logFilePathNew)
        # endIf
        os.rename(logFilePath, logFilePathNew)
        
    else:
        
        if(not atl03FileExists):
            messagebox.showinfo('Error','Please select valid input file(s).')
        # endIf
        
        if(not outPathExists):
            messagebox.showinfo('Error','Please select valid output directory.')
        # endIf
        
        if(useTruthSection and (not truthDataExists)):
            messagebox.showinfo('Error','Please select valid reference data file/directory.')
        # endIf
        
    # endIf
    
    # Enable run button
    RunButton.config(state=tk.NORMAL)
    window.update()

# endDef
    
### MEASUREMENT ERROR Panel label
runButtonLabelframe = tk.LabelFrame(tab1, width=545, height=60, text='', font=('Arial Bold', 14))
runButtonLabelframe.place(x=580, y=400)

# Run Button
RunButton = tk.Button(runButtonLabelframe, text='RUN', font=('Arial Bold', 16), width = 18, command=runAtl03) 
RunButton.place(x=50, y=5)

### Make status bar
lbl = tk.Label(runButtonLabelframe, text='Progress:', font=('Arial Bold', 10))
lbl.place(x=320, y=2)
statusBar = Progressbar(runButtonLabelframe, length=190)
statusBar['value'] = 0
statusBar.place(x=320, y=25)


###############################################################################
#
# TAB 2: STATS
#
###############################################################################


### Panel label
statsLabelframeRH = tk.LabelFrame(tab2, width=1110, height=450, text='ICESat-2 Stats', font=('Arial Bold', 14))
statsLabelframeRH.place(x=15, y=10)

def segmentByCallbackRH(event):
    
    # Get  filter value
    segmentByChoice = segmentByBoxRH.get()
    
    # Set increment units text next to increment entry box on GUI
    if('   ' in segmentByChoice):
        incrementTextRH.config(text = '')
    elif('Segment ID' in segmentByChoice):
        incrementTextRH.config(text = '')
    elif('Time (sec)' in segmentByChoice):
        incrementTextRH.config(text = 'sec')
    elif('Latitude (deg)' in segmentByChoice):
        incrementTextRH.config(text = 'deg')
    elif('UTM Northing (m)' in segmentByChoice):
        incrementTextRH.config(text = 'm')
    # endIf
    
    # Update GUI window
    window.update()
    
# endDef

# Function to compute (and export) relative height stats
def computeStatsRH(outputCsvName, verbose=False):
    
    # Only execute if there are selected stats to export
    if(len(statsNumList)>0):
    
        global statsRH_DF
        
        # Refresh stats
        statsRH_DF = []
        
        # Update status bar
        statsStatusBarRH['value'] = 0
        statuslblRH.config(text='Loading')
        
        # Disable stats button
        statsButtonRH.config(state=tk.DISABLED)
        
        # Update window
        window.update()
        
        # Get status inputs
        segmentBy = segmentByBoxRH.get()
        increment = incrementBoxRH.get()
        yval = heightBoxRH.get()
        
        # Set height tag for column name
        if('hae' in yval.lower()):
            ht_tag = 'HAE'
        else:
            ht_tag = 'MSL'
        # endIf
        
        # Get file numbers to plot
        indsToPlotTuple = [gtNumPlotBoxRH.current()]
        
        # Test if inputs are valid
        segmentByValid = segmentBy!=''
        incrementValid = increment!=''
        indsToPlotTupleValid = indsToPlotTuple!=()
        atl03DF_allValid = len(atl03DF_all)>=1
        
        # Test if all inputs are valid
        allValid = segmentByValid and incrementValid and indsToPlotTupleValid and atl03DF_allValid
        
        # Continue code if all inputs are valid, else send message box error
        if(allValid):
            
            if(len(indsToPlotTuple)==1):
            
                # Get correct data frame to use from user
                dfNum = indsToPlotTuple[0]
                atl03DF = atl03DF_all.copy()
                atl03DF = atl03DF[dfNum]
                
                # Get the segment key to bin data by            
                if('Segment ID' == segmentBy):
                    segmentKey = 'Segment ID'
                elif('Time (sec)' == segmentBy):
                    segmentKey = 'Time (sec)'
                elif('Latitude (deg)' == segmentBy):
                    segmentKey = 'Latitude (deg)'
                elif('UTM Northing (m)' == segmentBy):
                    segmentKey = 'UTM Northing (m)'
                elif('Along-Track (m)' == segmentBy):
                    segmentKey = 'Along-Track (m)'
                # endIf
                
                # Convert increment to float
                increment = float(increment)
                
                # Create aggregate list for binning function
                agg_listAll = np.array(['ATL03;atl03_ground_count;' + yval + ';count;[1]',
                                        'ATL03;atl03_ground_count_unique;' + yval + ';count_unique;[1]',
                                        'ATL03;atl03_ground_range (m ' + ht_tag + ');' + yval + ';range;[1]',
                                        'ATL03;atl03_ground_median (m ' + ht_tag + ');' + yval + ';median;[1]',
                                        'ATL03;atl03_ground_mode (m ' + ht_tag + ');' + yval + ';mode;[1]',
                                        'ATL03;atl03_ground_mean (m ' + ht_tag + ');' + yval + ';mean;[1]',
                                        'ATL03;atl03_ground_std (m ' + ht_tag + ');' + yval + ';std;[1]',
                                        'ATL03;atl03_ground_min (m ' + ht_tag + ');' + yval + ';min;[1]',
                                        'ATL03;atl03_ground_max (m ' + ht_tag + ');' + yval + ';max;[1]',
                                        'ATL03;atl03_ground_max98 (m ' + ht_tag + ');' + yval + ';max98;[1]',
                                        'ATL03;atl03_ground_max75 (m ' + ht_tag + ');' + yval + ';max75;[1]',
                                        'ATL03;atl03_ground_max50 (m ' + ht_tag + ');' + yval + ';max50;[1]',
                                        'ATL03;atl03_ground_max25 (m ' + ht_tag + ');' + yval + ';max25;[1]',
                                        'ATL03;atl03_ground_radiometry;' + yval + ';radiometry;[1]'
                                        'ATL03;atl03_all_canopy_count;' + yval + ';count;[2,3]',
                                        'ATL03;atl03_all_canopy_count_unique;' + yval + ';count_unique;[2,3]',
                                        'ATL03;atl03_all_canopy_range (m ' + ht_tag + ');' + yval + ';range;[2,3]',
                                        'ATL03;atl03_all_canopy_median (m ' + ht_tag + ');' + yval + ';median;[2,3]',
                                        'ATL03;atl03_all_canopy_mode (m ' + ht_tag + ');' + yval + ';mode;[2,3]',
                                        'ATL03;atl03_all_canopy_mean (m ' + ht_tag + ');' + yval + ';mean;[2,3]',
                                        'ATL03;atl03_all_canopy_std (m ' + ht_tag + ');' + yval + ';std;[2,3]',
                                        'ATL03;atl03_all_canopy_min (m ' + ht_tag + ');' + yval + ';min;[2,3]',
                                        'ATL03;atl03_all_canopy_max (m ' + ht_tag + ');' + yval + ';max;[2,3]',
                                        'ATL03;atl03_all_canopy_max98 (m ' + ht_tag + ');' + yval + ';max98;[2,3]',
                                        'ATL03;atl03_all_canopy_max75 (m ' + ht_tag + ');' + yval + ';max75;[2,3]',
                                        'ATL03;atl03_all_canopy_max50 (m ' + ht_tag + ');' + yval + ';max50;[2,3]',
                                        'ATL03;atl03_all_canopy_max25 (m ' + ht_tag + ');' + yval + ';max25;[2,3]',
                                        'ATL03;atl03_all_canopy_relative_height (m ' + ht_tag + ');' + yval + ';rh_canopy;[2,3]',
                                        'ATL03;atl03_all_canopy_radiometry;' + yval + ';radiometry;[2,3]'])
                                
                agg_list = agg_listAll[statsNumList]
                
                try:
                    
                    # Bin data into dataframe
                    atl03DF_binned = get_bin_df(atl03DF, segmentKey, increment, agg_list)
                        
                    # Pull subset of data into smaller dataframe
                    statsDF_orig = atl03DF_binned.copy()
                    statsDF_orig = statsDF_orig.drop(statsDF_orig.columns[[1,3,4]], axis=1)
                    
                    # Rename 'beg_id' column name to Segment By name
                    newColName = 'segment_start_' + segmentBy.lower()
                    statsDF_orig.rename(columns={'beg_id':newColName}, inplace=True)
                    
                    # Rename 'end_id' column name to Segment By name
                    newColName = 'segment_end_' + segmentBy.lower()
                    statsDF_orig.rename(columns={'end_id':newColName}, inplace=True)
                    
                    # Interpolate stats dataframe midpoints
                    statsRH_DF = interpStats(atl03DF, statsDF_orig, ht_tag, rh=True)
                    
                    # Ask user for name if it doesn't exist
                    if(not outputCsvName):
    
                        # Prompt user to save CSV file somewhere        
                        files = [('CSV File', '*.csv')] 
                        initialDir = 'C:/'
                        outputCsvObj = asksaveasfile(initialdir = initialDir, \
                                                      initialfile = 'outputData', title = 'Save CSV File', \
                                                      filetypes = files, defaultextension = files)
                        outputCsvName = outputCsvObj.name
                    # endIf
                    
                    # Write stats to csv file
                    writeCsv(statsRH_DF, outputCsvName, verbose)
                    
                    # Update status bar
                    statsStatusBarRH['value'] = 100
                    statuslblRH.config(text='Complete')
        
                    # Update window
                    window.update()
                
                except:
                    
                    messagebox.showinfo('Error','Could not compute stats. Please check inputs.')
                    
                # endTry
                
            else:
                messagebox.showinfo('Error','Please select only 1 file to compute stats.') 
            # endIf
        else:
            messagebox.showinfo('Error','Missing data to compute stats.') 
        # endIf
        
        # Enable button
        statsButtonRH.config(state=tk.NORMAL)
        window.update()
    
    # endIf
    
# endDef

# Function to compute relative height stats
def computeStatsRH_callback():
    
    # Call compute stats function
    outputCsvName = False
    
    try:
        if(len(atl03Data)>0):
            if(atl03Data[0].dataIsMapped):
                computeStatsRH(outputCsvName, verbose=True)
            else:
                messagebox.showinfo('Error','ATL03 data needs to be run with an ATL08 file to create stats.')
            # endIf
        # endIf
    except:
        pass
    # endTry
    
# endDef

# Function to move stats variables over from one listbox to another
def moveToExport():
            
    global statsNumList 
    
    # Get selection indices
    statsNums = statsListBoxRH.curselection()
    
    # For each index, insert into Stats to Export listbox
    try:
        if(np.size(statsNums)>0):
            
            currDataIsMapped = False
            try:
                # Get file numbers to plot
                indsToPlotTuple = [gtNumPlotBoxRH.current()]
                currDataIsMapped = atl03Data[indsToPlotTuple].dataIsMapped
            except:
                pass
            # endTry
            
            if((not atl08FileExists) and (not currDataIsMapped)):
                messagebox.showinfo('WARNING','An ATL08 file is required to calculate and export stats.')
            # endIf
    
            for i in statsNums:
                statsParam = statsListBoxRH.get(i)
                statsCompListBoxRH.insert('end', statsParam)
                window.update()
            # endFor
            statsListBoxRH.selection_clear(0, 'end')
            statsNumList = np.append(statsNumList,statsNums)
        # endIf
    except:
        pass
    # endTry
    
# endDef
    
# Function to delete stats variables from listbox
def deleteStatsRH():
            
    global statsNumList 
    
    # Get selection indices
    statsNums = statsCompListBoxRH.curselection()
    
    # For each index, insert into Stats to Export listbox
    if(np.size(statsNums)>0):
        for i in reversed(statsNums):
            statsCompListBoxRH.delete(i)
        # endFor
        window.update()
        statsNumList = np.delete(statsNumList,statsNums)
    # endIf
    
# endDef

# Function called when GT num is selected
def gtNumPlotCallback(event):
    
    # Refresh stats
    refreshStats()
 
# endDef
    
# Function to load new pickle file
def loadNewPklFile():
    
    # Check output file path
    startDir = outPath_textBox.get().strip()
    if('' == startDir or '.' == startDir):
        startDir = cwd
    else:
        outPathExists = os.path.isdir(startDir)
        if(not outPathExists):
            os.mkdir(outFilePath)
            startDir = cwd
        # endIf
    # endIf
    
    # Ask user for input .pkl files
    pklFile = filedialog.askopenfilename(initialdir = startDir, title = 'Select File to Load', filetypes = [('*data.pkl files','*data.pkl')], multiple=False)

    # If the user selection is not empty, then continue
    try:
        if(pklFile != '.'):
                
            global atl03Data, atl03DF_all, atl08Data, atlTruthDataFiltered, \
                   atlCorrections, gtNumsGood, beamNumsGood, beamStrengthGood
            
            # Update status bar
            statsStatusBarRH['value'] = 0
            statuslblRH.config(text='   Status:')
            
            # Delete contents in Selected Stats listbox
            statsCompListBoxRH.delete(0,'end')
        
            # Load pickle file
            with open(pklFile, 'rb') as f:
                atl03Data, atl03DF_all, atl08Data, atlTruthDataFiltered, atlCorrections, gtNumsGood, beamNumsGood, beamStrengthGood = pkl.load(f)
            # endWith
            
            # Update ATL03/ATL08 entries
            loadAtl03_info()
            if(atl03Data[0].dataIsMapped):
                loadAtl08_info()
            # endIf
        # endIf
    except:
        pass
    # endTry
    
# endDef
    
### File Name Entry Box
lbl = tk.Label(statsLabelframeRH, text='File Name:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=10)
fileName_textBoxRH = tk.Entry(statsLabelframeRH, width=50)
fileName_textBoxRH.place(x=10, y=40)

### GT Num Plot Combo Box
lbl = tk.Label(statsLabelframeRH, text='Ground Track:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=350, y=10)
gtNumPlotBoxRH = Combobox(statsLabelframeRH, width=26)
gtNumPlotBoxRH.place(x=350, y=40)
gtNumPlotBoxRH.bind("<<ComboboxSelected>>", gtNumPlotCallback)

### Load New File Button
atl03BrowseButtonRH = tk.Button(statsLabelframeRH, text='Load File', font=('Arial Bold', 16), width=12, command=loadNewPklFile) 
atl03BrowseButtonRH.place(x=560, y=25)

### Segment By section
lbl = tk.Label(statsLabelframeRH, text='Segment Stats by:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=75)
segmentByBoxRH = ttk.Combobox(statsLabelframeRH, width=20)
segmentByBoxRH.place(x=10, y=105)
segmentByTypesRH = ['Time (sec)','Latitude (deg)','UTM Northing (m)','Along-Track (m)']
segmentByBoxRH['values'] = segmentByTypesRH
segmentByBoxRH.current(0)
segmentByBoxRH.bind("<<ComboboxSelected>>", segmentByCallbackRH)

### Increment section
lbl = tk.Label(statsLabelframeRH, text='Increment:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=190, y=75)
incrementBoxRH = tk.Entry(statsLabelframeRH, width=14)
incrementBoxRH.insert(0,'1.0')
incrementBoxRH.place(x=190, y=105)
incrementTextRH = tk.Label(statsLabelframeRH, text='sec', font=('Arial', 12), anchor = 'w', justify='left')
incrementTextRH.place(x=280, y=102)

### Height Selection section
lbl = tk.Label(statsLabelframeRH, text='Height Type:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=350, y=75)
heightBoxRH = ttk.Combobox(statsLabelframeRH, width=26)
heightBoxTypesRH = ['Height (m HAE)','Height (m MSL)']
heightBoxRH['values'] = heightBoxTypesRH
heightBoxRH.current(0)
heightBoxRH.place(x=350, y=105)
    
### Auto Create Output Stats File Check Box
createStatsFileChkState = tk.BooleanVar()
createStatsFileChkState.set(True)
statsFileText = 'Auto Save Selected Stats to .csv File After PhoREAL Runs'
createStatsFile_checkBox = tk.Checkbutton(statsLabelframeRH, text=statsFileText, 
                                          font=('Arial', 12), wraplength=180, 
                                          justify='left', var=createStatsFileChkState) 
createStatsFile_checkBox.place(x=560, y=85)

### Available Stats List Box
lbl = tk.Label(statsLabelframeRH, text='Available Stats:', font=('Arial Bold', 12), anchor='w', justify='left')
lbl.place(x=10, y=150)
statsListBoxRH = tk.Listbox(statsLabelframeRH, width=36, height=14, selectmode='multiple') 
statsListBoxRH.place(x=10, y=180)
yScroll1 = tk.Scrollbar(statsLabelframeRH, orient=tk.VERTICAL)
yScroll1.config(command=statsListBoxRH.yview)
yScroll1.place(x=230, y=180)
statsListBoxRH.config(yscrollcommand=yScroll1)
statsColNamesRH = ['Ground Height Count','Ground Height Count (Unique)','Ground Height Range',
                   'Ground Height Median','Ground Height Mode','Ground Height Mean',
                   'Ground Height Standard Deviation','Ground Height Min','Ground Height Max',
                   'Ground Height Max98','Ground Height Max75','Ground Height Max50',
                   'Ground Height Max25','Ground Radiometry',
                   'All Canopy Height Count','All Canopy Height Count (Unique)','All Canopy Height Range',
                   'All Canopy Height Median','All Canopy Height Mode','All Canopy Height Mean',
                   'All Canopy Height Standard Deviation','All Canopy Height Min','All Canopy Height Max',
                   'All Canopy Height Max98','All Canopy Height Max75','All Canopy Height Max50',
                   'All Canopy Height Max25','All Canopy Relative Height','All Canopy Radiometry']
statsListBoxRH.delete(0, 'end')
statsListBoxRH.insert(0, *statsColNamesRH)
    
### Move to Export Stats Button
MoveToYaxisButton = tk.Button(statsLabelframeRH, text='>>', font=('Arial Bold', 12), width = 4, command=moveToExport) 
MoveToYaxisButton.place(x=265, y=250)

### Stats to Compute Listbox
lbl = tk.Label(statsLabelframeRH, text='Selected Stats:', font=('Arial Bold', 12), anchor='w', justify='left')
lbl.place(x=350, y=150)
statsCompListBoxRH = tk.Listbox(statsLabelframeRH, width=36, height=14, selectmode='multiple') 
statsCompListBoxRH.place(x=350, y=180)
yScroll2 = tk.Scrollbar(statsLabelframeRH, orient=tk.VERTICAL)
yScroll2.config(command=statsListBoxRH.yview)
yScroll2.place(x=570, y=180)
statsCompListBoxRH.config(yscrollcommand=yScroll2)

### Export Stats button
statsButtonRH = tk.Button(statsLabelframeRH, text='Export Stats', font=('Arial Bold', 16), width = 13, command=computeStatsRH_callback) 
statsButtonRH.place(x=610, y=210)

### Delete Stats button
removeStatsButtonRH = tk.Button(statsLabelframeRH, text='Remove Stats', font=('Arial Bold', 16), width = 13, command=deleteStatsRH) 
removeStatsButtonRH.place(x=610, y=280)

### Compute Stats Status Bar
statuslblRH = tk.Label(statsLabelframeRH, text='   Status:', font=('Arial Bold', 10))
statuslblRH.place(x=695, y=180)
statsStatusRH = int()
statsStatusBarRH = Progressbar(statsLabelframeRH, variable=statsStatusRH, length=20)
statsStatusBarRH['value'] = 0
statsStatusBarRH.place(x=765, y=180)

statsInfo = \
'User Instructions:\n' \
'-------------------------------------------------\n' \
'Stats can be exported in three ways:\n' \
'\n' \
'1) Execute PhoREAL with an ATL03\n' \
'    and ATL08 file, then select the stats\n' \
'    to export.\n' \
'\n' \
'2) Load a file previously executed with\n' \
'    PhoREAL (*_data.pkl), then select\n' \
'    the stats to export.\n' \
'\n' \
'3) To export a stats .csv file after\n' \
'    each run in batch mode: select the\n' \
'    stats to export, check the box at left\n' \
'    to auto save a stats .csv file after\n' \
'    each run, then execute PhoREAL\n' \
'    with ATL03 and ATL08 file(s).\n' \
'    Output .csv files will be stored in\n' \
'    the output file directory with the\n' \
'    name <baseFileName>_stats.csv.'

lbl = tk.Label(statsLabelframeRH, text=statsInfo, font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=820, y=10)


###############################################################################
#
# TAB 3: PLOT DATA
#
###############################################################################

### Panel label
labelframe = tk.LabelFrame(tab3, width=545, height=270, text='Plot ICESat-2 ATL03 Data', font=('Arial Bold', 14))
labelframe.place(x=15, y=10)

## Plot text
#lbl = tk.Label(tab3, text='ATL03 Plotting Options:', font=('Arial Bold', 12), anchor = 'w', justify='left')
#lbl.place(x=30, y=50)

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
    
    # Update status bar
    statsStatusBar['value'] = 0
    statuslbl.config(text='')

    # Clear Add Stats listbox in Stats Section
    addStatsTuple = ('')
    addStatsBox['values'] = addStatsTuple
    addStatsBox.delete(0,'end')
    
    # Update window
    window.update()
                
# endDef  
    
# Filter Number Checkboxes
filter0ChkState = tk.BooleanVar()
filter0ChkState.set(False)
filter0_checkBox = tk.Checkbutton(tab3, text = '0', font=('Arial', 12), var = filter0ChkState) 

filter1ChkState = tk.BooleanVar()
filter1ChkState.set(False)
filter1_checkBox = tk.Checkbutton(tab3, text = '1', font=('Arial', 12), var = filter1ChkState) 

filter2ChkState = tk.BooleanVar()
filter2ChkState.set(False)
filter2_checkBox = tk.Checkbutton(tab3, text = '2', font=('Arial', 12), var = filter2ChkState) 

filter3ChkState = tk.BooleanVar()
filter3ChkState.set(False)
filter3_checkBox = tk.Checkbutton(tab3, text = '3', font=('Arial', 12), var = filter3ChkState) 

filter4ChkState = tk.BooleanVar()
filter4ChkState.set(False)
filter4_checkBox = tk.Checkbutton(tab3, text = '4', font=('Arial', 12), var = filter4ChkState) 
     
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
        filter0_checkBox.place(x=30, y=230) 
        filter0_checkBox.config(text = 'Unclassified')
        filter0ChkState.set(False)
        filter1_checkBox.place(x=160, y=230)
        filter1_checkBox.config(text = 'Ground')
        filter1ChkState.set(False)
        filter2_checkBox.place(x=250, y=230)
        filter2_checkBox.config(text = 'Canopy') 
        filter2ChkState.set(False)
        filter3_checkBox.place(x=350, y=230)  
        filter3_checkBox.config(text = 'Top of Canopy') 
        filter3ChkState.set(False)
        filter4_checkBox.place_forget()
        filter4ChkState.set(False)
        
    elif('signal' in filterChoice.lower()):
        
        # Filter Number Checkboxes
        filter0_checkBox.place(x=30, y=230) 
        filter0_checkBox.config(text = '0')
        filter0ChkState.set(False)
        filter1_checkBox.place(x=80, y=230)  
        filter1_checkBox.config(text = '1')
        filter1ChkState.set(False)
        filter2_checkBox.place(x=130, y=230)
        filter2_checkBox.config(text = '2') 
        filter2ChkState.set(False)        
        filter3_checkBox.place(x=180, y=230)
        filter3_checkBox.config(text = '3')
        filter3ChkState.set(False)
        filter4_checkBox.place(x=230, y=230)
        filter4_checkBox.config(text = '4')
        filter4ChkState.set(False)
            
    # endIf
# endDef
 
# File Name Entry Box
lbl = tk.Label(tab3, text='File Name:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=50)
fileName_textBox = tk.Entry(tab3, width=50)
fileName_textBox.place(x=120, y=52)

# Load New File Button
atl03BrowseButton = tk.Button(tab3, text='Load', font=('Arial Bold', 16), width=6, command=loadNewPklFile) 
atl03BrowseButton.place(x=450, y=40)

# GT Num Plot Combo Box
lbl = tk.Label(tab3, text='Ground Track:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=80)
gtNumPlotBox = Combobox(tab3, width=26)
gtNumPlotBox.place(x= 140, y = 82)
gtNumPlotBox.bind("<<ComboboxSelected>>", gtNumPlotCallback)

# X Axis Combo Box
lbl = tk.Label(tab3, text='X Axis:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=110)
xValsBox = Combobox(tab3)
xValsBox.place(x= 90, y = 112)
xValsBox.bind("<<ComboboxSelected>>", xAxisCallback)

# X Label Entry Box
lbl = tk.Label(tab3, text='X Label:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=260, y=110)
xlabel_textBox = tk.Entry(tab3, width=33)
xlabel_textBox.place(x=330, y=112)

# Y Axis Combo Box
lbl = tk.Label(tab3, text='Y Axis:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=140)
yValsBox = Combobox(tab3)
yValsBox.place(x= 90, y = 142)
yValsBox.bind("<<ComboboxSelected>>", yAxisCallback)

# Y Label Entry Box
lbl = tk.Label(tab3, text='Y Label:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=260, y=140)
ylabel_textBox = tk.Entry(tab3, width=33)
ylabel_textBox.place(x=330, y=142)

# Filter On Combo Box
lbl = tk.Label(tab3, text='Filter On:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=190)
filterBox = Combobox(tab3, width=27)
filterBox.place(x= 110, y = 192)
filterBox.bind("<<ComboboxSelected>>", filterChoiceCallback)

# Itialize lists
plotList = ('time', 'deltaTime', 'lat', 'lon', 'easting', 'northing', \
            'crossTrack', 'alongTrack', 'z', 'zMsl', 'classification', 'signalConf')
filterState = np.array([0,1,2,3,4])

# Plot Button Callback
def plotAtl03():
    
    # Try plot code
    try:
        
        # Get
        gtNumToPlot = gtNumPlotBox.current()
        
        # Get x,y combo box number selections
        xVarNum = xValsBox.current()
        yVarNum = yValsBox.current()
        
        # Get x,y combo bxx text selections
        xData = eval('atl03Data[' + str(gtNumToPlot) + '].' + plotList[xVarNum])
        yData = eval('atl03Data[' + str(gtNumToPlot) + '].' + plotList[yVarNum])
        fileName = eval('atl03Data[' + str(gtNumToPlot) + '].atl03FileName')
        gtNum = eval('atl03Data[' + str(gtNumToPlot) + '].gtNum')
        
        # Get labels
        xLabel = xlabel_textBox.get()
        yLabel = ylabel_textBox.get()
        title = fileName + ' (' + gtNum + ')'
        
        # Get Filter data type (classification or signal confidence) and filter numbers
        filterChoice = filterBox.get()
        if('  ' in filterChoice.lower()):
            filterType = []
            filterData = []
            filterNum = []
        elif('class' in filterChoice.lower()):
            filterType = filterChoice
            filterData = eval('atl03Data[' + str(gtNumToPlot) + '].classification')
            filterTF = [filter0ChkState.get(), filter1ChkState.get(), filter2ChkState.get(), filter3ChkState.get(), filter4ChkState.get()]    
            filterNum = filterState[filterTF]
        elif('signal' in filterChoice.lower()):
            filterType = filterChoice
            filterData = eval('atl03Data[' + str(gtNumToPlot) + '].signalConf')
            filterTF = [filter0ChkState.get(), filter1ChkState.get(), filter2ChkState.get(), filter3ChkState.get(), filter4ChkState.get()]    
            filterNum = filterState[filterTF]
        # endIf
        
    
        # Get output path
        outPath = outPath_textBox.get().strip()
        
        # Get original plot title
        origTitle = title
        
        # Get gtNum
        gtNum = gtNumPlotBox.current()
        
        # Call getPlot function
        getPlot(xData, yData, xLabel, yLabel, title, outPath, origTitle,
                filterType, filterData, filterNum)
    
    except:
        
        print('Cannot plot data. Please check inputs')
        
    # endTry
# endDef
    
# Plot ATL03 Button
btn = tk.Button(tab3, text='Plot', font=('Arial Bold', 16), width = 15, command=plotAtl03) 
btn.place(x=330, y=178)  


###############################################################################
#
# TAB 3: PLOT DATA - STATS
#
###############################################################################

# Function to determine increment units and set on GUI
def segmentByCallback(event):
    
    # Get  filter value
    segmentByChoice = segmentByBox.get()
    
    # Set increment units text next to increment entry box on GUI
    if('   ' in segmentByChoice):
        incrementText.config(text = '')
    elif('Segment ID' in segmentByChoice):
        incrementText.config(text = '')
    elif('Time (sec)' in segmentByChoice):
        incrementText.config(text = 'sec')
    elif('Latitude (deg)' in segmentByChoice):
        incrementText.config(text = 'deg')
    elif('UTM Northing (m)' in segmentByChoice):
        incrementText.config(text = 'm')
    elif('Along-Track (m)' in segmentByChoice):
        incrementText.config(text = 'm')
    # endIf
    
    # Update GUI window
    window.update()
    
# endDef

# Interpolate stats dataframe midpoints
def interpStats(atl03DF, statsDF_orig, ht_tag, rh=False):
        
    # Get binned X data - start
    binned_x_colName = statsDF_orig.columns[0]
    binned_x_start = statsDF_orig[binned_x_colName]
    
    # Get binned X data - end
    binned_x_colName = statsDF_orig.columns[1]
    binned_x_end = statsDF_orig[binned_x_colName]
    
    # Get status inputs
    if(rh):
        segmentBy = segmentByBoxRH.get()
    else:
        segmentBy = segmentByBox.get()
    # endIf
    
    # Get the segment key to bin data by            
    if('Segment ID' == segmentBy):
        segmentKey = 'Segment ID'
    elif('Time (sec)' == segmentBy):
        segmentKey = 'Time (sec)'
    elif('Latitude (deg)' == segmentBy):
        segmentKey = 'Latitude (deg)'
    elif('UTM Northing (m)' == segmentBy):
        segmentKey = 'UTM Northing (m)'
    elif('Along-Track (m)' == segmentBy):
        segmentKey = 'Along-Track (m)'
    # endIf
            
    # Drop duplicate rows in segmentKey column
    atl03DF_in = atl03DF.drop_duplicates(subset=segmentKey, keep='first', inplace=False, ignore_index=False)

    # Get full X data
    full_x = atl03DF_in[segmentKey]
    
    # Interpolate delta_time
    full_delta_time = atl03DF_in['Delta Time (sec)']
    f1 = interpolate.interp1d(full_x, full_delta_time, kind='linear', fill_value='extrapolate')
    binned_delta_time_start = f1(binned_x_start)
    binned_delta_time_end = f1(binned_x_end)
    
    # Interpolate time
    full_time = atl03DF_in['Time (sec)']
    f2 = interpolate.interp1d(full_x, full_time, kind='linear', fill_value='extrapolate')
    binned_time_start = f2(binned_x_start)
    binned_time_end = f2(binned_x_end)

    # Interpolate lat_ph
    full_lat_ph = atl03DF_in['Latitude (deg)']
    f3 = interpolate.interp1d(full_x, full_lat_ph, kind='linear', fill_value='extrapolate')
    binned_lat_ph_start = f3(binned_x_start)
    binned_lat_ph_end = f3(binned_x_end)
    
    # Interpolate lon_ph
    full_lon_ph = atl03DF_in['Longitude (deg)']
    f4 = interpolate.interp1d(full_x, full_lon_ph, kind='linear', fill_value='extrapolate')
    binned_lon_ph_start = f4(binned_x_start)
    binned_lon_ph_end = f4(binned_x_end)
    
    # Interpolate easting
    full_easting = atl03DF_in['UTM Easting (m)']
    f5 = interpolate.interp1d(full_x, full_easting, kind='linear', fill_value='extrapolate')
    binned_easting_start = f5(binned_x_start)
    binned_easting_end = f5(binned_x_end)
    
    # Interpolate northing
    full_northing = atl03DF_in['UTM Northing (m)']
    f6 = interpolate.interp1d(full_x, full_northing, kind='linear', fill_value='extrapolate')
    binned_northing_start = f6(binned_x_start)
    binned_northing_end = f6(binned_x_end)
    
    # Interpolate crossTrack
    full_crossTrack = atl03DF_in['Cross-Track (m)']
    f7 = interpolate.interp1d(full_x, full_crossTrack, kind='linear', fill_value='extrapolate')
    binned_crossTrack_start = f7(binned_x_start)
    binned_crossTrack_end = f7(binned_x_end)
    
    # Interpolate alongTrack
    full_alongTrack = atl03DF_in['Along-Track (m)']
    f8 = interpolate.interp1d(full_x, full_alongTrack, kind='linear', fill_value='extrapolate')
    binned_alongTrack_start = f8(binned_x_start)
    binned_alongTrack_end = f8(binned_x_end)
    
    # Put columns into one array                      
    addedArray = np.column_stack([binned_delta_time_start, binned_delta_time_end,
                                  binned_time_start, binned_time_end,
                                  binned_lat_ph_start, binned_lat_ph_end,
                                  binned_lon_ph_start, binned_lon_ph_end,
                                  binned_easting_start, binned_easting_end,
                                  binned_northing_start, binned_northing_end,
                                  binned_crossTrack_start, binned_crossTrack_end,
                                  binned_alongTrack_start, binned_alongTrack_end])
    
    # Convert array into Pandas dataframe                    
    addedDF = pd.DataFrame(data=addedArray, columns=['seg_start_delta_time_interp (sec)',
                                                     'seg_end_delta_time_interp (sec)',
                                                     'seg_start_time_interp (sec)',
                                                     'seg_end_time_interp (sec)',
                                                     'seg_start_lat_interp (deg)',
                                                     'seg_end_lat_interp (deg)',
                                                     'seg_start_lon_interp (deg)',
                                                     'seg_end_lon_interp (deg)',
                                                     'seg_start_easting_interp (m)',
                                                     'seg_end_easting_interp (m)',
                                                     'seg_start_northing_interp (m)',
                                                     'seg_end_northing_interp (m)',
                                                     'seg_start_crossTrack_interp (m)',
                                                     'seg_end_crossTrack_interp (m)',
                                                     'seg_start_alongTrack_interp (m)',
                                                     'seg_end_alongTrack_interp (m)'])
    
    # Concatenate new dataframe onto original one
    statsNewDF = pd.concat([statsDF_orig, addedDF], axis=1)
    
    # Return new dataframe
    return statsNewDF

# endDef
    
# Function to compute stats
def computeStats():
        
    global statsDF
    
    # Refresh stats
    refreshStats()
    
    # Update status bar
    statsStatusBar['value'] = 0
    statuslbl.config(text='Loading')
    
    # Disable stats button
    statsButton.config(state=tk.DISABLED)
    
    # Update window
    window.update()
    
    # Get status inputs
    segmentBy = segmentByBox.get()
    increment = incrementBox.get()
    yval = yValsBox.get()
    
    # Set height tag for column name
    if('hae' in yval.lower()):
        ht_tag = 'HAE'
    else:
        ht_tag = 'MSL'
    # endIf
    
    # Only compute stats when HAE or MSL heights are on the y-axis
    if('height' in yval.lower()):
    
        # Get file numbers to plot
        indsToPlotTuple = [gtNumPlotBox.current()]
        
        # Test if inputs are valid
        segmentByValid = segmentBy!=''
        incrementValid = increment!=''
        indsToPlotTupleValid = indsToPlotTuple!=()
        atl03DF_allValid = len(atl03DF_all)>=1
        
        # Test if all inputs are valid
        allValid = segmentByValid and incrementValid and indsToPlotTupleValid and atl03DF_allValid
    
        # Continue code if all inputs are valid, else send message box error
        if(allValid):
            
            if(len(indsToPlotTuple)==1):
            
                # Get correct data frame to use from user
                dfNum = indsToPlotTuple[0]
                atl03DF = atl03DF_all.copy()
                atl03DF = atl03DF[dfNum]
                
                # Get the segment key to bin data by            
                if('Segment ID' == segmentBy):
                    segmentKey = 'Segment ID'
                elif('Time (sec)' == segmentBy):
                    segmentKey = 'Time (sec)'
                elif('Latitude (deg)' == segmentBy):
                    segmentKey = 'Latitude (deg)'
                elif('UTM Northing (m)' == segmentBy):
                    segmentKey = 'UTM Northing (m)'
                elif('Along-Track (m)' == segmentBy):
                    segmentKey = 'Along-Track (m)'
                # endIf
                
                # Convert increment to float
                increment = float(increment)
                
                # Create aggregate list for binning function
                agg_list = ['ATL03;atl03_ground_min (m ' + ht_tag + ');' + yval + ';min;[1]',
                            'ATL03;atl03_ground_max (m ' + ht_tag + ');' + yval + ';max100;[1]',
                            'ATL03;atl03_ground_median (m ' + ht_tag + ');' + yval + ';median;[1]',
                            'ATL03;atl03_ground_mean (m ' + ht_tag + ');' + yval + ';mean;[1]',
                            'ATL03;atl03_ground_std (m ' + ht_tag + ');' + yval + ';std;[1]',
                            'ATL03;atl03_all_canopy_min (m ' + ht_tag + ');' + yval + ';min;[2,3]',
                            'ATL03;atl03_all_canopy_max (m ' + ht_tag + ');' + yval + ';max100;[2,3]',
                            'ATL03;atl03_all_canopy_median (m ' + ht_tag + ');' + yval + ';median;[2,3]',
                            'ATL03;atl03_all_canopy_mean (m ' + ht_tag + ');' + yval + ';mean;[2,3]',
                            'ATL03;atl03_all_canopy_std (m ' + ht_tag + ');' + yval + ';std;[2,3]',
                            'ATL03;atl03_all_height_min (m ' + ht_tag + ');' + yval + ';min;[1,2,3]',
                            'ATL03;atl03_all_height_max (m ' + ht_tag + ');' + yval + ';max100;[1,2,3]',
                            'ATL03;atl03_all_height_median (m ' + ht_tag + ');' + yval + ';median;[1,2,3]',
                            'ATL03;atl03_all_height_mean (m ' + ht_tag + ');' + yval + ';mean;[1,2,3]',
                            'ATL03;atl03_all_height_std (m ' + ht_tag + ');' + yval + ';std;[1,2,3]']
                
                try:
                
                    # Set Add Stats listbox in Stats Section
                    addStatsBox.set('')
                    window.update()
                    
                    # Bin data into dataframe
                    atl03DF_binned = get_bin_df(atl03DF, segmentKey, increment, agg_list)
                        
                    # Pull subset of data into smaller dataframe
                    statsDF_orig = atl03DF_binned.copy()
                    statsDF_orig = statsDF_orig[statsDF_orig.columns[[0,2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]]]
                    
                    # Rename 'beg_id' column name to Segment By name
                    newColName = 'segment_start_' + segmentBy.lower()
                    statsDF_orig.rename(columns={'beg_id':newColName}, inplace=True)
                    
                    # Rename 'end_id' column name to Segment By name
                    newColName = 'segment_end_' + segmentBy.lower()
                    statsDF_orig.rename(columns={'end_id':newColName}, inplace=True)
                    
                    # Interpolate stats dataframe midpoints
                    statsDF = interpStats(atl03DF, statsDF_orig, ht_tag)
                    
                    # Update status bar
                    statsStatusBar['value'] = 100
                    statuslbl.config(text='Complete')
                    
                    # Set Add Stats listbox in Stats Section
                    addStatsTuple = ('Ground Min','Ground Max','Ground Median','Ground Mean','Ground Mean + 3*Std','Ground Mean - 3*Std',
                            'All Canopy Min','All Canopy Max','All Canopy Median','All Canopy Mean','All Canopy Mean + 3*Std','All Canopy Mean - 3*Std',
                            'All Height Min','All Height Max','All Height Median','All Height Mean','All Height Mean + 3*Std','All Height Mean - 3*Std')
                    addStatsBox['values'] = addStatsTuple
                    addStatsBox.current(0)
        
                    # Update window
                    window.update()
                
                except:
                    
                    messagebox.showinfo('Error','Could not compute stats. Please check inputs.')
                    
                # endTry
                
            else:
                messagebox.showinfo('Error','Please select only 1 file to compute stats.') 
            # endIf
        else:
            messagebox.showinfo('Error','Missing data to compute stats.') 
        # endIf
    
    else:
        messagebox.showinfo('Error','Can only compute stats with HAE or MSL Heights in Y Axis above.') 
    # endIf
    
    # Enable button
    statsButton.config(state=tk.NORMAL)
    window.update()
                
# endDef
    
# Add stats to figure callback
def addStatsCallback():
    
    try:
        
        # Get file numbers to plot
        indsToPlotTuple = [gtNumPlotBox.current()]
            
        # Get x parameter to plot
        xParam = plotList[xValsBox.current()]
        
        # Get y parameter to plot
        yParam = addStatsBox.get()
        
        # Get y-axis plot parameter
        yVar = plotList[yValsBox.current()]
        
        # Get HAE/MSL
        yHt = yValsBox.get()
    
        # Add stats to figure    
        addStatsToPlot(indsToPlotTuple,xParam,yParam,yVar,yHt,statsDF)
    
    except:
        
        messagebox.showinfo('Error','Missing data to plot stats.') 
        
    # endTry
    
# endDef
    
# Write stats callback button
def writeStatsToCsv():
    
    try:
        if(len(statsDF)>0):
            
            # Disable button
            writeStatsButton.config(state=tk.NORMAL)
            window.update()
            
            # Try to write .csv file
            try:
            
                # Get dataframe for file
                csvDF = statsDF
                    
                # Prompt user to save CSV file somewhere        
                files = [('CSV File', '*.csv')] 
                initialDir = 'C:/'
                outputCsvName = asksaveasfile(initialdir = initialDir, \
                                              initialfile = 'outputData', title = 'Save CSV File', \
                                              filetypes = files, defaultextension = files)
            
                # Write CSV file
                writeCsv(csvDF, outputCsvName)
                    
            except:
                
                messagebox.showinfo('Error','Could not write .csv file.')
                
            # end Try
            
            # Enable button
            writeStatsButton.config(state=tk.NORMAL)
            window.update()
        # endIf
    except:
        pass
    # endTry
    
# endDef    
    
# Write CSV Parameter Button Callback
def writeCsv(csvDF, outputCsvName, verbose=True):
    
    # Try to write CSV file
    try:
        
        if(verbose):
            print('Writing to CSV file...')
        # endIf
        csvDF.to_csv(outputCsvName, index=False, line_terminator='\n')
        if(verbose):
            print('File Complete!')
        # endIf
        
        # Send completion message
        if(verbose):
            messagebox.showinfo('Success','File Complete!')
        # endIf
        
    except:
        
        # Send error message
        if(verbose):
            messagebox.showinfo('Error','Could not write .csv file.')
        # endIf
    # endTry
# endDef
    
### Stats Data Panel label
statsLabelframe = tk.LabelFrame(tab3, width=545, height=175, text='Compute and Plot Stats', font=('Arial Bold', 14))
statsLabelframe.place(x=15, y=290)
    
### "Segment by:" section
lbl = tk.Label(statsLabelframe, text='Segment Stats by:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=20)
segmentByBox = ttk.Combobox(statsLabelframe, width=20)
segmentByBox.place(x=10, y=50)
segmentByBox.bind("<<ComboboxSelected>>", segmentByCallback)

### "Increment" section
lbl = tk.Label(statsLabelframe, text='Segment Increment:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=85)
incrementBox = tk.Entry(statsLabelframe, width=14)
incrementBox.place(x=10, y=115)
incrementText = tk.Label(statsLabelframe, text='', font=('Arial', 12), anchor = 'w', justify='left')
incrementText.place(x=100, y=112)
    
### Compute Stats button
statsButton = tk.Button(statsLabelframe, text='Compute Stats', font=('Arial Bold', 16), width = 13, command=computeStats) 
statsButton.place(x=180, y=30)

### Compute Stats Status Bar
statuslbl = tk.Label(statsLabelframe, text='   Status:', font=('Arial Bold', 10))
statuslbl.place(x=270, y=0)
statsStatus = int()
statsStatusBar = Progressbar(statsLabelframe, variable=statsStatus, length=20)
statsStatusBar['value'] = 0
statsStatusBar.place(x=340, y=0)

### "Add Stats to Plot:" section
lbl = tk.Label(statsLabelframe, text='Add Stats to Plot:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=180, y=85)
addStatsBox = ttk.Combobox(statsLabelframe, width=26)
addStatsBox.place(x=180, y=115)
#addStatsBox.bind("<<ComboboxSelected>>", addStatsCallback)

### Export Stats button
lbl = tk.Label(statsLabelframe, text='Export Stats to CSV:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=385, y=0)
writeStatsButton = tk.Button(statsLabelframe, text='Export', font=('Arial Bold', 16), width = 10, command=writeStatsToCsv) 
writeStatsButton.place(x=385, y=30)

### Add Stats button
addStatsButton = tk.Button(statsLabelframe, text='Plot Stats', font=('Arial Bold', 16), width = 10, command=addStatsCallback) 
addStatsButton.place(x=385, y=95)


###############################################################################
#
# TAB 3: PLOT DATA - ATL08 DATA
#
###############################################################################

### Plot ATL08 Data Panel label
dataLayersLabelframe = tk.LabelFrame(tab3, width=545, height=270, text='Add Layers to Plot', font=('Arial Bold', 14))
dataLayersLabelframe.place(x=580, y=10)

# ATL08 Plot text
lbl = tk.Label(dataLayersLabelframe, text='ATL08 Data to Add:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=15)

# ATL08 Y Axis Combo Box
#lbl = tk.Label(atl08DataLabelframe, text='Y Axis:', font=('Arial', 12), anchor = 'w', justify='left')
#lbl.place(x=10, y=100)
yValsBox_atl08 = Combobox(dataLayersLabelframe, width=30)
yValsBox_atl08.place(x=10, y=45)

# Itialize ATL08 plot lists
plotList_atl08 = ('maxCanopy', 'teBestFit', 'teMedian', 
                  'maxCanopyMsl', 'teBestFitMsl', 'teMedianMsl')

# Plot ATL08 Button Callback
def plotAtl08():
        
    # Try plot code
    try:
        
        if(len(atl08Data)>0):
            
            # Get
            gtNumToPlot = gtNumPlotBox.current()
            
            # Get x,y combo box number selections
            xVarNum = xValsBox.current()
            yVarNum = yValsBox_atl08.current()
            
            # Get x,y combo bxx text selections
            xData = eval('atl08Data[' + str(gtNumToPlot) + '].' + plotList[xVarNum])
            yData = eval('atl08Data[' + str(gtNumToPlot) + '].' + plotList_atl08[yVarNum])
            fileName = eval('atl03Data[' + str(gtNumToPlot) + '].atl03FileName')
            gtNum = eval('atl03Data[' + str(gtNumToPlot) + '].gtNum')
            
            # Remove bad data from ATL08
            indsToKeep = yData<=1e20
            xData = xData[indsToKeep]
            yData = yData[indsToKeep]
            
            # Get labels
            xLabel = xlabel_textBox.get()
            yLabel = ylabel_textBox.get()
            title = fileName + ' (' + gtNum + ')'
            
            yName = plotList_atl08[yVarNum]
            
            # Call getPlot function
            getPlot_atl08(xData, yData, xLabel, yLabel, title, yName)
        
        else:
            
            messagebox.showinfo('Error','No ATL08 data to plot.')
            
        # endIf
    
    except:
        
        messagebox.showinfo('Error','Cannot plot data. Please check inputs.')
        
    # endTry
# endDef
    
# Plot ATL08 Button
btn = tk.Button(dataLayersLabelframe, text='Add ATL08', font=('Arial Bold', 16), width = 15, command=plotAtl08) 
btn.place(x=320, y=25)


###############################################################################
#
# TAB 3: PLOT DATA - REFERENCE DATA
#
###############################################################################

#### Plot Reference Data Panel label
#truthPlotLabelframe = tk.LabelFrame(tab3, width=545, height=130, text='Add Reference Data', font=('Arial Bold', 14))
#truthPlotLabelframe.place(x=580, y=10)

# Reference Data text
truthText = ['Add reference data to plot, note:\n' \
             'Reference Time and Delta Time are\n' \
             'calculated using linear interpolation.']
    
# Add text        
lbl = tk.Label(dataLayersLabelframe, text=truthText[0], font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=90)

# Itialize Reference data plot lists
plotList_truth = ('time', 'deltaTime', \
                  'lat', 'lon', \
                  'easting', 'northing', \
                  'crossTrack', 'alongTrack', \
                  'z', 'z', 'class', 'confidence')

plotList_truthNames = ('Reference Time', 'Reference Delta Time', \
                       'Reference Lat', 'Reference Lon', \
                       'Reference Easting', 'Reference Northing', \
                       'Reference Cross-Track', 'Reference Along-Track', \
                       'Reference Height', 'Reference Height', \
                       'Reference Class', 'Reference Conf') 

# Plot Reference Button Callback
def plotTruth():
    
    # Try plot code
    try:
                
        if(len(atlTruthDataFiltered)>0):
        
            # Get
            gtNumToPlot = gtNumPlotBox.current()
            
            # Get x,y combo box number selections
            xVarNum = xValsBox.current()
            yVarNum = yValsBox.current()
            
            # x,y param names
            xParam = plotList_truth[xVarNum]
            yParam = plotList_truth[yVarNum]
            
            if('class' in xParam or 'class' in yParam):
                messagebox.showinfo('Error','Reference data cannot plot classification on x/y axis.')
            elif('confidence' in xParam or 'confidence' in yParam):
                messagebox.showinfo('Error','Reference data cannot plot signal confidence on x/y axis.')
            else:    
                
                # Get x,y combo bxx text selections
                xData = eval('atlTruthDataFiltered[' + str(gtNumToPlot) + '].' + xParam)
                yData = eval('atlTruthDataFiltered[' + str(gtNumToPlot) + '].' + yParam)
                fileName = eval('atl03Data[' + str(gtNumToPlot) + '].atl03FileName')
                gtNum = eval('atl03Data[' + str(gtNumToPlot) + '].gtNum')
            
                # Get labels
                xLabel = xlabel_textBox.get()
                yLabel = ylabel_textBox.get()
                title = fileName + ' (' + gtNum + ')'
                
                yName = plotList_truthNames[yVarNum]
                
                # Call getPlot function
                getPlot_truth(xData, yData, xLabel, yLabel, title, yName)
                
            # endIf
        
        else:
            
            messagebox.showinfo('Error','No Reference data to plot.')
            
        # endIf
    
    except:
        
        messagebox.showinfo('Error','Cannot plot data. Please check inputs.')
        
    # endTry
# endDef

# Plot Truth Button
btn = tk.Button(dataLayersLabelframe, text='Add Reference', font=('Arial Bold', 16), width = 15, command=plotTruth) 
btn.place(x=320, y=105)


###############################################################################
#
# TAB 3: PLOT DATA - CORRECTED MEASURED
#
###############################################################################

#### Plot Corrected Measured Panel label
#measCorrPlotLabelframe = tk.LabelFrame(tab3, width=545, height=130, text='Add Shifted ICESat-2 Data', font=('Arial Bold', 14))
#measCorrPlotLabelframe.place(x=580, y=150)

# Corrected Measured text
measCorrText = ['Plot shifted ICESat-2 track that has\n'\
                'been corrected in XYZ relative to the\n' \
                'reference data.']

# Corrected Measured Data Y Axis Combo Box
lbl = tk.Label(dataLayersLabelframe, text=measCorrText[0], font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=170)

plotList_measCorr = ('time', 'deltaTime', \
                     'latArray', 'lonArray', \
                     'eastingArray', 'northingArray', \
                     'crossTrackArray', 'alongTrackArray', \
                     'zArray', 'zMslArray', 'classification', 'signalConf')

plotList_measCorrNames = ('Shifted ATL03 Time', 'Shifted ATL03 Delta Time',\
                          'Shifted ATL03 Lat', 'Shifted ATL03 Lon', \
                          'Shifted ATL03 Easting', 'Shifted ATL03 Northing', \
                          'Shifted ATL03 Cross-Track', 'Shifted ATL03 Along-Track', \
                          'Shifted ATL03 Height HAE', 'Shifted ATL03 Height MSL', \
                          'Shifted ATL03 Classification', 'Shifted ATL03 Signal Confidence') 

# Plot Corrected Measured Button Callback
def plotMeasCorr():
    
    # Try plot code
    try:
                
        if(len(atlCorrections)>0):
            
            # Get
            gtNumToPlot = gtNumPlotBox.current()
            
            # Get x,y combo box number selections
            xVarNum = xValsBox.current()
            yVarNum = yValsBox.current()
            
            # x,y param names
            xParam = plotList_measCorr[xVarNum]
            yParam = plotList_measCorr[yVarNum]
                
            # Get x,y combo bxx text selections
            xData = eval('atlCorrections[' + str(gtNumToPlot) + '].' + xParam)
            yData = eval('atlCorrections[' + str(gtNumToPlot) + '].' + yParam)
            fileName = eval('atl03Data[' + str(gtNumToPlot) + '].atl03FileName')
            gtNum = eval('atl03Data[' + str(gtNumToPlot) + '].gtNum')
            
            # Get labels
            xLabel = xlabel_textBox.get()
            yLabel = ylabel_textBox.get()
            title = fileName + ' (' + gtNum + ')'            
            
            # Get Filter data type (classification or signal confidence) and filter numbers
            filterChoice = filterBox.get()
            if('  ' in filterChoice.lower()):
                filterType = []
                filterData = []
                filterNum = []
            elif('class' in filterChoice.lower()):
                filterType = filterChoice
                filterData = eval('atlCorrections[' + str(gtNumToPlot) + '].classification')
                filterTF = [filter0ChkState.get(), filter1ChkState.get(), filter2ChkState.get(), filter3ChkState.get(), filter4ChkState.get()]    
                filterNum = filterState[filterTF]
            elif('signal' in filterChoice.lower()):
                filterType = filterChoice
                filterData = eval('atlCorrections[' + str(gtNumToPlot) + '].signalConf')
                filterTF = [filter0ChkState.get(), filter1ChkState.get(), filter2ChkState.get(), filter3ChkState.get(), filter4ChkState.get()]    
                filterNum = filterState[filterTF]
            # endIf
            
            # Get output path
            outPath = outPath_textBox.get().strip()
            
            # Get original plot title
            origTitle = title
            
            # Get gtNum
            gtNum = gtNumPlotBox.current()
            
            # Call getPlot_measCorr function
            getPlot_measCorr(xData, yData, xLabel, yLabel, title, outPath, 
                             origTitle, filterType, filterData, filterNum)
    
        else:
        
            messagebox.showinfo('Error','No Shifted ATL03 data to plot.')
            
        # endIf
    
    except:
        
        messagebox.showinfo('Error','Cannot plot data. Please check inputs.')
        
    # endTry
# endDef
    
# Plot Corrected Measured Button
btn = tk.Button(dataLayersLabelframe, text='Add Shifted Data', font=('Arial Bold', 16), width = 15, command=plotMeasCorr) 
btn.place(x=320, y=185)


###############################################################################
#
# TAB 3: PLOT DATA - PHO SHOW
#
###############################################################################

### Plot Corrected Measured Panel label
phoShowPlotLabelframe = tk.LabelFrame(tab3, width=545, height=85, text='PhoSHOW', font=('Arial Bold', 14))
phoShowPlotLabelframe.place(x=580, y=290)

# Truth Data Y Axis Combo Box
phoShowText = ['View output data in HTML format:']

# Corrected Measured Data Y Axis Combo Box
lbl = tk.Label(phoShowPlotLabelframe, text=phoShowText[0], font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=12)

# PhoShow Button Callback
def phoShow():
    
    # Try plot code
    try:
        
        # Get GT Num to plot
        gtNumToPlot = gtNumPlotBox.current()
        
        # Get x,y combo bxx text selections
        ytrack = eval('atl03Data[' + str(gtNumToPlot) + '].alongTrack')
        h_ph = eval('atl03Data[' + str(gtNumToPlot) + '].z')
        classification = eval('atl03Data[' + str(gtNumToPlot) + '].classification')
        direction = eval('atl03Data[' + str(gtNumToPlot) + '].trackDirection')
        lat = eval('atl03Data[' + str(gtNumToPlot) + '].lat')
        lon = eval('atl03Data[' + str(gtNumToPlot) + '].lon')
        
        # Input/output paths
        inFilePath = cwd
        outFilePath = os.path.normpath(outPath_textBox.get().strip())
                
        # Update message
        print('Writing HTML File...')
        
        # Call PhoShow code
        html_file = createHTMLChart(ytrack, h_ph, classification, lat, lon, 
                    direction,
                    online = None,
                    classification_list = [1,2,3],
                    input_folder = inFilePath,
                    output_folder = outFilePath,
                    in_file03_name = 'ATL03', 
                    blank = viewerBlank_html, 
                    blanki = viewerBlankOnline_html)
        
        # Update message
        print('HTML File Complete!')
        
        # Open PhoShow in HTML browser
        webbrowser.open_new_tab(html_file)
    
    except:
        
        # Print error message
        print('Cannot open PhoSHOW. Please check inputs')
        
    # endTry
# endDef
    
# PhoShow Button
btn = tk.Button(phoShowPlotLabelframe, text='PhoSHOW', font=('Arial Bold', 16), width = 15, command=phoShow) 
btn.place(x=320, y=0)


###############################################################################
#
# TAB 3: PLOT DATA - LOAD FIGURE
#
###############################################################################

### Panel label
loadPlotLabelframe = tk.LabelFrame(tab3, width=545, height=80, text='Load Figure', font=('Arial Bold', 14))
loadPlotLabelframe.place(x=580, y=385)

# Plot text
lbl = tk.Label(loadPlotLabelframe, text='Load a saved figure (.pkl format):', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=12)

def loadPlot():
    
    # Check output file path
    startDir = outPath_textBox.get().strip()
    if('' == startDir or '.' == startDir):
        startDir = cwd
    else:
        outPathExists = os.path.isdir(startDir)
        if(not outPathExists):
            os.mkdir(outFilePath)
            startDir = cwd
        # endIf
    # endIf
    
    # Ask user for input .pkl files
    pklFiles = filedialog.askopenfilename(initialdir = startDir, title = 'Select File(s) to Plot', filetypes = [('*fig*.pkl files','*fig*.pkl')], multiple=True)

    # If the user selection is not empty, then continue
    if(pklFiles != '.'):
            
        # Loop through all .pkl files selected
        for pklFile in pklFiles:
            
            # Plot .pkl file
            getPklPlot(os.path.normpath(pklFile))        
    
        # endFor
    # endIf
# endDef
    
# PhoShow Button
btn = tk.Button(loadPlotLabelframe, text='Browse', font=('Arial Bold', 16), width = 15, command=loadPlot) 
btn.place(x=320, y=0)



###############################################################################
#
# TAB 4: HELP PAGE
#
###############################################################################

### Plot Corrected Measured Panel label
helpLabelframe1 = tk.LabelFrame(tab4, width=545, height=450, text='PhoREAL Help Information', font=('Arial Bold', 14))
helpLabelframe1.place(x=15, y=10)

phoreal_help_info1 = \
'Get ICESat-2 Data Input Section:\n' \
'--------------------------------------------\n' \
'This section handles the input ATL03/ATL08 files from ICESat-2\n' \
'\n' \
'+ ATL03 File: Path to input ATL03 .h5 file(s) or directory\n' \
'+ ATL08 File: Path to input ATL08 .h5 file(s) or directory\n' \
'+ Output Directory: Path to directory for output data\n' \
'+ Ground Track Numbers: Option to select ICESat-2 ground tracks\n' \
'+ Trim ICESat-2 Data Options: Methods to trim ICESat-2 data\n' \
'     - None: No trimming\n' \
'     - Manual: Trim to user-specified latitude or time min/max bounds\n' \
'     - Auto: Trim to reference region bounds (uses kmlBounds.txt file)\n' \
'+ Create Output Files: Option to create output files for ICESat-2 data\n' \
'\n' \
'Get Reference Data Input Section:\n' \
'-----------------------------------------------\n' \
'This section handles the reference data used to find ICESat-2 offsets\n' \
'\n' \
'+ Use Existing Data: Option to use existing reference buffer data\n' \
'+ Buffer Size: Buffer to create around ICESat-2 track in reference data\n' \
'+ Save Reference File: Option to create output reference buffer .las file\n' \
'+ Reference File(s): Path to reference file(s) or directory (.las or .tif)' \


lbl = tk.Label(helpLabelframe1, text=phoreal_help_info1, font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=10)

### Plot Corrected Measured Panel label
helpLabelframe2 = tk.LabelFrame(tab4, width=545, height=440, text='', font=('Arial Bold', 14))
helpLabelframe2.place(x=580, y=20)

phoreal_help_info2 = \
'Find ICESat-2 Offsets Relative to Reference Data Section:\n' \
'-------------------------------------------------------------------------------\n' \
'This section slides the ICESat-2 data over the reference data and finds the\n' \
'offset with the minimum Mean Absolute Error in the Z direction (location\n' \
'of best fit relative to reference data)\n' \
'\n' \
'+ Cross-Track Bounds (m): Cross-track search area [min, max] or one value\n' \
'+ Along-Track Bounds (m): Along-track search area [min, max] or one value\n' \
'+ Use Fixed Vertical Shift: Option to use a fixed vertical shift value\n' \
'+ Vertical Shift (m): Vertical shift value if previous option is selected\n' \
'+ Grid Resolution(s) (m): Raster resolution(s) to grid reference data\n' \
'+ Use ICESat-2 Heights: Option to use Ellipsoidal (HAE) heights or \n' \
'   Orthometric (MSL) heights when comparing to reference data\n' \
'   Note: Geoidal heights have been calculated on the ATL03\n' \
'         data product using the EGM2008 ellipsoid\n' \
'+ Use ICESat-2 Signal Confidence Value(s): Option to use ICESat-2\n' \
'   signal confidence values to filter measured data\n' \
'     - Input signal confidence values to filter ICESat-2 data\n' \
'+ Use Reference Ground Index: Option to use reference ground index value\n' \
'   to filter reference data (Requires ATL08 file)\n' \
'     - Reference ground index value (ASPRS ground class = 2)\n' \
'+ Save Shifted ICESat-2 File: Option to save ICESat-2 shifted XYZ file\n' \
'+ Make Output Plots: Option to create output plots'

lbl = tk.Label(helpLabelframe2, text=phoreal_help_info2, font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=10)


###############################################################################
#
# TAB 5: ABOUT PAGE
#
###############################################################################

### Panel label
aboutLabelFrame = tk.LabelFrame(tab5, width=1110, height=450, text='', font=('Arial Bold', 14))
aboutLabelFrame.place(x=15, y=10)

aboutInfo = \
'PhoREAL Description\n' \
'-------------------------------------------------------------------------------------------------------------\n' \
'PhoREAL is a comprehensive tool for analyzing and plotting ATL03 and ATL08\n' \
'data from ICESat-2. All of the capabilities of PhoREAL %s\n' \
'have been compiled into this Windows Graphical User Interface (GUI).\n' \
'\n' \
'PhoREAL Licensing Information\n' \
'------------------------------------------------------------------------------------------\n' \
'This package is free software; the copyright holder gives unlimited\n' \
'permission to copy and/or distribute, with or without modification, as\n' \
'long as this notice is preserved.\n' \
'\n' \
'Authors\n' \
'--------------------\n' \
'Mike Alonzo\n' \
'Eric Guenther\n' \
'\n' \
'Contact Information\n' \
'------------------------------------------------------------------------------\n' \
'Any feedback can be sent to phoreal@arlut.utexas.edu\n' \
'\n\n' \
'Copyright 2019 Applied Research Laboratories, The University of Texas at Austin' %phoRealVersion

lbl = tk.Label(aboutLabelFrame, text=aboutInfo, font=('Arial', 12), anchor = 'w', justify='center')
lbl.place(x=260, y=10)

# Open GUI window
window.mainloop()