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
from tkinter import ttk
from tkinter import messagebox
import numpy as np
import os

from getAtlMeasuredSwath_auto import getAtlMeasuredSwath
from getAtlTruthSwath_auto import getAtlTruthSwath
from getMeasurementError_auto import getMeasurementError, offsetsStruct

from icesatPlot import (getPlot, getPlot_atl08, getPlot_truth, getPlot_measCorr)
from icesatIO import createHTMLChart
from icesatUtils import (superFilter, getNameParts)
from gui_logo import images
from gui_addins import (viewerBlank_html, viewerBlankOnline_html)

import webbrowser


# Print opening message
print('\n')
print('*************************************************')
print('PhoREAL GUI is running...')
print('*************************************************')
print('\n')


window = tk.Tk()

# GUI title
window.title('PhoREAL v3.0 - Applied Research Labs (The University of Texas at Austin)')

# GUI size
window.geometry('1140x500')
window.resizable(width=False, height=False)

# Set icon image for GUI
if os.name == 'nt':
    window.wm_iconbitmap(images)


# Create top level tabs
tabControl = ttk.Notebook(window)  
       
# Create tab1
tab1 = ttk.Frame(tabControl)            
tabControl.add(tab1, text='Home')      
tabControl.pack(expand=1, fill='both') 

# Create tab2
tab2 = ttk.Frame(tabControl)            
tabControl.add(tab2, text='Plot Data')      
tabControl.pack(expand=1, fill='both') 

# Create tab3
tab3 = ttk.Frame(tabControl)            
tabControl.add(tab3, text='Help')      
tabControl.pack(expand=1, fill='both')

# Create tab4
tab4 = ttk.Frame(tabControl)            
tabControl.add(tab4, text='About')      
tabControl.pack(expand=1, fill='both')

# Set tab font sytle
s = ttk.Style()
s.configure('TNotebook.Tab', font=('Arial','10','bold'))

# Get current working directory
cwd = os.path.normpath(os.getcwd())


###############################################################################
#
# TAB 1: GET MEASURED DATA
#
###############################################################################

# Define global variable
global atl03FileExists, atl08FileExists, truthDataExists
atl03FileExists = False
atl08FileExists = False
outPathExists = True
truthDataExists = False

### Panel label
measLabelframe = tk.LabelFrame(tab1, width=545, height=450, text='Get ICESat-2 Data Input', font=('Arial Bold', 14))
measLabelframe.place(x=15, y=10)

### ATL03 INPUT FILE

# Check ATL08 Button Callback
def checkATL03():
    
    global atl03FileExists
    
    atl03File = os.path.normpath(atl03_textBox.get().strip())
    atl03FileExists = os.path.exists(atl03File)
    
# endDef
    
# ATL03 File Entry Box
lbl = tk.Label(measLabelframe, text='ATL03 File:', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=10)
atl03_textBox = tk.Entry(measLabelframe, width=70)
atl03_textBox.place(x=10, y=40)
atl03_textBox.bind('<Return>', (lambda event: checkATL03()))

# Browse ATL03 File Button Callback
def browseAtl03():
    
    global atl03FileExists
    
    currentData = atl03_textBox.get()
    atl03File = os.path.normpath(filedialog.askopenfilename(title = 'Select ATL03 .h5 file', filetypes = [('h5 files','ATL03_*.h5')]))
    if(atl03File != '.'):
        atl03FileExists = os.path.exists(atl03File)
        if(atl03FileExists):
            atl03_textBox.delete(0,len(currentData))
            atl03_textBox.insert(0,atl03File) 
        # endIf
    # endIf
# endDef
    
# ATL03 File Browse Button
atl03BrowseButton = tk.Button(measLabelframe, text='Browse', font=('Arial Bold', 12), command=browseAtl03) 
atl03BrowseButton.place(x=450, y=30)


### ATL08 INPUT FILE

# Check ATL08 Button Callback
def checkATL08():
    
    global atl08FileExists
    
    atl08File = os.path.normpath(atl08_textBox.get().strip())
    atl08FileExists = os.path.exists(atl08File)
    useMeasErrorSection = useMeasErrorSectionChkState.get()
    if(atl08FileExists and useMeasErrorSection):
        useGroundIndex_checkBox.config(state = 'normal')
    # endIf
# endDef
    
# ATL08 File Entry Box
lbl = tk.Label(measLabelframe, text='ATL08 File (Optional):', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=70)
atl08_textBox = tk.Entry(measLabelframe, width=70) 
atl08_textBox.place(x=10, y=100)
atl08_textBox.bind('<Return>', (lambda event: checkATL08()))

# Browse ATL08 File Button Callback
def browseAtl08():
    
    global atl08FileExists
    
    try:
        if(atl03FileExists):
            
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
    
    currentData = atl08_textBox.get()
    atl08File = os.path.normpath(filedialog.askopenfilename(title = 'Select ATL08 .h5 file', filetypes = [('h5 files', filterText)]))
    if(atl08File != '.'):
        atl08FileExists = os.path.exists(atl08File)
        useMeasErrorSection = useMeasErrorSectionChkState.get()
        if(atl08FileExists):
            atl08_textBox.delete(0,len(currentData))
            atl08_textBox.insert(0,atl08File)
            if(useMeasErrorSection):
                useGroundIndex_checkBox.config(state = 'normal')
            # endIf
        else:
            useGroundIndex_checkBox.config(state = 'disabled')
        # endIf
    # endIf
# endDef
    
# ATL08 File Browse Button
atl08BrowseButton = tk.Button(measLabelframe, text='Browse', font=('Arial Bold', 12), command=browseAtl08) 
atl08BrowseButton.place(x=450, y=90)


### Output Path

# Output Path Entry Box
lbl = tk.Label(measLabelframe, text='Output Directory:', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=130)
outPath_textBox = tk.Entry(measLabelframe, width=70) 
outPath_textBox.place(x=10, y=160)
outPath_textBox.insert(0,cwd)

# Browse Output File Button Callback
def browseOutput():
    
    global outPathExists
    
    currentData = outPath_textBox.get()
    outPath = os.path.normpath(filedialog.askdirectory(title = 'Select Output File Path'))
    if(outPath != '.'):
        outPath_textBox.delete(0,len(currentData))
        outPath_textBox.insert(0,outPath) 
    # endIf
    
# endDef
    
# Output Directory Browse Button
outPathBrowseButton = tk.Button(measLabelframe, text='Browse', font=('Arial Bold', 12), command=browseOutput) 
outPathBrowseButton.place(x=450, y=150)


### GT Number Text Label
lbl = tk.Label(measLabelframe, text='Ground Track Numbers:', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=190)

gtNumsAll = np.array(['GT1R','GT1L','GT2R','GT2L','GT3R','GT3L'])

### GT Number Checkboxes
gtNum1rChkState = tk.BooleanVar()
gtNum1rChkState.set(True)
gtNum1r_checkBox = tk.Checkbutton(measLabelframe, text = 'GT1R', font=('Arial', 12), var = gtNum1rChkState) 
gtNum1r_checkBox.place(x=10, y=220)

gtNum2rChkState = tk.BooleanVar()
gtNum2rChkState.set(True)
gtNum2r_checkBox = tk.Checkbutton(measLabelframe, text = 'GT2R', font=('Arial', 12), var = gtNum2rChkState) 
gtNum2r_checkBox.place(x=90, y=220)

gtNum3rChkState = tk.BooleanVar()
gtNum3rChkState.set(True)
gtNum3r_checkBox = tk.Checkbutton(measLabelframe, text = 'GT3R', font=('Arial', 12), var = gtNum3rChkState) 
gtNum3r_checkBox.place(x=170, y=220)

gtNum1lChkState = tk.BooleanVar()
gtNum1lChkState.set(False)
gtNum1l_checkBox = tk.Checkbutton(measLabelframe, text = 'GT1L', font=('Arial', 12), var = gtNum1lChkState) 
gtNum1l_checkBox.place(x=250, y=220)

gtNum2lChkState = tk.BooleanVar()
gtNum2lChkState.set(False)
gtNum2l_checkBox = tk.Checkbutton(measLabelframe, text = 'GT2L', font=('Arial', 12), var = gtNum2lChkState) 
gtNum2l_checkBox.place(x=330, y=220)

gtNum3lChkState = tk.BooleanVar()
gtNum3lChkState.set(False)
gtNum3l_checkBox = tk.Checkbutton(measLabelframe, text = 'GT3L', font=('Arial', 12), var = gtNum3lChkState) 
gtNum3l_checkBox.place(x=410, y=220)


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
createCsvChkState.set(True)
createCsv_checkBox = tk.Checkbutton(measLabelframe, text = 'ATL03 .csv File', font=('Arial', 12), var = createCsvChkState) 
createCsv_checkBox.place(x=190, y=390)

### Create ATL08 .kml Text Entry
createATL08KmlChkState = tk.BooleanVar()
createATL08KmlChkState.set(False)
createATL08Kml_checkBox = tk.Checkbutton(measLabelframe, text = 'ATL08 .kml File', font=('Arial', 12), var = createATL08KmlChkState) 
createATL08Kml_checkBox.place(x=370, y=360)

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
    
### Activate Truth section
useTruthSectionChkState = tk.BooleanVar()
useTruthSectionChkState.set(True)
useTruthSection_checkBox = tk.Checkbutton(tab1, text = '', font=('Arial', 12), var = useTruthSectionChkState, command = useTruthSectionCallback) 
useTruthSection_checkBox.place(x=1110, y=10)

# Browse Truth Data Directory File Button Callback
def browseTruthDir():
    
    global truthDataExists
    
    currentData = truthDataDir_textBox.get()
    useTruth = useExistingTruthChkState.get()
    
    if(useTruth):
        truthData = os.path.normpath(filedialog.askopenfilename(title = 'Select Reference File', filetypes = [('reference files','*.las *.laz *.tif')]))
        truthDataExists = os.path.exists(truthData)
    else:
        truthData = os.path.normpath(filedialog.askdirectory(title = 'Select Reference Directory'))
        truthDataExists = os.path.isdir(truthData)
    # endIf
    
    if(truthDataExists and truthData!='.'):
        truthDataDir_textBox.delete(0,len(currentData))
        truthDataDir_textBox.insert(0,truthData) 
    else:
        truthDataExists = False
    # endIf
# endDef
    
# Use Existing Truth button callback
def checkUseExistingTruth():
    # useExistingTruthChkState.set(True) # CHANGE AFTER DEMO !!!
    useTruth = useExistingTruthChkState.get()
    if(useTruth):
        truthBuffer_textBox.config(state = 'disabled')
        truthDir_lbl.config(text = 'Reference File Name:')
    else:
        truthBuffer_textBox.config(state = 'normal')
        truthDir_lbl.config(text = 'Reference Directory:')
    # endIf
# endDef

### Use Existing Truth Check Box
useExistingTruthChkState = tk.BooleanVar()
useExistingTruthChkState.set(True)
useExistingTruth_checkBox = tk.Checkbutton(truthLabelframe, text = 'Use Existing Data', font=('Arial', 12), var = useExistingTruthChkState, command = checkUseExistingTruth) 
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
createTruthFileChkState.set(False)
createTruthFile_checkBox = tk.Checkbutton(truthLabelframe, text = 'Create Reference File', font=('Arial', 12), var = createTruthFileChkState) 
createTruthFile_checkBox.place(x=350, y=10)

# Browse Truth Data Directory File Button Callback
def checkTruthDir():
    
    global truthDataExists
    
    truthData = os.path.normpath(truthDataDir_textBox.get().strip())
    useTruth = useExistingTruthChkState.get()
    if(useTruth):
        truthDataExists = os.path.exists(truthData)
    else:
        truthDataExists = os.path.isdir(truthData)
    # endIf
    
    # endIf
# endDef
    
### Truth Data Directory Entry Box
truthDir_lbl = tk.Label(truthLabelframe, text='Reference File Name:', font=('Arial', 12), anchor = 'w', justify='left')
truthDir_lbl.place(x=10, y=50)
truthDataDir_textBox = tk.Entry(truthLabelframe, width=43)
truthDataDir_textBox.place(x=170, y=55)
truthDataDir_textBox.bind('<Return>', (lambda event: checkTruthDir()))

### Truth Data Dir Browse Button
truthDataDirBrowseButton = tk.Button(truthLabelframe, text='Browse', font=('Arial Bold', 12), command=browseTruthDir) 
truthDataDirBrowseButton.place(x=450, y=50)


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
    useTruthSection = useTruthSectionChkState.get()
    useMeasErrorSection = useMeasErrorSectionChkState.get()
    if(useTruthSection and useMeasErrorSection):
        for child in measErrorLabelframe.winfo_children():
            child.configure(state='normal')
        # endFor
        useVertShiftCallback()
        useMeasSigConfCallback()
        useGroundIndexConfCallback()
        if(not atl08FileExists):
            useGroundIndex_checkBox.config(state = 'disabled')
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
useMeasErrorSectionChkState.set(True)
useMeasErrorSection_checkBox = tk.Checkbutton(tab1, text = '', font=('Arial', 12), var = useMeasErrorSectionChkState, command = useMeasErrorSectionCallback) 
useMeasErrorSection_checkBox.place(x=1110, y=140)

### Cross-Track Bounds Entry Box
lbl = tk.Label(measErrorLabelframe, text='Cross-Track Bounds (m):', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=10)
crossTrackBounds_textBox = tk.Entry(measErrorLabelframe, width=10)
crossTrackBounds_textBox.place(x=195, y=15)
crossTrackBounds_textBox.insert(0,'-48, 48')

### Along-Track Bounds Entry Box
lbl = tk.Label(measErrorLabelframe, text='Along-Track Bounds (m):', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=270, y=10)
alongTrackBounds_textBox = tk.Entry(measErrorLabelframe, width=10)
alongTrackBounds_textBox.place(x=455, y=15)
alongTrackBounds_textBox.insert(0,'-48, 48')

### Multi-Resolutional Stepdown Entry Box
lbl = tk.Label(measErrorLabelframe, text='Grid Resolution(s) (m):', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=80)
multiresStepdown_textBox = tk.Entry(measErrorLabelframe, width=15)
multiresStepdown_textBox.place(x=175, y=85)
multiresStepdown_textBox.insert(0,'8, 4, 2, 1')

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
useGroundIndex_checkBox.config(state = 'disabled')

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
createMeasCorrFile_checkBox = tk.Checkbutton(measErrorLabelframe, text = 'Create Shifted ICESat-2 File', font=('Arial', 12), var=createMeasCorrFileChkState) 
createMeasCorrFile_checkBox.place(x=10, y=185)

### Make Plots Check Box
makePlotsChkState = tk.BooleanVar()
makePlotsChkState.set(True)
makePlots_checkBox = tk.Checkbutton(measErrorLabelframe, text = 'Make Output Plots', font=('Arial', 12), var=makePlotsChkState) 
makePlots_checkBox.place(x=280, y=185)


# Run Button Callback
def runAtl03():
    
    # Update status bar
    statusBar['value'] = 0
    window.update()
        
    # Make atlMeasuredData a global variable
    global atl03Data, atl08Data, atlTruthData, atlTruthDataFiltered, atlCorrections, atl03FileExists, atl08FileExists, truthDataExists
        
    # Initialize variables
    atl03Data = []
    atl03DataSingle = []
    atl08Data = []
    atl08DataSingle = []
    atlTruthData = []
    atlTruthDataSingle = []
    atlTruthDataFiltered = []
    atlTruthDataFilteredSingle = []
    atlCorrections = []
    atlCorrectionsSingle = []
    
    # Check output file path
    outFilePath = outPath_textBox.get().strip()
    if('' == outFilePath):
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
          
    if(atl03FileExists and outPathExists):
    
        # Disable run button
        RunButton.config(state=tk.DISABLED)
        
        # Try code
        try:
            
            # Get Measured Data inputs
            atl03FilePath = atl03_textBox.get().strip()
            if(atl08FileExists):
                atl08FilePath = atl08_textBox.get().strip()
            else:
                atl08FilePath = ''
            # endIf
            gtNumsTF = [gtNum1rChkState.get(), gtNum1lChkState.get(), gtNum2rChkState.get(), gtNum2lChkState.get(), gtNum3rChkState.get(), gtNum3lChkState.get()]    
            gtNums = gtNumsAll[gtNumsTF]  
        
            # Get Truth Data inputs
            if(useTruthSection):
                useExistingTruth = useExistingTruthChkState.get()
                truthSwathDir = truthDataDir_textBox.get().strip()
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
            createATL08KmlFile = createATL08KmlChkState.get()
            createATL08CsvFile = createATL08CsvChkState.get()
                
            # Loop through all gt nums
            headerData = False
            rotationData = False
            for i in range(0,len(gtNums)):
                
                # Update status bar
                statusBar['value'] = 0
                window.update()
            
                # Get GT Nums
                gtNum = gtNums[i].lower()
                    
                # Run getAtlMeasuredSwath function
                print('Getting Measured Data...\n')
                atl03DataSingle, atl08DataSingle, headerData, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, gtNum, trimInfo, createLasFile, createKmlFile, createATL08KmlFile, createCsvFile, createATL08CsvFile)
                         
                # Update status bar
                progressNum = 30
                statusBar['value'] = progressNum
                window.update()
                
                # Append measured data for each run
                atl03Data.append(atl03DataSingle)
                if(atl08FileExists):
                    atl08Data.append(atl08DataSingle)
                # endIf

                if(useTruthSection and truthDataExists):
                    print('Getting Reference Data...\n')
                    atlTruthDataSingle = getAtlTruthSwath(atl03DataSingle, headerData, rotationData, useExistingTruth, truthSwathDir, buffer, outFilePath, createTruthFile)
                
                    # Run superfilter on data
                    print('   Running superfilter on data...')
                    atlTruthDataFilteredSingle, _ = superFilter(atl03DataSingle, atlTruthDataSingle, xBuf = 1, classCode = [])
                    print('\n')
            
                    # Append reference data for each run
                    atlTruthData.append(atlTruthDataSingle)
                    atlTruthDataFiltered.append(atlTruthDataFilteredSingle)
                    
                    # Update status bar
                    progressNum = 60
                    statusBar['value'] = progressNum
                    window.update()
                
                # endIf
                
                if(useMeasErrorSection):
                    print('Finding Measurement Offsets...\n')
                    atlCorrectionsSingle = getMeasurementError(atl03DataSingle, atlTruthDataSingle, rotationData, outFilePath, useMeasSigConf, filterData, offsets, createMeasCorrFile, makePlots, showPlots)
                
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
                            
            # endFor
            
            if(len(gtNums)==0):
                
                # Reset status bar and print warning
                statusBar['value'] = 0
                print('\nWARNING: No Ground Track Selected.\n')
            # endIf
                    
            if(atl03Data):
                
                # Set Ground Tracks to plot
                gtNumsTuple = tuple(gtNums)
                gtNumPlotBox['values'] = gtNumsTuple
                gtNumPlotBox.current(0)
                
                if(atl03Data[0].zone=='3413' or atl03Data[0].zone=='3976'):
                    
                    plotVarsTuple = ('Time (sec)', 'Latitude (deg)', 'Longitude (deg)', \
                                'Polar Stereo X (m)', 'Polar Stereo Y (m)', \
                                'Cross-Track (m)', 'Along-Track (m)', \
                                'Height (m)', \
                                'Classification', 'Signal Confidence')
                    
                else:
                    
                    plotVarsTuple = ('Time (sec)', 'Latitude (deg)', 'Longitude (deg)', \
                                'UTM Easting (m)', 'UTM Northing (m)', \
                                'Cross-Track (m)', 'Along-Track (m)', \
                                'Height (m)', \
                                'Classification', 'Signal Confidence')
                    
                # endIf
                    
                # Set X Vals to plot
                xValsBox['values'] = plotVarsTuple
                xValsBox.current(0)
                
                # Set Y Vals to plot
                yValsBox['values'] = plotVarsTuple
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
                
#                # Set Title label
#                gtNumVal = gtNumPlotBox.get()
#                currentData = title_textBox.get()
#                title_textBox.delete(0,len(currentData))
#                title_textBox.insert(0,gtNumVal) 
                
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
            
            if(atl08Data):
                
                # Set Y Vals to plot
                yValsTuple_atl08 = ('Max Canopy (m)', 'Terrain Best Fit (m)', 'Terrain Median (m)')
                yValsBox_atl08['values'] = yValsTuple_atl08
                yValsBox_atl08.current(0)
                
            # endIf
            
            # Print completion message
            print('RUN COMPLETE.')
        
        except:
            
            # Enable run button
            RunButton.config(state=tk.NORMAL)
        
            print('Could not process data. Please check inputs.')
        #endTry
        
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
# TAB 2: PLOT OPTIONS
#
###############################################################################

### Panel label
labelframe = tk.LabelFrame(tab2, width=545, height=380, text='Plot ICESat-2 Data Options', font=('Arial Bold', 14))
labelframe.place(x=15, y=10)

# Plot text
lbl = tk.Label(tab2, text='ATL03 Plotting Options:', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=50)

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
    
## GT Num Plot Callback
#def gtNumPlotCallback(event):
#    gtNumVal = gtNumPlotBox.get()
#    currentData = title_textBox.get()
#    title_textBox.delete(0,len(currentData))
#    title_textBox.insert(0,gtNumVal)   
## endDef
  
    
# Filter Number Checkboxes
filter0ChkState = tk.BooleanVar()
filter0ChkState.set(False)
filter0_checkBox = tk.Checkbutton(tab2, text = '0', font=('Arial', 12), var = filter0ChkState) 

filter1ChkState = tk.BooleanVar()
filter1ChkState.set(False)
filter1_checkBox = tk.Checkbutton(tab2, text = '1', font=('Arial', 12), var = filter1ChkState) 

filter2ChkState = tk.BooleanVar()
filter2ChkState.set(False)
filter2_checkBox = tk.Checkbutton(tab2, text = '2', font=('Arial', 12), var = filter2ChkState) 

filter3ChkState = tk.BooleanVar()
filter3ChkState.set(False)
filter3_checkBox = tk.Checkbutton(tab2, text = '3', font=('Arial', 12), var = filter3ChkState) 

filter4ChkState = tk.BooleanVar()
filter4ChkState.set(False)
filter4_checkBox = tk.Checkbutton(tab2, text = '4', font=('Arial', 12), var = filter4ChkState) 
     
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
    
# GT Num Plot Combo Box
lbl = tk.Label(tab2, text='Ground Track:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=80)
gtNumPlotBox = Combobox(tab2, width=10)
gtNumPlotBox.place(x= 150, y = 82)
#gtNumPlotBox.bind("<<ComboboxSelected>>", gtNumPlotCallback)

# X Axis Combo Box
lbl = tk.Label(tab2, text='X Axis:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=110)
xValsBox = Combobox(tab2)
xValsBox.place(x= 90, y = 112)
xValsBox.bind("<<ComboboxSelected>>", xAxisCallback)

# X Label Entry Box
lbl = tk.Label(tab2, text='X Label:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=260, y=110)
xlabel_textBox = tk.Entry(tab2, width=30)
xlabel_textBox.place(x=330, y=112)

# Y Axis Combo Box
lbl = tk.Label(tab2, text='Y Axis:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=140)
yValsBox = Combobox(tab2)
yValsBox.place(x= 90, y = 142)
yValsBox.bind("<<ComboboxSelected>>", yAxisCallback)

# Y Label Entry Box
lbl = tk.Label(tab2, text='Y Label:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=260, y=140)
ylabel_textBox = tk.Entry(tab2, width=30)
ylabel_textBox.place(x=330, y=142)

## Title Entry Box
#lbl = tk.Label(tab2, text='Title:', font=('Arial', 12), anchor = 'w', justify='left')
#lbl.place(x=30, y=170)
#title_textBox = tk.Entry(tab2, width=35)
#title_textBox.place(x=80, y=175)

# Filter On Combo Box
lbl = tk.Label(tab2, text='Filter On:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=190)
filterBox = Combobox(tab2, width=27)
filterBox.place(x= 110, y = 192)
filterBox.bind("<<ComboboxSelected>>", filterChoiceCallback)

# Itialize lists
plotList = ('time', 'lat', 'lon', 'easting', 'northing', \
            'crossTrack', 'alongTrack', 'z', 'classification', 'signalConf')
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
        getPlot(xData, yData, xLabel, yLabel, title, outPath, origTitle, atl03Data[gtNum], filterType, filterData, filterNum)
    
    except:
        
        print('Cannot plot data. Please check inputs')
        
    # endTry
# endDef
    
# Plot ATL03 Button
btn = tk.Button(tab2, text='PLOT', font=('Arial Bold', 16), width = 15, command=plotAtl03) 
btn.place(x=335, y=178)  

# ATL08 Plot text
lbl = tk.Label(tab2, text='Add ATL08 Data:', font=('Arial Bold', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=305)

# ATL08 Y Axis Combo Box
lbl = tk.Label(tab2, text='Y Axis:', font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=30, y=335)
yValsBox_atl08 = Combobox(tab2, width=30)
yValsBox_atl08.place(x=90, y=337)

# Itialize ATL08 plot lists
plotList_atl08 = ('maxCanopy', 'teBestFit', 'teMedian')

# Plot ATL08 Button Callback
def plotAtl08():
    
    # Try plot code
    try:
        
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
    
    except:
        
        print('Cannot plot data. Please check inputs')
        
    # endTry
# endDef
    
# Plot ATL08 Button
btn = tk.Button(tab2, text='ADD TO PLOT', font=('Arial Bold', 16), width = 15, command=plotAtl08) 
btn.place(x=335, y=330)


###############################################################################
#
# TAB 2: PLOT OPTIONS - TRUTH
#
###############################################################################

### Plot Truth Panel label
truthPlotLabelframe = tk.LabelFrame(tab2, width=545, height=130, text='Add Reference Data', font=('Arial Bold', 14))
truthPlotLabelframe.place(x=580, y=10)

# Truth Data Y Axis Combo Box
truthText = ['Reference data can only be plotted\n' \
             'when ATL03 x/y values are set to:\n' \
             'easting, northing, cross-track,\n' \
             'along-track, or height.']
    
        
lbl = tk.Label(truthPlotLabelframe, text=truthText[0], font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=8)

# Itialize Truth data plot lists
plotList_truth = ('time', 'latitude', 'longitude', \
                  'easting', 'northing', \
                  'crossTrack', 'alongTrack', \
                  'z', 'class', 'confidence')

plotList_truthNames = ('Reference Time', 'Reference Lat', 'Reference Lon', \
                       'Reference Easting', 'Reference Northing', \
                       'Reference Cross-Track', 'Reference Along-Track', \
                       'Reference Height', 'Reference Class', 'Reference Conf') 

# Plot Truth Button Callback
def plotTruth():
    
    # Try plot code
    try:
                
        # Get
        gtNumToPlot = gtNumPlotBox.current()
        
        # Get x,y combo box number selections
        xVarNum = xValsBox.current()
        yVarNum = yValsBox.current()
        
        # x,y param names
        xParam = plotList_truth[xVarNum]
        yParam = plotList_truth[yVarNum]
        
        if('time' in xParam or 'time' in yParam):
            messagebox.showinfo('Error','Reference data cannot plot time.')
        elif('latitude' in xParam or 'latitude' in yParam):
            messagebox.showinfo('Error','Reference data cannot plot latitude.')
        elif('longitude' in xParam or 'longitude' in yParam):
            messagebox.showinfo('Error','Reference data cannot plot longitude.')
        elif('class' in xParam or 'class' in yParam):
            messagebox.showinfo('Error','Reference data cannot plot classification.')
        elif('confidence' in xParam or 'confidence' in yParam):
            messagebox.showinfo('Error','Reference data cannot plot signal confidence.')
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
    
    except:
        
        print('Cannot plot data. Please check inputs')
        
    # endTry
# endDef

# Plot Truth Button
btn = tk.Button(truthPlotLabelframe, text='ADD TO PLOT', font=('Arial Bold', 16), width = 15, command=plotTruth) 
btn.place(x=320, y=25)


###############################################################################
#
# TAB 2: PLOT OPTIONS - CORRECTED MEASURED
#
###############################################################################

### Plot Corrected Measured Panel label
measCorrPlotLabelframe = tk.LabelFrame(tab2, width=545, height=130, text='Add Shifted ICESat-2 Data', font=('Arial Bold', 14))
measCorrPlotLabelframe.place(x=580, y=150)

# Truth Data Y Axis Combo Box
measCorrText = ['Shifted ICESat-2 data can only be\n' \
                'plotted when ATL03 x/y values are set to:\n' \
                'easting, northing, cross-track,\n' \
                'along-track, or height.']

# Corrected Measured Data Y Axis Combo Box
lbl = tk.Label(measCorrPlotLabelframe, text=measCorrText[0], font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=8)

plotList_measCorr = ('time', 'latitude', 'longitude', \
                  'easting', 'northing', \
                  'crossTrack', 'alongTrack', \
                  'z', 'class', 'confidence')

plotList_measCorrNames = ('Shifted ATL03 Time', 'Shifted ATL03 Lat', 'Shifted ATL03 Lon', \
                       'Shifted ATL03 Easting', 'Shifted ATL03 Northing', \
                       'Shifted ATL03 Cross-Track', 'Shifted ATL03 Along-Track', \
                       'Shifted ATL03 Height', 'Shifted ATL03 Class', 'Shifted ATL03 Conf') 

# Plot Corrected Measured Button Callback
def plotMeasCorr():
    
    # Try plot code
    try:
                
        # Get
        gtNumToPlot = gtNumPlotBox.current()
        
        # Get x,y combo box number selections
        xVarNum = xValsBox.current()
        yVarNum = yValsBox.current()
        
        # x,y param names
        xParam = plotList_truth[xVarNum]
        yParam = plotList_truth[yVarNum]
        
        if('time' in xParam or 'time' in yParam):
            messagebox.showinfo('Error','Corrected measured data cannot plot time.')
        elif('latitude' in xParam or 'latitude' in yParam):
            messagebox.showinfo('Error','Corrected measured data cannot plot latitude.')
        elif('longitude' in xParam or 'longitude' in yParam):
            messagebox.showinfo('Error','Corrected measured data cannot plot longitude.')
        elif('class' in xParam or 'class' in yParam):
            messagebox.showinfo('Error','Corrected measured data cannot plot classification.')
        elif('confidence' in xParam or 'confidence' in yParam):
            messagebox.showinfo('Error','Corrected measured data cannot plot signal confidence.')
        else: 
            
            # Get x,y combo bxx text selections
            xDataOrig = eval('atl03Data[' + str(gtNumToPlot) + '].' + xParam)
            yDataOrig = eval('atl03Data[' + str(gtNumToPlot) + '].' + yParam)
            xDataCorr = eval('atlCorrections[' + str(gtNumToPlot) + '].' + xParam)
            yDataCorr = eval('atlCorrections[' + str(gtNumToPlot) + '].' + yParam)
            fileName = eval('atl03Data[' + str(gtNumToPlot) + '].atl03FileName')
            gtNum = eval('atl03Data[' + str(gtNumToPlot) + '].gtNum')
        
            xData = xDataOrig + xDataCorr
            yData = yDataOrig + yDataCorr
            
            # Get labels
            xLabel = xlabel_textBox.get()
            yLabel = ylabel_textBox.get()
            title = fileName + ' (' + gtNum + ')'
            
            yName = plotList_measCorrNames[yVarNum]
            
            # Call getPlot function
            getPlot_measCorr(xData, yData, xLabel, yLabel, title, yName)
    
    except:
        
        print('Cannot plot data. Please check inputs')
        
    # endTry
# endDef
    
# Plot Corrected Measured Button
btn = tk.Button(measCorrPlotLabelframe, text='ADD TO PLOT', font=('Arial Bold', 16), width = 15, command=plotMeasCorr) 
btn.place(x=320, y=25)


###############################################################################
#
# TAB 2: PLOT OPTIONS - PHO SHOW
#
###############################################################################

### Plot Corrected Measured Panel label
phoShowPlotLabelframe = tk.LabelFrame(tab2, width=545, height=100, text='PhoSHOW', font=('Arial Bold', 14))
phoShowPlotLabelframe.place(x=580, y=290)

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
        
        outFilePath = outPath_textBox.get().strip()
                
        # Update message
        print('Writing HTML File...')
        
        # Call PhoShow code
        createHTMLChart(ytrack, h_ph, classification, lat, lon, 
                    direction,
                    online = None,
                    classification_list = [1,2,3],
                    output_folder = outFilePath,
                    in_file03_name = 'ATL03', 
                    blank = viewerBlank_html, 
                    blanki = viewerBlankOnline_html)
        
        # Update message
        print('HTML File Complete!')
        
        # Open PhoShow in HTML browser
        webbrowser.open_new_tab('Viewer_ATL03.html')
    
    except:
        
        # Print error message
        print('Cannot open PhoSHOW. Please check inputs')
        
    # endTry
# endDef
    
# PhoShow Button
btn = tk.Button(phoShowPlotLabelframe, text='PhoSHOW', font=('Arial Bold', 16), width = 15, command=phoShow) 
btn.place(x=30, y=15)


###############################################################################
#
# TAB 3: HELP PAGE
#
###############################################################################

### Plot Corrected Measured Panel label
helpLabelframe1 = tk.LabelFrame(tab3, width=545, height=450, text='PhoREAL Help Information', font=('Arial Bold', 14))
helpLabelframe1.place(x=15, y=10)

phoreal_help_info1 = \
'Get ICESat-2 Data Input Section:\n' \
'--------------------------------------------\n' \
'This section handles the input ATL03/ATL08 files from ICESat-2\n' \
'\n' \
'+ ATL03 File: Path to input ATL03 .h5 file\n' \
'+ ATL08 File: Path to input ATL08 .h5 file\n' \
'+ Output Directory: Path to directory for output data\n' \
'+ Ground Track Numbers: Option to select ICESat-2 ground tracks\n' \
'+ Trim ICESat-2 Data Options: Methods to trim ICESat-2 data\n' \
'     - None: No trimming\n' \
'     - Manual: Trim to user-specified latitude or time min/max bounds\n' \
'     - Auto: Trim to reference region bounds (ARL ONLY)\n' \
'+ Create Output Files: Option to create output files for measured data\n' \
'\n' \
'Get Reference Data Input Section:\n' \
'-----------------------------------------------\n' \
'This section handles the reference data used to find ICESat-2 offsets\n' \
'\n' \
'+ Use Existing Data: Option to use existing or new reference data\n' \
'+ Buffer Size: For creating new reference data (ARL Only)\n' \
'+ Reference File Name: Path to reference file (.las, .laz, or .tif files)\n' \
'+ Create Reference File: Option to create output reference file'


lbl = tk.Label(helpLabelframe1, text=phoreal_help_info1, font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=10)

### Plot Corrected Measured Panel label
helpLabelframe2 = tk.LabelFrame(tab3, width=545, height=440, text='', font=('Arial Bold', 14))
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
'+ Grid Resolution(s) (m): Raster resolution(s) to grid reference data\n' \
'+ Use Fixed Vertical Shift: Option to use a fixed vertical shift value\n' \
'+ Vertical Shift (m): Vertical shift value if previous option is selected\n' \
'+ Use ICESat-2 Signal Confidence Value(s): Option to use ICESat-2\n' \
'   signal confidence values to filter measured data\n' \
'     - Input signal confidence values to filter ICESat-2 data\n' \
'+ Use Reference Ground Index: Option to use reference ground index value\n' \
'   to filter reference data (Requires ATL08 file)\n' \
'     - Reference ground index value (Texpert ground class = 2)\n' \
'+ Create Shifted ICESat-2 File: Option to create corrected file\n' \
'+ Make Output Plots: Option to create output plots'

lbl = tk.Label(helpLabelframe2, text=phoreal_help_info2, font=('Arial', 12), anchor = 'w', justify='left')
lbl.place(x=10, y=10)


###############################################################################
#
# TAB 4: ABOUT PAGE
#
###############################################################################

### Panel label
aboutLabelFrame = tk.LabelFrame(tab4, width=1110, height=450, text='', font=('Arial Bold', 14))
aboutLabelFrame.place(x=15, y=10)

aboutInfo = \
'PhoREAL Description\n' \
'-------------------------------------------------------------------------------------------------------------\n' \
'PhoREAL is a comprehensive tool for analyzing and plotting ATL03 and ATL08\n' \
'data from ICESat-2. All of the capabilities of PhoREAL v3.0\n' \
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
'Copyright 2019 Applied Research Laboratories, The University of Texas at Austin'

lbl = tk.Label(aboutLabelFrame, text=aboutInfo, font=('Arial', 12), anchor = 'w', justify='center')
lbl.place(x=260, y=10)

# Open GUI window
window.mainloop()