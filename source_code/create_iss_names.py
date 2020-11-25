# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:15:52 2020

@author: malonzo
"""

# Get all necessary source files and directories for installer .exe and put them into 
# format for Inno Setup to wrap them into Windows installer. 
# Inno Setup uses a .iss file format with files and dirs in the format below

# Import libraries
import os
import glob
import numpy as np

# File extensions to include in .iss file for Windows Installer
fileExts = ['.exe','.pyd','.dll','.ico','.zip','.manifest','.so','.html','.mat','.txt']

# Source path for files and directories
#sourcePath = 'C:/Users/malonzo/GLAM/Python/PhoREAL_GUI/GUI_3.18/dist/PhoReal_v3.18'
#sourcePath = 'C:/Users/malonzo/GLAM/Python/PhoREAL_GUI/GUI_3.19/dist/PhoReal_v3.19'
#sourcePath = 'C:/Users/malonzo/GLAM/Python/PhoREAL_GUI/GUI_3.20/dist/PhoReal_v3.20'
#sourcePath = 'C:/Users/malonzo/GLAM/Python/PhoREAL_GUI/GUI_3.21/dist/PhoReal_v3.21'
#sourcePath = 'C:/Users/malonzo/GLAM/Python/PhoREAL_GUI/GUI_3.22/dist/PhoReal_v3.22'
#sourcePath = 'C:/Users/malonzo/GLAM/Python/PhoREAL_GUI/GUI_3.23/dist/PhoReal_v3.23'
sourcePath = 'C:/Users/malonzo/GLAM/Python/PhoREAL_GUI/GUI_3.24/dist/PhoReal_v3.24'

# Get all source files
sourceFiles = []
for fileExt in fileExts:
    sourceFiles.extend(glob.glob(os.path.normpath(sourcePath + '/*' + fileExt)))
# endFor

# Get all source directories
sourceDirs = [os.path.normpath(f.path) for f in os.scandir(sourcePath) if f.is_dir()]

# Initialize output variable
newStringAll = []

# Loop through all source files and create .iss string
for i in range(0,len(sourceFiles)):
    
    newString = 'Source: "' + sourceFiles[i] + '"; DestDir: "{app}"; Flags: ignoreversion'
    newStringAll = np.append(newStringAll,newString)
    
# endFor
    
# Loop through all source directories and create .iss string
for i in range(0,len(sourceDirs)):
    
    dirName = sourceDirs[i].split('\\')[-1]
    newString = 'Source: "' + sourceDirs[i] + '\\*"; DestDir: "{app}\\' + dirName + '"; Flags: ignoreversion recursesubdirs createallsubdirs'
    newStringAll = np.append(newStringAll,newString)
    
# endFor