# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 09:07:22 2021

This script creates the PhoREAL GUI exe along with the Windows Installer file.

Steps:
1) Place all necessary PhoREAL files into the path given by the 
   "phorealDirPath" parameter
   
2) Set the PhoREAL version number for the naming of the exe
    
3) Run this script
    
    
If Pyinstaller fails to compile the Python code, check the generated .spec file:
    + hiddenimports - you may need to include packages here, the error will tell you
    + setrecursionlimit - you may need to up this limit

@author: malonzo
"""

import os
import subprocess
import glob
import shutil
import numpy as np
from datetime import datetime

### SETUP INPUT PARAMETERS

# Set your Conda environment name which has all PhoREAL Python libraries
condaEnvName = 'py36'

# Set path with all PhoREAL files to be compiled into Windows Installer EXE
phorealDirPath = 'C:/Users/malonzo/GLAM/Python/PhoREAL_GUI/GUI_3.28'

# Set PhoREAL Version # for EXE
phorealVersion = '3.28'

# Path to InnoSetup EXE to build Windows Installer EXE
issExe = 'C:/Program Files (x86)/Inno Setup 6/ISCC.exe'

# Set non-Python file extensions to include in Windows Installer build
fileExts = ['.exe','.pyd','.dll','.ico','.zip','.manifest','.so','.html','.txt']

# Set non-Python directories to include in Windows Installer build
dirs = ['d3']

# Auto-run Windows Installer setup .exe after build (requires smart card credentials)
runSetupExe = True


### CODE BELOW ----------------------------------------------------------------

print('\nScript is running...\n')

# Set Hidden Imports needed by Pyinstaller
hidden_imports = [
        'pyproj._compat',
        'pyproj._datadir',
        'pyproj.datadir',
        'rasterio.crs'
        'rasterio._shim',
        'rasterio.compat',
        'rasterio.coords',
        'rasterio.crs',
        'rasterio.dtypes',
        'rasterio.enums',
        'rasterio.env',
        'rasterio.errors',
        'rasterio.features',
        'rasterio.fill',
        'rasterio.mask',
        'rasterio.merge',
        'rasterio.plot',
        'rasterio.profiles',
        'rasterio.sample',
        'rasterio.tool',
        'rasterio.transform',
        'rasterio.vfs',
        'rasterio.warp',
        'rasterio.windows',
        'rasterio.__init__'
        ]


cd = os.getcwd()
phorealDirPathSys = os.path.normpath(phorealDirPath)
specFileName = 'PhoReal_v' + phorealVersion + '.spec'
specFilePath = os.path.normpath(phorealDirPathSys + '/' + specFileName)
batchFileName = 'run_pyinstaller.bat'
batchFilePath = os.path.normpath(cd + '/' + batchFileName) 

# Create pyinstaller .bat file
print('Step 1) Creating Windows .bat file...', end='')

with open(batchFileName,'w') as f:
    f.write('cd ' + phorealDirPathSys + '\n')
    f.write('call activate ' + condaEnvName + '\n')
    f.write('pyinstaller --noconfirm ' + specFileName + '\n')
    f.write('call deactivate')
# endWith
print('Complete.\n') 
print('        Windows .bat file created here: ')
print('        %s\n' %batchFilePath)

# Create Pyinstaller .spec file
print('Step 2) Creating Pyinstaller .spec file...', end='')
with open(specFilePath,'w') as f:
    f.write('# -*- mode: python ; coding: utf-8 -*-\n\n')
    f.write('import sys\n')
    f.write('sys.setrecursionlimit(5000)\n\n')
    f.write('block_cipher = None\n\n\n')
    f.write("a = Analysis(['getAtl03_GUI.py'],\n")
    f.write("             pathex=[" + repr(phorealDirPathSys) + "],\n")
    f.write("             binaries=[],\n")
    f.write("             datas=[('pho_image.ico', '.')],\n")
    f.write("             hiddenimports=" + repr(hidden_imports) + ",\n")
    f.write("             hookspath=[],\n")
    f.write("             runtime_hooks=[],\n")
    f.write("             excludes=[],\n")
    f.write("             win_no_prefer_redirects=False,\n")
    f.write("             win_private_assemblies=False,\n")
    f.write("             cipher=block_cipher,\n")
    f.write("             noarchive=False)\n")
    f.write("pyz = PYZ(a.pure, a.zipped_data,\n")
    f.write("             cipher=block_cipher)\n")
    f.write("exe = EXE(pyz,\n")
    f.write("          a.scripts,\n")
    f.write("          [],\n")
    f.write("          exclude_binaries=True,\n")
    f.write("          name='PhoReal_v" + phorealVersion + "',\n")
    f.write("          debug=False,\n")
    f.write("          bootloader_ignore_signals=False,\n")
    f.write("          strip=False,\n")
    f.write("          upx=True,\n")
    f.write("          console=True , icon='pho_image.ico')\n")
    f.write("coll = COLLECT(exe,\n")
    f.write("               a.binaries,\n")
    f.write("               a.zipfiles,\n")
    f.write("               a.datas,\n")
    f.write("               strip=False,\n")
    f.write("               upx=True,\n")
    f.write("               upx_exclude=[],\n")
    f.write("               name='PhoReal_v" + phorealVersion + "')\n")
# endWith
print('Complete.\n') 
print('        Spec file created here:')
print('        %s\n' %specFilePath)

# Execute .bat file to call Pyinstaller and compile GUI
print('Step 3) Running Pyinstaller .spec file...\n')
process = subprocess.run([batchFilePath])
#process = subprocess.Popen([batchFilePath], shell=True)
print()
    
# Check if Pyinstaller was successful
pythonExeDir = os.path.normpath(phorealDirPathSys + '\\dist\\PhoReal_v' + phorealVersion)
pythonExeDirExists = os.path.isdir(pythonExeDir)
if(pythonExeDirExists):
    print('        Complete.\n') 
    print('        Python EXE files located here:')
    print('        %s\n' %pythonExeDir)
else:
    print('ERROR: Pyinstaller did not complete. Output files not created.\n')
    print('Open a Windows Command shell and run this command:')
    print('pyinstalller PhoReal_v#.##.spec')
    print('For more information on Pyinstaller error messages.\n')
# endIf
   
if(pythonExeDirExists):
    
    # Copy files to Python EXE directory to wrap into Windows Installer
    print('Step 4) Copying necessary files to Python EXE directory...')
    
    # Get files to copy from Python source directory to Python EXE directory  
    sourceFiles = []
    for fileExt in fileExts:
        sourceFiles.extend(glob.glob(os.path.normpath(phorealDirPathSys + '/*' + fileExt)))
    # endFor
    
    sourceDirs = []
    for curDir in dirs:
        sourceDirs.append(os.path.normpath(phorealDirPathSys + '/' + curDir))
    # endForj
       
    # Copy files
    print()
    if(len(sourceFiles)>0):
        print('        Copying Files...')
        for i in range(0,len(sourceFiles)):
            fileName = sourceFiles[i]
            print('        %d) Copying %s to %s' %(i+1,fileName,pythonExeDir))
            shutil.copy(fileName, pythonExeDir)
        # endFor
        print()
    # endIf
    
    # Copy directories
    if(len(sourceDirs)>0):
        print('        Copying directories')
        for i in range(0,len(sourceDirs)):
            sourceDir = sourceDirs[i]
            baseName = os.path.basename(sourceDir)
            print('        %d) Copying %s to %s' %(i+1,baseName,pythonExeDir))
            destDir = os.path.normpath(pythonExeDir + '/' + baseName)
            if(os.path.isdir(destDir)):
                shutil.rmtree(destDir)
            # endIf
            shutil.copytree(sourceDir, destDir)
        # endFor
        print()
    # endIf
    print('Complete.\n')
    
    # Call InnoSetup to wrap all files into a Windows Installer exe
    print('Step 5) Creating InnoSetup .iss file...', end='')
        
    # Get all source files to bundle into Windows Installer
    sourceFiles = []
    for fileExt in fileExts:
        sourceFiles.extend(glob.glob(os.path.normpath(pythonExeDir + '/*' + fileExt)))
    # endFor
    
    # Get all source directories
    sourceDirs = []
    sourceDirs = [os.path.normpath(f.path) for f in os.scandir(pythonExeDir) if f.is_dir()]

    # Initialize output variable
    allSourceFiles = []

    # Loop through all source files and create .iss string
    for i in range(0,len(sourceFiles)):
        newString = 'Source: "' + sourceFiles[i] + '"; DestDir: "{app}"; Flags: ignoreversion'
        allSourceFiles = np.append(allSourceFiles,newString)
    # endFor
    
    # Loop through all source directories and create .iss string
    for i in range(0,len(sourceDirs)):
        dirName = sourceDirs[i].split('\\')[-1]
        newString = 'Source: "' + sourceDirs[i] + '\\*"; DestDir: "{app}\\' + dirName + '"; Flags: ignoreversion recursesubdirs createallsubdirs'
        allSourceFiles = np.append(allSourceFiles,newString)
    # endFor
        
    # Write InnoSetup .iss file
    issFileName = 'PhoREAL_setup.iss'
    issFilePath = os.path.normpath(phorealDirPathSys + '/dist/' + issFileName)
    with open(issFilePath,'w') as f:
        f.write('; Script generated by the Inno Setup Script Wizard.\n')
        f.write('; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!\n\n')

        myAppName = 'PhoREAL_v' + phorealVersion
        f.write('#define MyAppName "' + myAppName + '"\n')
        f.write('#define MyAppVersion "' + phorealVersion + '"\n')
        f.write('#define MyAppPublisher "The Applied Research Laboratories at The University of Texas at Austin"\n')
        f.write('#define MyAppURL "phoreal@arlut.utexas.edu"\n')
        f.write('#define MyAppExeName "PhoReal_v' + phorealVersion + '.exe"\n\n')

        # Create unique AppID (GUID) to identify EXE Build for InnoSetup
        currentDateTime = (datetime.now()).strftime('%m-%d-%Y-%H-%M-%S')
        myAppId = myAppName + '_' + currentDateTime
        f.write('[Setup]\n')
        f.write('; NOTE: The value of AppId uniquely identifies this application. Do not use the same AppId value in installers for other applications.\n')
        f.write('; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)\n')
        f.write('AppId={{' + myAppId + '}\n')
        f.write('AppName={#MyAppName}\n')
        f.write('AppVersion={#MyAppVersion}\n')
        f.write(';AppVerName={#MyAppName} {#MyAppVersion}\n')
        f.write('AppPublisher={#MyAppPublisher}\n')
        f.write('AppPublisherURL={#MyAppURL}\n')
        f.write('AppSupportURL={#MyAppURL}\n')
        f.write('AppUpdatesURL={#MyAppURL}\n')
        f.write('DefaultDirName={autopf}\\{#MyAppName}\n')
        f.write('DisableProgramGroupPage=yes\n')
        f.write('; Uncomment the following line to run in non administrative install mode (install for current user only.)\n')
        f.write(';PrivilegesRequired=lowest\n')
        f.write('OutputBaseFilename=PhoREAL_v' + phorealVersion + '_init\n')
        f.write('SetupIconFile=' + pythonExeDir + '\\pho_image.ico\n')
        f.write('Compression=lzma\n')
        f.write('SolidCompression=yes\n')
        f.write('WizardStyle=modern\n\n')
        
        f.write('[Languages]\n')
        f.write('Name: "english"; MessagesFile: "compiler:Default.isl"\n\n')
        
        f.write('[Tasks]\n')
        f.write('Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked\n\n')
        
        f.write('[Files]\n')
        for i in range(0,len(allSourceFiles)):
            f.write(allSourceFiles[i] + '\n')
        # endFor
        
        f.write('; NOTE: Don\'t use "Flags: ignoreversion" on any shared system files\n\n')
        
        f.write('[Icons]\n')
        f.write('Name: "{autoprograms}\\{#MyAppName}"; Filename: "{app}\\{#MyAppExeName}"\n')
        f.write('Name: "{autodesktop}\\{#MyAppName}"; Filename: "{app}\\{#MyAppExeName}"; Tasks: desktopicon\n\n')
        
        f.write('[Run]\n')
        f.write('Filename: "{app}\\{#MyAppExeName}"; Description: "{cm:LaunchProgram,{#StringChange(MyAppName, \'&\', \'&&\')}}"; Flags: nowait postinstall skipifsilent\n\n')

    # endWith
    print('Complete.\n')    
    print('        InnoSetup .iss file created here:')
    print('        %s\n' %issFilePath)
        
    # Call InnoSetup to wrap all files into a Windows Installer exe
    if(runSetupExe):
        print('Step 6) Calling InnoSetup to create Windows Installer...')
        issExeSys = os.path.normpath('"' + issExe + '"')
        issFilePathSys = '"' + issFilePath + '"'
        issCmd = issExeSys + ' ' + issFilePathSys
        process = subprocess.run(issCmd)
        issExeFile = 'PhoREAL_v' + phorealVersion + '_init.exe'
        issExeDir = os.path.normpath(phorealDirPathSys + '/dist/Output')
        issExePath = os.path.normpath(issExeDir + '/' + issExeFile)
        issFileExists = os.path.isfile(issExePath)
        if(issFileExists):
            print('Complete.\n')
            print()
            print('Running the PhoREAL v' + phorealVersion + ' Windows Installer Setup EXE file...')
            process = subprocess.run([issExePath]) 
            print('\n')
            print('Script is Complete!')
        else:
            print('ERROR: Could not compile with Inno Setup.')
            print('Please use the Inno Setup .iss file through the Inno Setup IDE.')
        # endIf
    else:
        print('User selected not to auto run Windows Installer setup.exe')
    # endIf
        
# endIf