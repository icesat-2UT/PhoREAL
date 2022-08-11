# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 08:42:41 2022

@author: ERICG
"""

import pandas as pd
import numpy as np
import ntpath
from osgeo import gdal, osr
from gdalconst import GA_ReadOnly
import os 
import glob
from icesatUtils import transform
# from icesatUtils import (identifyEPSG, getCoordRotFwd, transform, getGeoidHeight, \
#                          getCoordRotRev, superFilter, getRaster, getUTM2LatLon)

### Function to get truth file paths from GUI input
def get_truth_file_paths(input_f, file_ext):
    
    # INPUTS:
    # input_f - string, list, tuple of path(s)
    # file_ext - '.las', '.tif', '.h5' etc.
    # logFileID - a file identifier to write to a file if desired
    
    # Initialize output
    file_list = []
    
    # Convert tuple to list if necessary
    if(isinstance(input_f, tuple)):
        input_f = list(input_f)
    
    elif(isinstance(input_f, list)):
        if(len(input_f)==1):
            input_f = input_f[0]

    # Get truth file paths based on file, files, or directory
    if(isinstance(input_f, str)):
        if(os.path.exists(input_f)):
            if(os.path.isfile(input_f)):
                # Get truth file name and extension
                file_list = [os.path.normpath(input_f)]
            else:
                # Get full paths to input files
                if os.name == 'nt':
                    strSearch = os.path.normpath(input_f + '\\*' + file_ext)
                else:
                    strSearch = os.path.normpath(input_f + '/*' + file_ext)
                file_list = glob.glob(strSearch, recursive = True)   

    elif(isinstance(input_f, list)):
        file_list_list = [os.path.normpath(x) for x in input_f]
        file_list_true = all([os.path.isfile(x) for x in file_list_list])
        if(file_list_true):
            file_list = file_list_list

    return file_list

### Function reads EPSG code from Tif File
def read_dem_epsg(file):
    ds = gdal.Open(file)
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    proj.AutoIdentifyEPSG()
    epsg = proj.GetAttrValue('AUTHORITY',1)
    return epsg

### Function to find and reproject min/max extents in header data
def reproject_header_data(truthHeaderDF, epsg_atl):
    '''
    Description:
        Reproject header data 
    
    Input:
        truthheader_df - Pandas DF that has the location for each truth file
        epsg_atl - str epsg code for ATL file, e.g., 'epsg:32618'
        
    Output:
        truthheader_df - Reprojected DF that has the location for each file
    
    '''
    # Copy dataframe
    # truthHeaderNewDF = truthHeaderDF.copy()
    
    # Input EPSG code
    #epsg_atl = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
    
    # Loop through all truth file EPSG codes
    for i in range(0,len(truthHeaderDF)):
#        writeLog(str(i) + ' out of ' + str(len(truthHeaderDF)))
        # Reproject EPSG code to match input EPSG if unmatching
        epsg_truth = truthHeaderDF['epsg'][i]
#        writeLog(epsg_truth, logFileID)
#        writeLog(epsg_atl, logFileID)
        if(epsg_truth!=epsg_atl):
            
            # Header extents
            xmin = truthHeaderDF['xmin'][i]
            xmax = truthHeaderDF['xmax'][i]
            ymin = truthHeaderDF['ymin'][i]
            ymax = truthHeaderDF['ymax'][i]
            x = [xmin, xmax]
            y = [ymin, ymax]
            
            # Reproject extents to input EPSG code
            try:

                xout, yout = transform(epsg_truth, epsg_atl, x, y)
                
                # Store new extents into output dataframe
                truthHeaderDF['xmin'][i] = xout[0]
                truthHeaderDF['xmax'][i] = xout[1]
                truthHeaderDF['ymin'][i] = yout[0]
                truthHeaderDF['ymax'][i] = yout[1]
                truthHeaderDF['epsg'][i] = epsg_atl
            
            except:
                
                # Store new extents into output dataframe
                print('WARNING: Cannot reproject data, skipping file: %s' %truthHeaderDF['fileName'][i])
                truthHeaderDF['xmin'][i] = 'None'
                truthHeaderDF['xmax'][i] = 'None'
                truthHeaderDF['ymin'][i] = 'None'
                truthHeaderDF['ymax'][i] = 'None'
                truthHeaderDF['epsg'][i] = 'None'
            
            # endTry
                
            
        # endIf
    # endFor
    
    return truthHeaderDF


### Function to check file inputs
def get_files_from_input(fileInput, truthFileType):
    
    # Set input file check
    inputFileCheck = False
    
    # Determine input type
    if(isinstance(fileInput, str)):
        fileInput = os.path.normpath(fileInput)
        if('*' in fileInput):
            # Wilcard case
            file_ext = fileInput[-4:]
            if(truthFileType in file_ext):
                inpFiles = glob.glob(fileInput, recursive = True) 
                inputFileCheck = True
            else:
                print('ERROR: File search must contain *%s in string.' %truthFileType)
        else:
            # File/Directory case
            if(os.path.exists(fileInput)):
                if(os.path.isdir(fileInput)):
                    # Directory case
                    strSearch = os.path.normpath(fileInput + '\\*' + truthFileType)
                    inpFiles = glob.glob(strSearch, recursive = True) 
                    inputFileCheck = True
                elif(os.path.isfile(fileInput)):
                    # File case
                    file_ext = fileInput[-4:]
                    if(truthFileType in file_ext):
                        inpFiles = [fileInput]
                        inputFileCheck = True
                    else:
                        print('ERROR: Input file must be %s format' %truthFileType)
            else:
                print('ERROR: File/Directory does not exist.')
    elif(isinstance(fileInput, list)):
        inpFiles = [os.path.normpath(x) for x in fileInput]
        inputFileCheck = True

    return inpFiles, inputFileCheck

### Function to read header info for LAS files
def read_las_header(lasFileInput):
    
    """
    This function mimics the MATLAB LASreadheader.m function and 
    retrieves .las file header info (for LAS Formats 1.1 - 1.4).
    
    INPUTS:
        1) lasFileInput - contains path to input .las file(s)
           This input can be a:
               - String of 1 file path to a .las file
               - String of 1 file path to a directory with .las files
               - String of 1 file path and a * wildcard for .las files
               - Python list of multiple file paths to .las files
    
    OUTPUTS:
        1) dataOutDF - Pandas dataframe containing las info metadata
           Contents of dataframe:
               - fileName = Input .las File Name
               - version = LAS Version
               - xmin = Point cloud xmin
               - xmax = Point cloud xmax
               - ymin = Point cloud ymin
               - ymax = Point cloud ymax
               - zmin = Point cloud zmin
               - zmax = Point cloud zmax 
               - nPoints = Point cloud number of points
               - epsg = EPSG code
                           
    """
    
    # Import libraries
    import struct
    import numpy as np
    import pandas as pd
    import os

    # Get input .las files        
    if isinstance(lasFileInput,str):
        if lasFileInput[-4:] == '.laz':
            lasFiles, inputFileCheck = get_files_from_input(lasFileInput, '.laz')
        elif lasFileInput[-4:] == '.las':
            lasFiles, inputFileCheck = get_files_from_input(lasFileInput, '.las')

    if isinstance(lasFileInput,list):
        if lasFileInput[0][-4:] == '.laz':
            lasFiles, inputFileCheck = get_files_from_input(lasFileInput, '.laz')
        elif lasFileInput[0][-4:] == '.las':
            lasFiles, inputFileCheck = get_files_from_input(lasFileInput, '.las')
        
    # Set dataframe column names
    columnNames = ['filename','version','xmin','xmax','ymin','ymax','zmin','zmax', 
                   'nPoints','epsg']
    
    # Set dataframe
    dataOutDF = pd.DataFrame(columns=columnNames)
    
    if(inputFileCheck):
        
        # Loop through all files
        for numFile in range(0,len(lasFiles)):
            
            # Get input file
            lasFile = lasFiles[numFile]
        
            # Get file name from path
            fileName = lasFile
            
            # Initialize output parameters
            version = np.nan
            xmin = np.nan
            xmax = np.nan
            ymin = np.nan
            ymax = np.nan
            zmin = np.nan
            zmax = np.nan
            Npoints = np.nan
            epsg = 'None'
            
            # Open binary .las file for reading
            with open(lasFile, mode='rb') as file:
                
                # Move forward in file from beginning
                file.seek(24,0)
                    
                # Get parameters
                VersionMajor = struct.unpack('B',file.read(1))[0] # uint8
                VersionMinor = struct.unpack('B',file.read(1))[0] # uint8
                
                # Store LAS version
                version = 'LAS v' + str(VersionMajor) + '.' + str(VersionMinor)
                
                # Move forward in file form beginning
                file.seek(94,0)
                
                # Get parameters
                headerSize = struct.unpack('H',file.read(2))[0] # uint16
                OffsetToPointData = struct.unpack('I',file.read(4))[0] # uint32
                nvlr = struct.unpack('I',file.read(4))[0] # uint32
                recordFormat = struct.unpack('B',file.read(1))[0] # uint8
                recordLength = struct.unpack('H',file.read(2))[0] # uint16
                Npoints = struct.unpack('I',file.read(4))[0] # uint32
                
                # Move forward in file from current position
                file.seek(20, 1)
                
                # Get parameters
                XScaleFactor = struct.unpack('d',file.read(8))[0] # double
                YScaleFactor = struct.unpack('d',file.read(8))[0] # double
                ZScaleFactor = struct.unpack('d',file.read(8))[0] # double
                XOffset = struct.unpack('d',file.read(8))[0] # double
                YOffset = struct.unpack('d',file.read(8))[0] # double
                ZOffset = struct.unpack('d',file.read(8))[0] # double
                xmax = struct.unpack('d',file.read(8))[0] # double
                xmin = struct.unpack('d',file.read(8))[0] # double
                ymax = struct.unpack('d',file.read(8))[0] # double
                ymin = struct.unpack('d',file.read(8))[0] # double
                zmax = struct.unpack('d',file.read(8))[0] # double
                zmin = struct.unpack('d',file.read(8))[0] # double
                    
                if(np.isin(VersionMinor,[0,1,2,3])):
                    
                    # LAS 1.0 - 1.3 Format
        
                    # Move forward in file from beginning
                    file.seek(headerSize, 0)
                    
                    # Initialize arrays
                    uID = np.zeros((1,16))
                    desc = np.zeros((1,32))
                    
                    for i in range(0,nvlr):
                        
                        # Move forward in file from current position
                        file.seek(2, 1)
                        
                        for j in range(0,16):
                            uID[0,j] = struct.unpack('B',file.read(1))[0] # uint8
                        
                        
                        # Convert userID to ASCII char
                        userID = ''.join([chr(int(item)) for item in uID[0]])
                        
                        recordID = struct.unpack('H',file.read(2))[0] # uint16
                        recordLen = struct.unpack('H',file.read(2))[0] # uint16
                    
                        for j in range(0,32):
                            desc[0,j] = struct.unpack('B',file.read(1))[0] # uint8
                        
                    
                        if('LASF_Projection' in userID):
                            
                            if(recordID == 34735):
                                
                                sGeoKeys_wKeyDirectoryVersion = struct.unpack('H',file.read(2))[0] # uint16
                                sGeoKeys_wKeyRevision = struct.unpack('H',file.read(2))[0] # uint16
                                sGeoKeys_wMinorRevision = struct.unpack('H',file.read(2))[0] # uint16
                                sGeoKeys_wNumKeys = struct.unpack('H',file.read(2))[0] # uint16
                                
                                sGeoKeys_wKeyID = np.zeros((sGeoKeys_wNumKeys,1))
                                sGeoKeys_wTIFFTagLocation = np.zeros((sGeoKeys_wNumKeys,1))
                                sGeoKeys_wCount = np.zeros((sGeoKeys_wNumKeys,1))
                                sGeoKeys_wValue_Offset = np.zeros((sGeoKeys_wNumKeys,1))
                                
                                for k in range(0,sGeoKeys_wNumKeys):
                                    
                                    sGeoKeys_wKeyID[k,0] = struct.unpack('H',file.read(2))[0] # uint16
                                    sGeoKeys_wTIFFTagLocation[k,0] = struct.unpack('H',file.read(2))[0] # uint16
                                    sGeoKeys_wCount[k,0] = struct.unpack('H',file.read(2))[0] # uint16
                                    sGeoKeys_wValue_Offset[k,0] = struct.unpack('H',file.read(2))[0] # uint16
                                    
                                    if(sGeoKeys_wKeyID[k,0] == 1024): # Lat/Lon file
                                        
                                        if(sGeoKeys_wTIFFTagLocation[k,0] == 0):
                                            
                                            coordType = int(sGeoKeys_wValue_Offset[k,0])
                                            keyId = struct.unpack('H',file.read(2))[0] # uint16
                                            tagLoc = struct.unpack('H',file.read(2))[0] # uint16
                                            nCount = struct.unpack('H',file.read(2))[0] # uint16
                                            geoCScode = struct.unpack('H',file.read(2))[0] # uint16
                                            
                                            # Get EPSG code
                                            epsg = 'epsg:' + str(int(geoCScode))
                                                                                
                                    elif(sGeoKeys_wKeyID[k,0] == 3072): # UTM file
                                        
                                        if(sGeoKeys_wTIFFTagLocation[k,0] == 0):
                                            projCScode = sGeoKeys_wValue_Offset[k,0]
                                        else :
                                            projCScode = sGeoKeys_wTIFFTagLocation[k,0]
                                        
                                        # Get EPSG code
                                        epsg = 'epsg:' + str(int(projCScode))                                       
                                                    
                        else:
                            
                            # Move to specified point in file
                            file.seek(recordLen, 1)
                    
                else:
                    
                    # LAS 1.4 Format
                    
                    # Move forward in file from current position
                    file.seek(20, 1)
                        
                    # Get parameters
                    Npoints = struct.unpack('Q',file.read(8))[0] # double
                    
                    # Move forward in file from beginning
                    file.seek(headerSize, 0);
                    
                    # Initialize arrays
                    uID = np.zeros((1,16))
                    desc = np.zeros((1,32))
                    
                    for i in range(0,nvlr):
                        
                        # Move forward in file from current position
                        file.seek(2, 1)
                        
                        for j in range(0,16):
                            uID[0,j] = struct.unpack('B',file.read(1))[0] # uint8
                        
                        
                        # Convert userID to ASCII char
                        userID = ''.join([chr(int(item)) for item in uID[0]])
                        
                        recordID = struct.unpack('H',file.read(2))[0] # uint16
                        recordLen = struct.unpack('H',file.read(2))[0] # uint16
                        
                        for j in range(0,32):
                            desc[0,j] = struct.unpack('B',file.read(1))[0] # uint8
                        
                        
                        # Convert description to ASCII char
                        description = ''.join([chr(int(item)) for item in desc[0]])
                        
                        wkt = np.zeros((1,recordLen))
                        for j in range(0,recordLen):
                            wkt[0,j] = struct.unpack('B',file.read(1))[0] # uint8
                        
                        
                        # Convert wkt to ASCII char
                        wktAll = ''.join([chr(int(item)) for item in wkt[0]])
                        
                        # Get EPSG code
                        proj = osr.SpatialReference(wkt=wktAll)
                        try:
                            epsg = 'epsg:' + proj.GetAttrValue('AUTHORITY',1)
                        except:
                            pass
    
#                        # Get ellipsoid
#                        ellipsoid = (findStr(wktAll,'PROJCS["','/ UTM zone')).strip()
#                        ellipsoid = ''.join(ellipsoid.split())
#                        
#                        # Get UTM Zone/Hemi
#                        utmInfo = (findStr(wktAll,'/ UTM zone','",GEOGCS["WGS')).strip()
#                        zone = str(int(utmInfo[0:-1])).zfill(2)
#                        hemi = utmInfo[-1:]
                        
                    
                
            
            
#            # Define EPSG code
#            if(hemi == 'N'):
#                epsg = 'epsg:326'
#            elif(hemi == 'S'):
#                epsg = 'epsg:327'
#            
#            epsg = epsg + zone
            
            # Set output array
            dataOut = [fileName, version, xmin, xmax, ymin, ymax, zmin, zmax,
                       Npoints, epsg]
            
            # Set dataframe
            dataOutDFsingle = pd.DataFrame(data=dataOut).T
            dataOutDFsingle.columns = columnNames
            
            # Append dataframe
            dataOutDF = dataOutDF.append(dataOutDFsingle, ignore_index=True)  
    
    return dataOutDF

### Function to read header info for TIF files
def read_tif_header(tifFileInput):
    
    """
    This function pulls the header info for .tif files. This uses Python
    binding for GDAL which have an input .tif file size limit. To avoid this,
    the GDAL exe will have to integrated into PhoREAL and used.
    
    INPUTS:
        1) tifFileInput - contains path to input .tif file(s)
           This input can be a:
               - String of 1 file path to a .tif file
               - String of 1 file path to a directory with .tif files
               - String of 1 file path and a * wildcard for .tif files
               - Python list of multiple file paths to .tif files
               
    
    OUTPUTS:
        1) dataOutDF - Pandas dataframe containing tif info metadata
           Contents of dataframe:
               - fileName = Input .tif File Name
               - version = NaN for .tif files (populated for .las files)
               - xmin = Tif file xmin
               - xmax = Tif file xmax
               - ymin = Tif file ymin
               - ymax = Tif file ymax
               - zmin = Tif file zmin
               - zmax = Tif file zmax 
               - nPoints = NaN for .tif files (populated for .las files)
               - epsg = EPSG code
               
=            
    """
    
    # Import libraries

    # Get input .tif files
    tifFiles, inputFileCheck = get_files_from_input(tifFileInput, '.tif')
        
    # Set dataframe column names
    columnNames = ['filename','version','xmin','xmax','ymin','ymax','zmin','zmax', 
                   'npoints','epsg']
    
    # Set dataframe
    dataOutDF = pd.DataFrame(columns=columnNames)
    
    if(inputFileCheck):
        
        # Loop through all files
        for numFile in range(0,len(tifFiles)):
            
            # Get input file
            tifFile = tifFiles[numFile]
        
            # Get file name from path
            fileName = tifFile
            
            # Initialize output parameters
            version = np.nan
            xmin = np.nan
            xmax = np.nan
            ymin = np.nan
            ymax = np.nan
            zmin = np.nan
            zmax = np.nan
            Npoints = np.nan
            epsg = 'None'
            
            # Get x/y min/max extents            
            data = gdal.Open(fileName, GA_ReadOnly)
            geoTransform = data.GetGeoTransform()
            xmin = geoTransform[0]
            ymax = geoTransform[3]
            xmax = xmin + geoTransform[1]*data.RasterXSize
            ymin = ymax + geoTransform[5]*data.RasterYSize
            
            # Get EPSG code
            epsg = 'epsg:' + read_dem_epsg(fileName)
            
            # Set output array
            dataOut = [fileName, version, xmin, xmax, ymin, ymax, zmin, zmax,
                       Npoints, epsg]
            
            # Set dataframe
            dataOutDFsingle = pd.DataFrame(data=dataOut).T
            dataOutDFsingle.columns = columnNames
            
            # Append dataframe
            dataOutDF = dataOutDF.append(dataOutDFsingle, ignore_index=True)
            
    return dataOutDF


def get_truth_headers(truth_filepath, truth_filetype):
    # Initialize output
    truthheader_df = []
    writenew = True
    appendnew = False
        
    # Create header file output name
    header_filename = 'phoReal_headers.csv'
    truth_filedir = ntpath.dirname(truth_filepath[0])
    header_filepath = os.path.join(truth_filedir, header_filename)

    # Check if header file exists
    if(os.path.exists(header_filepath)):
        stored_truthheader_df = pd.read_csv(header_filepath)
        stored_truth_filepath = (stored_truthheader_df['filename']).to_list()
        storedfiles_true = np.isin(truth_filepath, stored_truth_filepath)
        storedfiles_same = all(storedfiles_true)
        
        # Determine if files are the same
        if(storedfiles_same):
            if(len(truth_filepath)==len(stored_truth_filepath)):
                truthheader_df = stored_truthheader_df
                truthheader_df.reset_index(inplace=True)
            elif(len(truth_filepath)<len(stored_truth_filepath)):
                indsTrue = np.isin(stored_truth_filepath, truth_filepath)
                truthheader_df = stored_truthheader_df[indsTrue]
                truthheader_df.reset_index(inplace=True)
            writenew = False
        else:
            newInds = storedfiles_true==False
            truth_filepath = (np.array(truth_filepath)[newInds]).tolist()
            appendnew = True
    
    if(writenew):
        # Message to user
        num_files = len(truth_filepath)
        print('   Reading header info for all %d reference files in: %s' %\
              (num_files, truth_filedir))

    
        if(('las' in truth_filetype.lower()) or ('laz' in truth_filetype.lower())):  
            # Get .las header info
            truthheader_df_new = read_las_header(truth_filepath)              
        elif('tif' in truth_filetype):
            # Load .tif file
            truthheader_df_new = read_tif_header(truth_filepath)
        

        if(os.path.exists(header_filepath)):
            print('   Updating output header file: %s' %header_filepath)
            truthheader_df_new.to_csv(header_filepath, mode='a', header=False, index=False)
        else:
            print('   Writing output header file: %s' %header_filepath)
            truthheader_df_new.to_csv(header_filepath, index=False)
        
        # Append dataframes together if necessary
        if(appendnew):
            frames = [stored_truthheader_df, truthheader_df_new]
            truthheader_df = pd.concat(frames, ignore_index=True)
            truthheader_df.reset_index(inplace=True)

        else:
            truthheader_df = truthheader_df_new
            
    return truthheader_df

def main():
    # Get Generate Truth Header
    truth_filepath = 'E:/data/2018spl_2959_6647'
    truth_filetype = '.laz'
    epsg_atl = 'epsg:32618'
    truth_file_list = get_truth_file_paths(truth_filepath, truth_filetype)
    truthheader_df = get_truth_headers(truth_file_list, truth_filetype)
    
    # Reproject Truth Header to be in the same projection as Measured Data 
    truthheader_df_new = reproject_header_data(truthheader_df, epsg_atl)
    
    # # Find truth files that intersect ICESat-2 track
    # print('   Determining which reference files intersect ground track...')
    # _, matchingTruthFileInds = findMatchingTruthFiles(truthheader_df_new, 
    #                                                   atlMeasuredData, rotationData, buffer)
    # matchingTruthFiles = np.array(truth_filepaths)[matchingTruthFileInds]
    
    
    # Loop over matching files, read files, put into Truth Pandas DF
    
    
    
    # Export as .las
    
    # Export as .csv/.hdf
    
    
    
if __name__ == "__main__":
    main()
   