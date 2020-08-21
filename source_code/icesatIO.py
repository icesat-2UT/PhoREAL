# -*- coding: utf-8 -*-
"""
Script that provides basic I/O functionality for ATL03

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
import os
import sys
import warnings
import csv
import numpy as np
import h5py
try:
  import laspy
  from laspy.file import File
except ImportError:
  laspy_func = ['writeLas', 'writeHeaderFile']
  print('warning: module laspy not found')
  print('affected functions:', laspy_func)

from scipy.io import loadmat
try:
  import simplekml
except ImportError:
  simplekml_func = ['writeKml']
  print('warning: module simplekml not found')
  print('affected functions:', simplekml_func)
  
from osgeo import gdal,osr,ogr
import socket
import shutil
from gui_addins import (viewerBlank_html, viewerBlankOnline_html)



# Object for readKmlBounds function
class kmlStruct:
        
    # Define class with designated fields
    def __init__(self, regionName, lonMin, lonMax, latMin, latMax,headerFilePath, truthFilePath):
            
        self.regionName = regionName
        self.lonMin = lonMin
        self.lonMax = lonMax
        self.latMin = latMin
        self.latMax = latMax
        self.headerFilePath = headerFilePath
        self.truthFilePath = truthFilePath


# Object for readHeaderFile function
class headerStruct:
        
    # Define class with designated fields
    def __init__(self, coordType, zone, hemi, ellipsoid, xmin, xmax, ymin, ymax, tileName):
            
        self.coordType = coordType
        self.zone = zone
        self.hemi = hemi
        self.ellipsoid = ellipsoid
        self.xmin = np.c_[xmin]
        self.xmax = np.c_[xmax]
        self.ymin = np.c_[ymin]
        self.ymax = np.c_[ymax]
        self.tileName = np.c_[tileName]
    
    
# Object for readLas function
class lasStruct:
        
    # Define class with designated fields
    def __init__(self, x, y, z, classification, intensity, headerData):
            
        self.x = np.c_[x]
        self.y = np.c_[y]
        self.z = np.c_[z]
        self.classification = np.c_[classification]
        self.intensity = np.c_[intensity]
        self.headerData = np.c_[headerData]
        
        
# Object for readGeoidFile function
class geoidStruct:
        
    # Define class with designated fields
    def __init__(self, lats, lons, geoidalHeights):
            
        self.lats = lats
        self.lons = lons
        self.geoidalHeights = geoidalHeights
        
        
##### Function to read kmlBounds.txt    
def readTruthRegionsTxtFile(kmlBoundsTextFile):
    
    # Initialize output parameters
    regionName = []
    lonMin = []
    lonMax = []
    latMin = []
    latMax = []
    headerFilePath = []
    truthFilePath = []

    # Open file for reading
    f = open(kmlBoundsTextFile,'r')
    
    # Get all lines of file into a list
    allLines = list(f)
    
    # Close file for reading
    f.close()
    
    # Read lines with path info
    headerPathLine = allLines[2]
    headerPath = headerPathLine.split(',')[1].strip()
    truthPathLine = allLines[3]
    truthPath = truthPathLine.split(',')[1].strip()
        
    # Read rest of lines, skip first six lines (comments)
    textLines = allLines[6:]
            
    # Loop through each line in text file
    for line in textLines:
            
        # If line text exists
        if line:
                
            # Split line by commas
            fields = line.split(',')
            
            # Store field names into lists
            regionName.append(fields[0].strip())
            lonMin.append(float(fields[1].strip()))          
            lonMax.append(float(fields[2].strip()))              
            latMin.append(float(fields[3].strip()))           
            latMax.append(float(fields[4].strip()))           
            headerFilePath.append(headerPath + fields[5].strip())             
            truthFilePath.append(truthPath + fields[6].strip())
            
        # Endif
    # Endfor
    
    # Call class to populate field names
    kmlInfo = kmlStruct(regionName, lonMin, lonMax, latMin, latMax, headerFilePath, truthFilePath)
    
    # Return output
    return kmlInfo


##### Function to read header truth .mat file
def readHeaderMatFile(headerFilePath):

    if headerFilePath.endswith('.mat'):

        # Initialize output data
        coordType = []
        zone = []
        hemi = [] 
        
        # Get coordinate type from header .mat file (1 = UTM, 2 = Lat/Lon)
        matData = loadmat(headerFilePath)
        coordNum= matData['headerData'][0][0][10][0][0]
        
        # Convert lat/lon data to UTM coordinates
        if(coordNum==1):
            
            coordType = 'UTM'
            
            # Get UTM zone and hemisphere
            zone = str(matData['headerData'][0][0][7][0][0])
            hemi = matData['headerData'][0][0][8][0][0]
            ellipsoid = matData['headerData'][0][0][9][0][0:]
            
        else:
            
            coordType = 'Lat/Lon'
            
        # Endif
        
        # Get x/y min/max data and truth tile name for each truth tile
        xmin = [matData['headerData'][0][i][0][0][0] for i in range(len(matData['headerData'][0]))]
        xmax = [matData['headerData'][0][i][1][0][0] for i in range(len(matData['headerData'][0]))]
        ymin = [matData['headerData'][0][i][2][0][0] for i in range(len(matData['headerData'][0]))]
        ymax = [matData['headerData'][0][i][3][0][0] for i in range(len(matData['headerData'][0]))]
        tileName = [matData['headerData'][0][i][12][0] for i in range(len(matData['headerData'][0]))]
        
        # Store data as object
        headerData = headerStruct(coordType, zone, hemi, ellipsoid, xmin, xmax, ymin, ymax, tileName)
        
        # Return data
        return headerData

    elif headerFilePath.endswith('.pkl'):

        from icesatReader import read_pickle
        from icesatUtils import mode, transform
        import pyproj

        df = read_pickle(headerFilePath)

        proj4_all = list(df['proj4'])
        proj4_mode = mode(proj4_all)
        zone0, hemi0, ellipsoid0, coordNum0 = get_params(proj4_mode)
        tileName_v = list(df['name'])

        coordType0 = 'UTM'
        if coordNum0 != 1:
            coordType0 = 'Lat/Lon'

        xmin = np.array(df['xmin'])
        xmax = np.array(df['xmax'])
        ymin = np.array(df['ymin'])
        ymax = np.array(df['ymax'])

        crs_mode = pyproj.CRS.from_proj4(proj4_mode)
        epsg_mode = crs_mode.to_epsg()

        for k, proj4 in enumerate(proj4_all):
            crs_las = pyproj.CRS.from_proj4(proj4)
            epsg_las = crs_las.to_epsg()

            epsg_in = int(epsg_las)
            epsg_out = int(epsg_mode)
            if epsg_in != epsg_out:
                xmin0, ymin0 = transform(epsg_in, epsg_out, xmin[k], ymin[k])
                xmax0, ymax0 = transform(epsg_in, epsg_out, xmax[k], ymax[k])
                xmin[k], ymin[k] = xmin0, ymin0
                xmax[k], ymax[k] = xmax0, ymax0

        xmin_v = list(xmin)
        xmax_v = list(xmax)
        ymin_v = list(ymin)
        ymax_v = list(ymax)

        headerData = headerStruct(coordType0, zone0, hemi0, ellipsoid0, xmin_v, xmax_v, ymin_v, ymax_v, tileName_v)

        return headerData


def get_params(proj4):
    from icesatUtils import identify_hemi_zone
    import pyproj
    crs_las = pyproj.CRS.from_proj4(proj4)
    epsg_las = crs_las.to_epsg()
    hemi, zone = identify_hemi_zone(epsg_las)

    proj = pyproj.Proj(proj4)
    d = proj.crs.to_json_dict()
    ellipsoid = d['base_crs']['datum']['ellipsoid']['name']
    coordType = d['conversion']['name']
    coordNum = -1
    if 'UTM' in coordType:
        coordNum = 1
    # else:
    #     coordType = 'Lat/Lon'
    return zone, hemi, ellipsoid, coordNum


def writeHeaderMatFile(regionName, PKL_DIR, LAS_DIR):
    writeHeaderFile(regionName, PKL_DIR, LAS_DIR)


def writeHeaderFile(regionName, PKL_DIR, LAS_DIR):

    regionName = regionName + '_HeaderData.pkl'
    file_pkl = os.path.join(PKL_DIR, regionName)

    from icesatReader import write_pickle
    from icesatIO import readLas
    import laspy as las
    import pandas as pd

    files_las = [fn for fn in os.listdir(LAS_DIR) if fn.endswith('.las')]
    files_las = sorted(files_las)
    files_las = [os.path.join(LAS_DIR, fn) for fn in files_las]

    df = pd.DataFrame(columns=['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax',
                                        'pointcount', 'zone', 'hemi', 'ellipsoid',
                                        'coordNum', 'numRegs', 'name', 'proj4'])

    for f, file_las in enumerate(files_las):

        file_las_sub = os.path.basename(file_las)
        print(file_las_sub)

        fp = las.file.File(file_las, mode='r')

        min_bounds = fp.get_header().get_min()
        max_bounds = fp.get_header().get_max()
        xb_las = np.array([min_bounds[0], max_bounds[0]])
        yb_las = np.array([min_bounds[1], max_bounds[1]])
        zb_las = np.array([min_bounds[2], max_bounds[2]])

        xmin, xmax = xb_las[0], xb_las[1]
        ymin, ymax = yb_las[0], yb_las[1]
        zmin, zmax = zb_las[0], zb_las[1]

        pointcount = len(fp.points) #readLas(file_las, metadata='count')
        proj4_las = readLas(file_las, metadata='proj4')
        zone, hemi, ellipsoid, coordNum = get_params(proj4_las)

        numRegs = 0
        name = file_las

        df.loc[f] = [xmin, xmax, ymin, ymax, zmin, zmax,
                            pointcount, zone, hemi, ellipsoid, coordNum, numRegs, name, proj4_las]

        fp.close()

    # save header file
    # convert_df_to_mat(df_mat, file_mat)
    write_pickle(df, file_pkl)

    # return df

# def getHeaderData():


##### Functions to read ATL03 .h5 files
def readAtl03H5(in_file03, fieldName, label):
    
    # fieldName Options:
    # /heights/lat_ph
    # /heights/lon_ph
    # /heights/h_ph
    # /heights/delta_time
    # /heights/crossing_time
    # /heights/signal_conf_ph

    # Initialize output
    dataOut = []
    
    if not os.path.isfile(in_file03):
      print('ATL03 file does not exist')
    try:
      with h5py.File(in_file03, 'r') as f:
          dsname=''.join([label, fieldName])
          if dsname in f:
              dataOut = np.array(f[dsname])
              if('signal_conf_ph' in fieldName.lower()):
                  dataOut = dataOut[:,0]
          else:
              dataOut = []
    except Exception as e:
        print('Python message: %s\n' % e)
    return dataOut


##### Functions to read ATL08 .h5 files
def readAtl08H5(in_file08, fieldName, label):
    
    # fieldName Options:
    # /land_segments/longitude
    # /land_segments/latitude
    # /land_segments/canopy/h_max_canopy_abs
    # /land_segments/terrain/h_te_best_fit
    # /land_segments/terrain/h_te_median

    # Initialize output
    dataOut = []

    if not os.path.isfile(in_file08):
      print('ATL08 file does not exist')
    try:
      with h5py.File(in_file08, 'r') as f:
          dsname=''.join([label, fieldName])
          if dsname in f:
              dataOut = np.array(f[dsname])
          else:
              dataOut = []
    except Exception as e:
        print('Python message: %s\n' % e)
    return dataOut

def readAtlH5(in_file, fieldName):
    # Initialize output
    dataOut = []
    
    if not os.path.isfile(in_file):
      print('ATL08 file does not exist')
    try:
      with h5py.File(in_file, 'r') as f:
          dsname=''.join([fieldName])
          if dsname in f:
              dataOut = np.array(f[dsname])
          else:
              dataOut = []
    except Exception as e:
        print('Python message: %s\n' % e)
    return dataOut

##### Function to read ATL03 .h5 files for mapping
def readAtl03DataMapping(in_file03, label):
#
# Reads the data from ATL03
#
# Open the file
#
  if not os.path.isfile(in_file03):
    print('File does not exist')

  try:
    f = h5py.File(in_file03, 'r')
  except Exception as e:
    print('Python message: %s\n' % e)
    return [], []
# endif
#
# segment_ph_cnt
#

#
# segment_id
#
  dsname=label+'/geolocation/segment_id'
  if dsname in f:
    segment_id=np.array(f[dsname])
  else:
    segment_id=[]
# endif
#
# segment_lat
#


#
# ph_index_beg
#
  dsname=label+'/geolocation/ph_index_beg'
  if dsname in f:
    ph_index_beg=np.array(f[dsname])
  else:
    ph_index_beg=[]
# endif
#

#
# Close the file
#
  f.close()

  return ph_index_beg, segment_id 


##### Function to read ATL08 .h5 files for mapping
def readAtl08DataMapping(in_file08, label):
#
# Reads the data from ATL08
#
# Open the file
#
  if not os.path.isfile(in_file08):
    print('File does not exist')
  try:
    f = h5py.File(in_file08, 'r')
  except Exception as e:
    print('Python message: %s\n' % e)
    return [], [], []
# endif
#
# classed_pc_indx
#
  dsname=label+'/signal_photons/classed_pc_indx'
  if dsname in f:
    classed_pc_indx=np.array(f[dsname])
  else:
    classed_pc_indx=[]
# endif
#
# classed_pc_flag
#
  dsname=label+'/signal_photons/classed_pc_flag'
  if dsname in f:
    classed_pc_flag=np.array(f[dsname])
  else:
    classed_pc_flag=[]
# endif
#
# d_flag
#
 
# seg08_id
#
  dsname=label+'/signal_photons/ph_segment_id'
  if dsname in f:
    seg08_id=np.array(f[dsname])
  else:
    seg08_id=[]
# endif
#
#
# Close the file
#
  f.close()

  return classed_pc_indx, classed_pc_flag, seg08_id
    

##### Function to read geoid .mat file
def readGeoidFile(geoidFile):
    
    # Read .mat file
    matFile = os.path.normpath(geoidFile)
    matData = loadmat(matFile)
    
    # Get lats, lons, and geoidal heights
    lats = matData['geoid']['lats'][0][0]
    lons = matData['geoid']['lons'][0][0]
    geoidalHeights = matData['geoid']['geoidalHeight'][0][0]
    
    # Store data as an object
    geoid = geoidStruct(lats, lons, geoidalHeights)
    
    # Return output
    return geoid


##### Function to read .las files
def readLas(lasFilePath, metadata=None):
    
    if metadata is None:
        # Read contents of .las file
        with File(lasFilePath, mode = 'r') as lasFile:
        
            # Store output from .las file
            x = lasFile.x
            y = lasFile.y
            z = lasFile.z
            classification = lasFile.classification
            intensity = lasFile.intensity
            headerData = lasFile.header
            
            # Store output into class structure
            lasData = lasStruct(x, y, z, classification, intensity, headerData)
            
        # EndWith
        
        return lasData

    else:
        # reading some metadata instead
        # need:
        #   pdal installed (sudo yum install pdal)
        #   subprocess, json
        import subprocess, json
        def search_dict(obj, key):
            # https://stackoverflow.com/questions/14962485/finding-a-key-recursively-in-a-dictionary
            # recursively search thru json dict
            if key in obj: return obj[key]
            for k, v in obj.items():
                if isinstance(v,dict):
                    item = search_dict(v, key)
                    if item is not None:
                        return item

        cmd = 'pdal info {} --metadata'.format(lasFilePath)
        output = subprocess.check_output(cmd, shell=True)
        output = output.decode()

        json_info = json.loads(output)
        if metadata != 'json':
            rtn = search_dict(json_info, metadata)
            if rtn is None:
                raise ValueError('%s metadata not found' % metadata)
            else:
                return rtn
        else:
            return json_info

##### Functions to write .las files
def selectwkt(proj,hemi=None,zone=None):
    if proj.lower() == "utm":
        if zone:
            zone = str(zone)
            if hemi.lower() == "n":
                if zone == "1" or zone == "01":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 1N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32601"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "2" or zone == "02":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 2N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-171],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32602"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "3" or zone == "03":
                    wkt = b'''ROJCS["WGS 84 / UTM zone 3N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32603"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "4" or zone == "04":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 4N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-159],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32604"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "5" or zone == "05":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 5N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-153],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32605"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "6" or zone == "06":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 6N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32606"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "7" or zone == "07":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 7N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32607"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "8" or zone == "08":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 8N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-135],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32608"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "9" or zone == "09":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 9N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32609"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "10": 
                    wkt = b'''PROJCS["WGS 84 / UTM zone 10N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-123],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32610"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''                        
                elif zone == "11":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 11N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32611"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "12":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 12N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32612"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "13":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 13N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32613"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "14":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 14N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-99],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32614"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "15":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 15N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32615"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "16":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 16N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-87],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32616"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "17":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 17N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-81],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32617"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "18":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 18N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32618"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "19":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 19N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-69],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32619"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "20":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 20N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-63],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32620"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "21":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 21N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-57],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32621"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "22":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 22N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32622"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "23":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 23N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32623"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "24":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 24N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-39],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32624"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "25":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 25N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-33],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32625"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "26":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 26N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-27],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32626"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "27":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 27N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-21],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32627"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "28":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 28N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32628"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "29":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 29N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32629"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "30":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 30N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32630"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "31":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 31N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32631"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "32":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 32N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32632"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "33":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 33N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32633"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "34":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 34N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",21],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32634"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "35":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 35N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",27],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32635"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "36":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 36N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",33],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32636"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "37":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 37N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",39],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32637"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "38":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 38N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",45],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32638"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "39":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 39N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32639"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "40":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 40N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",57],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32640"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "41":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 41N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",63],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32641"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "42":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 42N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",69],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32642"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "43":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 43N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32643"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "44":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 44N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",81],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32644"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "45":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 45N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",87],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32645"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "46":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 46N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32646"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "47":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 46N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32646"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "48":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 48N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32648"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "49":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 49N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32649"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "50":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 50N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32650"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "51":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 51N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",123],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32651"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "52":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 52N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32652"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "53":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 53N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",135],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32653"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "54":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 54N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32654"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "55":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 55N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32655"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "56":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 56N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",153],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32656"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "57":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 57N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",159],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32657"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "58":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 58N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32658"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "59":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 59N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",171],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32659"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "60":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 60N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32660"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
            elif hemi.lower() == "s":
                if zone == "1" or zone == "01":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 1S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32701"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "2" or zone == "02":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 2S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-171],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32702"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "3" or zone == "03":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 3S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32703"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "4" or zone == "04":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 4S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-159],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32704"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "5" or zone == "05":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 5S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-153],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32705"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "6" or zone == "06":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 6S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32706"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "7" or zone == "07":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 7S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32707"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "8" or zone == "08":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 8S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-135],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32708"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "9" or zone == "09":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 9S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32709"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "10":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 10S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-123],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32710"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "11":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 11S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32711"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "12":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 12S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32712"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "13":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 13S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32713"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "14":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 14S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-99],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32714"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "15":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 15S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32715"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "16":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 16S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-87],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32716"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "17":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 17S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-81],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32717"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "18":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 18S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32718"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "19":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 19S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-69],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32719"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "20":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 20S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-63],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32720"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "21":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 21S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-57],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32721"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "22":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 22S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32722"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "23":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 23S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32723"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "24":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 24S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-39],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32724"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "25":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 25S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-33],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32725"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "26":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 26S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-27],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32726"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "27":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 27S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-21],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32727"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "28":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 28S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32728"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "29":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 29S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32729"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "30":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 30S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32730"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "31":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 31S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32731"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "32":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 32S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32732"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "33":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 33S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32733"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "34":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 34S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",21],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32734"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "35":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 35S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",27],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32735"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "36":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 36S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",33],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32736"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "37":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 37S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",39],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32737"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "38":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 38S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",45],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32738"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "39":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 39S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32739"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "40":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 40S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",57],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32740"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "41":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 41S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",63],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32741"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "42":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 42S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",69],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32742"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "43":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 43S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32743"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "44":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 44S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",81],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32744"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "45":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 45S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",87],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32745"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "46":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 46S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32746"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "47":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 47S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",99],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32747"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "48":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 48S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32748"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "49":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 49S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32749"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "50":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 50S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32750"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "51":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 51S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",123],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32751"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "52":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 52S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32752"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "53":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 53S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",135],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32753"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "54":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 54S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32754"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "55":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 55S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32755"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "56":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 56S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",153],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32756"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "57":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 57S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",159],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32757"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "58":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 58S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32758"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "59":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 59S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",171],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32759"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "60":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 60S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32760"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
    elif proj.lower() == "arctic" or proj.lower() == "north polar":
        wkt = b'''PROJCS["NSIDC Sea Ice Polar Stereographic North",GEOGCS["Unspecified datum based upon the Hughes 1980 ellipsoid",DATUM["Not_specified_based_on_Hughes_1980_ellipsoid",SPHEROID["Hughes 1980",6378273,298.279411123064,AUTHORITY["EPSG","7058"]],AUTHORITY["EPSG","6054"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4054"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],AUTHORITY["EPSG","3411"],AXIS["X",UNKNOWN],AXIS["Y",UNKNOWN]]'''
    elif proj.lower() == "antarctic" or proj.lower() == "south polar": 
        wkt = b'''PROJCS["NSIDC Sea Ice Polar Stereographic South",GEOGCS["Unspecified datum based upon the Hughes 1980 ellipsoid",DATUM["Not_specified_based_on_Hughes_1980_ellipsoid",SPHEROID["Hughes 1980",6378273,298.279411123064,AUTHORITY["EPSG","7058"]],AUTHORITY["EPSG","6054"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4054"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-70],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],AUTHORITY["EPSG","3412"],AXIS["X",UNKNOWN],AXIS["Y",UNKNOWN]]'''
    elif proj.lower() == "wgs84":
        wkt = b'''GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]''' 
    else:
        print("No defined Projected Coordinate System Selected")
    return wkt


def writeLas(xx,yy,zz,proj,output_file,classification,intensity,signalConf=None,hemi=None,zone=None):
    
    import laspy
    # from laspy.file import File

    wkt = selectwkt(proj,hemi,zone)
    
    #Create new VLR
    new_vlr = laspy.header.VLR(user_id = "LASF_Projection",
                           record_id = 2112,
                           VLR_body = wkt,
                           description = "OGC Coordinate System WKT")
    inVLRs = []
    inVLRs.append(new_vlr)

    #Create new Header
    hdr = laspy.header.Header(file_version=1.4)
    hdr.file_sig = 'LASF'
    
    #Create new las file with Header and VLR  
    outfile = laspy.file.File(output_file, mode="w", header=hdr)
    outfile.header.vlrs = inVLRs
    outfile.header.set_wkt = 1
    
    #Establish offset
    xmin = np.min(xx)
    ymin = np.min(yy)
    zmin = np.min(zz)
    
    xmax = np.max(xx)
    ymax = np.max(yy)
    zmax = np.max(zz)
    
    #Calculate x, y, and z scale factor
    if xmax == xmin:
        xscale = 1
    else:
        xscale = (xmax - xmin) / 100000000;    
    if ymax == ymin:
        yscale = 1
    else:
        yscale = (ymax - ymin) / 100000000;
        
    if zmax == zmin:
        zscale = 1
    else:
        zscale = (zmax - zmin) / 100000000;
    
    in_scale = [xscale, yscale, zscale]
    
    outfile.header.offset = [xmin,ymin,zmin]
    
    #Establish scale
    outfile.header.scale = in_scale
    
    #Write x, y, z data and if available classification and intensity data
    outfile.x = xx
    outfile.y = yy
    outfile.z = zz
    if classification is not None:
        outfile.raw_classification = classification
    if intensity is not None:
        outfile.intensity = intensity
    if signalConf is not None:
        outfile.scan_angle_rank = signalConf
    
    #Close las
    outfile.close()
    
    
##### Functions to write .kml files
def writeKml(lat, lon, time, kmlName):

    import simplekml
    # Suppress warnings that may come from simple kml
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    # endif
    
    # Open Simple KML
    kml = simplekml.Kml()
    
    # Plot line for ground track
    latLon = np.column_stack((lon, lat))
    latLonTuple = tuple(latLon)
    ls = kml.newlinestring(name = 'test', coords=latLonTuple)
    ls.extrude = 1
    ls.altitudemode = simplekml.AltitudeMode.clamptoground
    ls.style.linestyle.width = 5
    ls.style.linestyle.color = simplekml.Color.blue
        
    # Open Simple KML style editor
    style = simplekml.Style()
    style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'

    # Loop through all lat/lon values and make KML markers
    for i in range(0,len(lon)):
        
        # Get time string
        timeRounded = str(np.round(time[i],1)) + ' sec'
        
        # Plot marker points
        pnt = kml.newpoint(name = timeRounded, coords=[(lon[i][0], lat[i][0])])
        pnt.style = style
        
    # EndFor

    # Save output KML file
    kml.save(kmlName)


def writeArrayToCSV(csv_out,namelist,datalist,header=True):
    if datalist:
        in_data = datalist[0]
        if len(datalist) > 0:
            for data in datalist[1:]:
               in_data = np.column_stack((in_data,data))
    with open(csv_out, 'w', newline = '') as csvFile:
        writer = csv.writer(csvFile)
        if(header):
            writer.writerow(namelist)
        # endIf
        writer.writerows(in_data)
    
    csvFile.close

def writeATL08toCSV(in_file08,groundtrack,csv_out):
    delta_time = readAtl08H5(in_file08, '/land_segments/delta_time', 
                                 groundtrack)
    lat = readAtl08H5(in_file08, '/land_segments/latitude', groundtrack)
    lon = readAtl08H5(in_file08, '/land_segments/longitude', groundtrack)
    h_max_canopy_abs = readAtl08H5(in_file08, 
                                   '/land_segments/canopy/h_max_canopy_abs', 
                                   groundtrack)
    h_te_best_fit = readAtl08H5(in_file08, 
                                '/land_segments/terrain/h_te_best_fit', 
                                groundtrack)
    h_te_median = readAtl08H5(in_file08, '/land_segments/terrain/h_te_median', 
                              groundtrack)
    namelist = ['Delta Time','Latitude','Longitude',
                'Absolute Max Canopy Height','Best Fit Ground Height',
                'Median Ground Height']
    datalist = [delta_time,lat,lon,h_max_canopy_abs,h_te_best_fit,h_te_median]
    writeArrayToCSV(csv_out,namelist,datalist)
    
    
def writeATL03toCSV(in_file03,groundtrack,csv_out):
    delta_time = readAtl03H5(in_file03, '/heights/delta_time', groundtrack)
    lat = readAtl03H5(in_file03, '/heights/lat_ph', groundtrack)
    lon = readAtl03H5(in_file03, '/heights/lon_ph', groundtrack)
    h_ph = readAtl03H5(in_file03, '/heights/h_ph', groundtrack)
    signal_conf_ph = readAtl03H5(in_file03, '/heights/signal_conf_ph', groundtrack)
    namelist = ['Delta Time','Latitude','Longitude','Height',
                'Signal Confidence']
    datalist = [delta_time,lat,lon,h_ph,signal_conf_ph]
    writeArrayToCSV(csv_out,namelist,datalist)
    

def isconnected():
    try:
        host = socket.gethostbyname("https://www.nasa.gov/")
        s = socket.create_connection((host, 80), 2)
        s.close()
        print("Internet Connection Established")
        return True
    except:
        print("No Internet Connection Established")
        return False
 
    
def createHTMLChart(ytrack, h_ph, classification, lat, lon, 
                    direction = 'Descending',
                    online = None,
                    classification_list = [1,2,3],
                    input_folder = "",
                    output_folder = "",
                    in_file03_name = "ATL03", 
                    blank = viewerBlank_html, 
                    blanki = viewerBlankOnline_html):
    
    if len(ytrack.shape) == 2:
        ytrack = ytrack[:,0]
    
    if len(h_ph.shape) == 2:
        h_ph = h_ph[:,0]

    if len(classification.shape) == 2:
        classification = classification[:,0]
        
    if len(lat.shape) == 2:
        lat = lat[:,0]
    
    if len(lon.shape) == 2:
        lon = lon[:,0]
    
    if online == None:
        online = isconnected()
        
    if online == True:
        viewer_output = (output_folder + "\Viewer_Online" + in_file03_name + ".html")
        blank = blanki
    else:
        viewer_output = (output_folder + "\Viewer_" + in_file03_name + ".html")
        
        # Copy d3 file to output directory
        srcPath = os.path.normpath(input_folder + '\\d3')
        tarPath = os.path.normpath(output_folder + '\\d3')
        shutil.copytree(srcPath, tarPath)

    #Copy the Blank Template into 
    with open(blank) as f:
        lines = f.readlines()
        lines = [l for l in lines]
        with open(viewer_output, "w") as f1:
            f1.writelines(lines)

    
    #Read file and load it to memory
    with open(viewer_output, "r") as in_file:
        buffer = in_file.readlines()

    #Write data into the HTML file
    with open(viewer_output, "w") as out_file:
        for line in buffer:
            if line == "var data = [\n":
                for j in range(0,len(classification)):
                    if online == True:
                        if ((classification[j] == 0) &
                            (0 in classification_list)):
                            line = line + ("{ytrack: " 
                                             + str(ytrack[j]) + ", zheight: " 
                                             + str(h_ph[j]) +
                                             ", color: '#C2C5CC', lat: " 
                                             + str(lat[j]) + ", lon: " 
                                             + str(lon[j]) + "}, \n") 
                        elif ((classification[j] == 1) &
                              (1 in classification_list)):
                            line = line + ("{ytrack: "
                                             + str(ytrack[j]) + ", zheight: "
                                             + str(h_ph[j]) +
                                             ", color: '#D2B826', lat: " 
                                             + str(lat[j]) + ", lon: " 
                                             + str(lon[j]) + "}, \n") 
                        elif ((classification[j] == 2) &
                              (2 in classification_list)):
                            line = line + ("{ytrack: " 
                                             + str(ytrack[j])  
                                             + ", zheight: " + str(h_ph[j]) 
                                             + ", color: '#45811A', lat: " 
                                             + str(lat[j]) + ", lon: " 
                                             + str(lon[j]) + "}, \n") 
                        elif ((classification[j] == 3) & 
                              (3 in classification_list)):
                            line = line + ("{ytrack: " 
                                             + str(ytrack[j]) + ", zheight: " 
                                             + str(h_ph[j]) 
                                             + ", color: '#85F334', lat: " 
                                             + str(lat[j]) + ", lon: " 
                                             + str(lon[j]) + "}, \n") 
                            j += 1
                    else:
                        if ((classification[j] == 0) &
                            (0 in classification_list)):
                            line = line + ("{ytrack: " 
                                             + str(ytrack[j]) + ", zheight: " 
                                             + str(h_ph[j]) +
                                             ", color: '#C2C5CC' }, \n") 
                        elif ((classification[j] == 1) &
                              (1 in classification_list)):
                            line = line + ("{ytrack: " 
                                             + str(ytrack[j]) + ", zheight: " 
                                             + str(h_ph[j]) +
                                             ", color: '#D2B826' }, \n")
                        elif ((classification[j] == 2) &
                              (2 in classification_list)):
                            line = line + ("{ytrack: " 
                                             + str(ytrack[j])  
                                             + ", zheight: " + str(h_ph[j]) 
                                             + ", color: '#45811A' }, \n")
                        elif ((classification[j] == 3) & 
                              (3 in classification_list)):
                            line = line + ("{ytrack: " 
                                             + str(ytrack[j]) + ", zheight: " 
                                             + str(h_ph[j]) 
                                             + ", color: '#85F334' }, \n")
                            j += 1
            if line == "//INCLUDE BASELINE\n":
                if direction == 'Ascending':
                    line = line + "baseline = L.polyline([[lat1, lon0],\n" + \
                    "[lat0,  lon1]],{opacity:0.5, color:'gray'})" + \
                    ".addTo(mymap).bindPopup('Icesat-2 Track' );\n"
                else:
                    line = line + "baseline = L.polyline([[lat1, lon1],\n" + \
                    "[lat0,  lon0]],{opacity:0.5, color:'gray'})" + \
                    ".addTo(mymap).bindPopup('Icesat-2 Track' );\n"     
                    
            if line == "//INCLUDE GUIDELINE\n":
                if direction == 'Ascending':
                    line = line + "guideline = L.polyline([[lat1, lon0],\n" + \
                    "[lat0,  lon1]],{opacity:1, color:'red'})" + \
                    ".addTo(mymap).bindPopup('Icesat-2 Track' );\n"
                else:
                    line = line + "guideline = L.polyline([[lat1, lon1],\n" + \
                    "[lat0,  lon0]],{opacity:1, color:'red'})" + \
                    ".addTo(mymap).bindPopup('Icesat-2 Track' );\n"                     
                    
                
                
            out_file.write(line)
        
    return viewer_output
                 
            
def getDEMArrays(file):
    # Open DEM
    ds = gdal.Open(file)
    
    # Read DEM as NP Array
    data = np.array(ds.GetRasterBand(1).ReadAsArray())
    
    # Get Geotransform Information
    gt = ds.GetGeoTransform()
    
    #Generate X Array
    row = np.arange(data.shape[1])
    for pos in range(0,data.shape[1]):
        row[pos] = gt[0] + (gt[1] * row[pos])
    xarr = np.array([row,]*data.shape[0])
    
    #Generate Y Array
    column = np.arange(data.shape[0])  
    for pos in range(0,data.shape[0]):
        column[pos] = gt[3] + (gt[5] * column[pos])  
    yarr = np.array([column,]*data.shape[1]).transpose()
    
    # Return Data (zarray), xarray, and zarray
    return data, xarr, yarr

def readDEMepsg(file):
    ds = gdal.Open(file)
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    epsg = proj.GetAttrValue('AUTHORITY',1)
    return epsg

def formatDEM(file, epsg_backup = 0):
    data, xarr, yarr = getDEMArrays(file)
    
    epsg = readDEMepsg(file)
    
    if epsg == None:
        if epsg_backup == 0:
            print('Missing Projection Information')
        else:
            epsg = epsg_backup
            
    
    data = data.flatten()
    xarr = xarr.flatten()
    yarr = yarr.flatten()
    
    xarr = xarr[data > -999]
    yarr = yarr[data > -999]
    zarr = data[data > -999]
    intensity = np.ones(len(zarr))
    classification = np.ones(len(zarr)) * 2
    
    return xarr, yarr, zarr, intensity, classification, epsg

def write_geotiff(data, epsg, x_min, y_max, x_pixel, y_pixel, outputfile, 
                    nodata = None, write_raster = False):

    x_offset = 0.0
    y_offset = 0.0
    if write_raster:
      x_offset = -x_pixel/2.0
      y_offset = y_pixel/2.0

    x_pixels = data.shape[1]
    y_pixels = data.shape[0]
    
    driver = gdal.GetDriverByName('GTiff')

    dataset = driver.Create(
        outputfile,
        x_pixels,
        y_pixels,
        1,
        gdal.GDT_Float64, )

    dataset.SetGeoTransform((
        x_min + x_offset,    
        x_pixel,  
        0,                      
        y_max + y_offset,    
        0,                      
        -y_pixel))  

    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    
    dataset.SetProjection(srs.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(data)
    if nodata:
        dataset.GetRasterBand(1).SetNoDataValue(nodata)
    dataset.FlushCache()


def write_mat(file_mat, data, datasets, debug=0):

  """
  Function to write .mat files from 1D numpy arrays
  Input:
    file_mat - full path of output .mat file
    data - list of numpy arrays of data
      i.e. [x, y, z], where x, y, and z are numpy arrays, not necessarily of same size
      x, y, z must be (1D arrays) of dimensions (n,1) or (n,) where n is number of elements
    datasets - list of descriptions for each data
      i.e. for x, y, and z, we have datasets = ['along', 'across', 'height']
      Note that datasets must contain variable names allowed in MATLAB

    Note that data types are maintained.

  Output:
    .mat file with file_mat name

  Example:
    # load in atl03Data, import icesatIO functions
    OUT_DIR = '.'
    file_mat = OUT_DIR + '/test.mat'
    datasets = ['class','signal_conf','along_track','z']
    data = [atl03Data.classification, atl03Data.signalConf, atl03Data.alongTrack, atl03Data.z]
    write_mat(file_mat, data, datasets, debug=1)

  """

  from scipy.io import savemat

  def valid_variable_name(name):

    import string
    allowed_char = list(string.ascii_lowercase) + \
                    list(string.ascii_uppercase) + \
                    list(string.digits) + ['_']

    matlab_keywords = [
    'break',
    'case',
    'catch',
    'classdef',
    'continue',
    'else',
    'elseif',
    'end',
    'for',
    'function',
    'global',
    'if',
    'otherwise',
    'parfor',
    'persistent',
    'return',
    'spmd',
    'switch',
    'try',
    'while']

    # must start with a letter or underscore
    # cannot start with a number
    name_lower = name.lower()
    if name_lower[0].islower() or name[0] == '_':
      # cannot be a matlab keyword
      if not (name in matlab_keywords):
        # must contain only ascii letters, digits, and underscores
        c = 0
        for char in name:
          if char in allowed_char:
            c += 1

        if c == len(name):
          return True

    return False

  err = False
  for name in datasets:
    if not valid_variable_name(name):
      err = True
      print('error: \"%s\" is an invalid MATLAB variable name' % name)
      break

  n_datasets = len(datasets)
  n_data = len(data)
  if n_datasets != n_data and not err:
    print('error: n_datasets != n_data')
    err = True

  if n_datasets == 0 and not err:
    print('error: len(datasets) == 0')
    err = True

  if type(datasets) != list or type(data) != list and not err:
    print('error: type(datasets) != list or type(data) != list')
    err = True

  valid = 0
  if not err:
    data_out = {}
    
    for k in range(n_datasets):
      d = data[k]
      
      if type(d) == list:
        u = np.unique(['%s' % type(val) for val in d])
        if len(u) == 1:
          d = np.array(d).astype(type(d[0]))
        else:
          print('error: inconsistent datatypes in \"%s\"' % datasets[k])
          break
        # endIf
      # endIf
      
      if(not(isinstance(d,str))):
          if len(d.shape) == 1:
            d = np.c_[d]
          elif len(d.shape) > 1:
            if not (d.shape[1] == 1):
              print('error: \"%s\" must be of (n>0,1) or (n>0,) dimensions' % datasets[k])
              break
          elif len(d.shape) == 0:
            print('error: \"%s\" must be of (n>0,1) or (n>0,) dimensions' % datasets[k])
            break

          if d.shape[0] < 1:
            print('error: \"%s\" must be of (n>0,1) or (n>0,) dimensions' % datasets[k])
            break
      # endIf
      
      data_out[datasets[k]] = data[k]
      valid += 1
    #endFor
    
  if valid == n_datasets:
    savemat(file_mat, data_out)
    if debug:
      print('file written:')
      print(file_mat)
  else:
    if debug:
      print('file not written:')
      print(file_mat)


def calculateangle(x1,x2,y1,y2):
    if (x2 - x1) == 0:
        slope = np.inf
    else:
        slope = (y2 - y1)/(x2 - x1)
    degree = np.rad2deg(np.arctan(slope))
    return degree

def calculategrounddirection(xx,yy):
    degree = np.zeros(len(xx))
    for i in range(0,len(xx)):
        if i == 0:
            degree[i] = calculateangle(xx[i], xx[i+1], yy[i], yy[i+1])
        elif i == (len(xx))-1:
            degree[i]  = calculateangle(xx[i-1], xx[i], yy[i-1], yy[i])
        else:
            degree[i]  = calculateangle(xx[i-1], xx[i+1], yy[i-1], yy[i+1])
    return degree
        
def rotatepoint(degree,xpos,ypos):
    angle = np.deg2rad(degree)
    xrot = (xpos * np.cos(angle)) - (ypos * np.sin(angle)) 
    yrot = (xpos * np.sin(angle)) + (ypos * np.cos(angle))
    return xrot, yrot

def calculatecorners(degree,xcenter,ycenter,width,height):
    # Set corner values
    xul = -width / 2
    yul = height / 2
    xur = width / 2
    yur = height / 2
    xll = -width / 2
    yll = -height / 2
    xlr = width / 2
    ylr = -height / 2
    
    # Rotate based on the angle degree
    xul, yul = rotatepoint((degree-90),xul,yul)
    xur, yur = rotatepoint((degree-90),xur,yur)
    xll, yll = rotatepoint((degree-90),xll,yll)
    xlr, ylr = rotatepoint((degree-90),xlr,ylr)
    
    # Add corner values to centeroid
    xul = xcenter + xul
    yul = ycenter + yul
    xur = xcenter + xur
    yur = ycenter + yur
    xll = xcenter + xll
    yll = ycenter + yll
    xlr = xcenter + xlr
    ylr = ycenter + ylr
    
    return xul, yul, xur, yur, xll, yll, xlr, ylr

def createShapefiles(xx, yy, width, height, epsg, outfile = "atl08.shp"):
    # Generate list of degrees
    degreelist = calculategrounddirection(xx,yy)
    
    # Define Esri Shapefile output
    driver = ogr.GetDriverByName('Esri Shapefile')
    
    # Name output shape file (foo.shp)
    ds = driver.CreateDataSource(outfile)
    
    # Define spatial reference based on EPSG code 
    # https://spatialreference.org/ref/epsg/
    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    
    # Create file with srs
    layer = ds.CreateLayer('', srs, ogr.wkbPolygon)
    
    # Create arbitary id field
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    defn = layer.GetLayerDefn()
    
    # Create a new feature (attribute and geometry)
    for i in range(0,len(xx)):
        # Generate the corner points
        xul, yul, xur, yur, xll, yll, xlr, ylr  = \
        calculatecorners(degreelist[i],xx[i],yy[i],width,height)     
        
        # Create rectangle corners
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(xul, yul)
        ring.AddPoint(xur, yur)
        ring.AddPoint(xlr, ylr)
        ring.AddPoint(xll, yll)
        ring.AddPoint(xul, yul)
        
        # Create polygon from corners
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        
        # Export well-known binary
        wkb = poly.ExportToWkb()
        
        # Assign arbitary number to field ID
        feat = ogr.Feature(defn)
        feat.SetField('id', i)
                
        # Make a geometry, from Shapely object
        geom = ogr.CreateGeometryFromWkb(wkb)
        feat.SetGeometry(geom)
        
        # Write out geometry
        layer.CreateFeature(feat)
        
        # Remove ring and poly
        ring = poly = None
    
    # Remove feat and geom
    feat = geom = None
    
    # Save and close everything
    ds = layer = feat = geom = None    
    
def read_geotiff(file, read_raster=False):
    # Open DEM
    ds = gdal.Open(file)
    
    # Read DEM as NP Array
    data = np.array(ds.GetRasterBand(1).ReadAsArray())
    
    # Get Geotransform Information
    gt = ds.GetGeoTransform()
    
    # Get upper-left coordinates and pixel resolution
    ulx = gt[0]
    resx = gt[1]
    uly = gt[3]
    resy = gt[5] # usually negative

    x_offset = 0.0
    y_offset = 0.0
    if read_raster:
      # opposite of write_geotiff
      x_offset = resx/2.0
      y_offset = -abs(resy)/2.0
      # y_offset = resy/2.0

    ulx += x_offset
    uly += y_offset
    
    # Check GT[2] and GT[4] are 0, if not, throw warning
    if gt[2] != 0 or gt[4] != 0:
        print('WARNING: CASE NOT ACCOUNTED')
    
    # Check if resolution for X and Y are the same
    # if np.abs(resx) != np.abs(resy):
    #     print('WARNING: Unequal X-Y Pixel Size')
    
    # Get EPSG Code
    proj = osr.SpatialReference(wkt=ds.GetProjection())
    epsg = proj.GetAttrValue('Authority',1)
    
    # Return Data (zarray), xarray, and zarray
    return data, epsg, ulx, uly, resx, resy


def find_intersecting_values(x,y,data,ulx,uly,resx,resy):
    # Calcualte opposite extents
    if not isinstance(x,np.ndarray):
        x = np.array([x])
        y = np.array([y])
    lly = uly + (resy * data.shape[0])
    urx = ulx + (resx * data.shape[1])
    # Create out of range (oar) filter
    oar = np.where((x > urx) | (x < ulx) | (y > uly) | (y < lly))    
    # Calcualte alpha-beta coordiantes
    alpha = np.floor((y/resy) - (uly/resy))
    beta = np.floor((x/resx) - (ulx/resx))
    # If there are out of bound coordinates, assign them to 0 for now    
    if oar[0].size > 0:    
        alpha[oar] = 0
        beta[oar] = 0    
    # Create results with alpha-beta coordinates of data
    result = np.array(
        [data[int(alpha[i]),int(beta[i])] for i in range(0,len(alpha))])
    # If there are out of range values, assign nan
    if oar[0].size > 0:    
        result[oar] = np.nan
    return result

if __name__ == "__main__":
    print("Test")
    
    
### Function to write console and .log messages
def writeLog(message, fileID=False):
    
    # Print message to screen
    print(message)
        
    # Print message to .log file
    if(fileID):   
        print(message, file=fileID)
    # endIf
    
# endDef