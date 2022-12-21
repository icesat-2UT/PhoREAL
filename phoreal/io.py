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

# Filter Runtime warnings
import warnings
warnings.filterwarnings('ignore')

# Import Python modules
import os
import sys
import csv
import numpy as np
import h5py
import ntpath
import glob
from scipy.io import loadmat
import socket
import shutil
import pandas as pd
import rasterio
import laspy
from laspy.file import File
import simplekml

# Import ICESat-2 modules
from phoreal.gui_addins import (viewerBlank_html, viewerBlankOnline_html)
from phoreal.utils import (identifyEPSG, getCoordRotFwd, transform, getGeoidHeight, \
                         getCoordRotRev, superFilter, getRaster, getUTM2LatLon)


# Object for readKmlBounds function
class kmlStruct:
        
    # Define class with designated fields
    def __init__(self, regionName, lonMin, lonMax, latMin, latMax):
            
        self.regionName = regionName
        self.lonMin = lonMin
        self.lonMax = lonMax
        self.latMin = latMin
        self.latMax = latMax
        
    # endDef
# endClass

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
    # endDef
# endClass
    
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
    # endDef
# endClass
    
# Object for readGeoidFile function
class geoidStruct:
        
    # Define class with designated fields
    def __init__(self, lats, lons, geoidalHeights):
            
        self.lats = lats
        self.lons = lons
        self.geoidalHeights = geoidalHeights
    # endDef
# endClass
        
class atl03Struct:
        
    # Define class with designated fields
    def __init__(self, atl03_lat, atl03_lon, atl03_easting, atl03_northing, 
                 atl03_crossTrack, atl03_alongTrack, atl03_z, atl03_zMsl, 
                 atl03_time, atl03_deltaTime,
                 atl03_signalConf, 
                 atl03_yapcConf, atl03_yapcSnr, atl03_yapcSnrNorm,
                 atl03_refDem,
                 atl03_classification, atl03_intensity,
                 atl03_solar_elev,
                 atl03_segment_id,
                 gtNum, beamNum, beamStrength, zone, hemi,
                 atl03FilePath, atl03FileName, trackDirection, alt03h5Info, dataIsMapped):
            
        self.lat = np.c_[atl03_lat]
        self.lon = np.c_[atl03_lon]
        self.easting = np.c_[atl03_easting]
        self.northing = np.c_[atl03_northing]
        self.crossTrack = np.c_[atl03_crossTrack]
        self.alongTrack = np.c_[atl03_alongTrack]
        self.z = np.c_[atl03_z]
        self.zMsl = np.c_[atl03_zMsl]
        self.time = np.c_[atl03_time]
        self.deltaTime = np.c_[atl03_deltaTime]
        self.signalConf = np.c_[atl03_signalConf]
        self.yapcConf = np.c_[atl03_yapcConf]
        self.yapcSnr = np.c_[atl03_yapcSnr]
        self.yapcSnrNorm = np.c_[atl03_yapcSnrNorm]
        self.refDem = np.c_[atl03_refDem]
        self.classification = np.c_[atl03_classification]
        self.intensity = np.c_[atl03_intensity]
        self.solar_elev = np.c_[atl03_solar_elev]
        self.segmentID = np.c_[atl03_segment_id]
        self.gtNum = gtNum
        self.beamNum = beamNum
        self.beamStrength = beamStrength
        self.zone = zone
        self.hemi = hemi
        self.atl03FilePath = atl03FilePath
        self.atl03FileName = atl03FileName
        self.trackDirection = trackDirection
        self.atlVersion = alt03h5Info.atlVersion
        self.year = alt03h5Info.year
        self.month = alt03h5Info.month
        self.day = alt03h5Info.day
        self.hour = alt03h5Info.hour
        self.minute = alt03h5Info.minute
        self.second = alt03h5Info.second
        self.trackNum = alt03h5Info.trackNum
        self.unknown = alt03h5Info.unknown
        self.releaseNum = alt03h5Info.releaseNum
        self.incrementNum = alt03h5Info.incrementNum
        self.dataIsMapped = dataIsMapped
    # endDef
# endClass 

class atl08Struct:
        
    # Define class with designated fields
    def __init__(self, atl08_lat, atl08_lon, atl08_easting, atl08_northing, 
                 atl08_crossTrack, atl08_alongTrack, 
                 atl08_maxCanopy, atl08_teBestFit, atl08_teMedian, 
                 atl08_maxCanopyMsl, atl08_teBestFitMsl, atl08_teMedianMsl,
                 atl08_time, atl08_deltaTime,
                 atl08_signalConf, atl08_classification, atl08_intensity, 
                 gtNum, beamNum, beamStrength, zone, hemi, 
                 atl08FilePath, atl08FileName, trackDirection, alt08h5Info, dataIsMapped):
            
        self.lat = np.c_[atl08_lat]
        self.lon = np.c_[atl08_lon]
        self.easting = np.c_[atl08_easting]
        self.northing = np.c_[atl08_northing]
        self.crossTrack = np.c_[atl08_crossTrack]
        self.alongTrack = np.c_[atl08_alongTrack]
        self.maxCanopy = np.c_[atl08_maxCanopy]
        self.teBestFit = np.c_[atl08_teBestFit]
        self.teMedian = np.c_[atl08_teMedian]
        self.maxCanopyMsl = np.c_[atl08_maxCanopyMsl]
        self.teBestFitMsl = np.c_[atl08_teBestFitMsl]
        self.teMedianMsl = np.c_[atl08_teMedianMsl]
        self.time = np.c_[atl08_time]
        self.deltaTime = np.c_[atl08_deltaTime]
        self.signalConf = np.c_[atl08_signalConf]
        self.classification = np.c_[atl08_classification]
        self.intensity = np.c_[atl08_intensity]
        self.gtNum = gtNum
        self.beamNum = beamNum
        self.beamStrength = beamStrength
        self.zone = zone
        self.hemi = hemi
        self.atl08FilePath = atl08FilePath
        self.atl08FileName = atl08FileName
        self.trackDirection = trackDirection
        self.atlVersion = alt08h5Info.atlVersion
        self.year = alt08h5Info.year
        self.month = alt08h5Info.month
        self.day = alt08h5Info.day
        self.hour = alt08h5Info.hour
        self.minute = alt08h5Info.minute
        self.second = alt08h5Info.second
        self.trackNum = alt08h5Info.trackNum
        self.unknown = alt08h5Info.unknown
        self.releaseNum = alt08h5Info.releaseNum
        self.incrementNum = alt08h5Info.incrementNum
        self.dataIsMapped = dataIsMapped
    # endDef
# endClass        
        
class atlRotationStruct:
    
    # Define class with designated fields
    def __init__(self, R_mat, xRotPt, yRotPt, desiredAngle, phi):
        
        self.R_mat = R_mat
        self.xRotPt = xRotPt
        self.yRotPt = yRotPt
        self.desiredAngle = desiredAngle
        self.phi = phi
    # endDef
# endClass
    
# Truth class 
class atlTruthStruct:
        
    # Define class with designated fields
    def __init__(self, easting, northing, crossTrack, alongTrack, lat, lon, z, 
                 classification, intensity, year, month, day, zone, hemi, epsg=False):
            
        self.easting = np.c_[easting]
        self.northing = np.c_[northing]
        self.crossTrack = np.c_[crossTrack]
        self.alongTrack = np.c_[alongTrack]
        self.lat = np.c_[lat]
        self.lon = np.c_[lon]
        self.z = np.c_[z]
        self.classification = np.c_[classification]
        self.intensity = np.c_[intensity]
        self.time = np.zeros(np.shape(self.easting))
        self.deltaTime = np.zeros(np.shape(self.easting))
        self.year = year
        self.month = month
        self.day = day
        self.zone = zone
        self.hemi = hemi
        self.epsg = epsg


    # endDef
    
    # Define append method
    def append(self, newClass):
        
        self.easting = np.concatenate((self.easting, newClass.easting), axis=0)
        self.northing = np.concatenate((self.northing, newClass.northing), axis=0)
        self.crossTrack = np.concatenate((self.crossTrack, newClass.crossTrack), axis=0)
        self.alongTrack = np.concatenate((self.alongTrack, newClass.alongTrack), axis=0)
        self.lat = np.concatenate((self.lat, newClass.lat), axis=0)
        self.lon = np.concatenate((self.lon, newClass.lon), axis=0)
        self.z = np.concatenate((self.z, newClass.z), axis=0)
        self.classification = np.concatenate((self.classification, newClass.classification), axis=0)
        self.intensity = np.concatenate((self.intensity, newClass.intensity), axis=0)
        self.time = np.concatenate((self.time, newClass.time), axis=0)
        self.deltaTime = np.concatenate((self.deltaTime, newClass.deltaTime), axis=0)
        self.year = np.concatenate((self.year, newClass.year), axis=0)
        self.month = np.concatenate((self.month, newClass.month), axis=0)
        self.day = np.concatenate((self.day, newClass.day), axis=0)
        self.zone = newClass.zone
        self.hemi = newClass.hemi
        self.epsg = newClass.epsg
        
    # endDef
# endClass
        
class offsetsStruct:
    
    # Define class with designated fields
    def __init__(self, crossTrackBounds, alongTrackBounds, rasterResolutions, 
                 useVerticalShift, verticalShift):
        
        self.crossTrackBounds = crossTrackBounds
        self.alongTrackBounds = alongTrackBounds
        self.rasterResolutions = rasterResolutions
        self.useVerticalShift = useVerticalShift
        self.verticalShift = verticalShift
    # endDef
# endClass

class atlMeasuredDataReducedStruct:
    
    # Define class with designated fields
    def __init__(self, easting, northing, z, crossTrack, alongTrack, lat, lon, 
                 time, classification, signalConf):
        
        self.easting = easting
        self.northing = northing
        self.z = z
        self.crossTrack = crossTrack
        self.alongTrack = alongTrack
        self.lat = lat
        self.lon = lon
        self.time = time
        self.classification = classification
        self.signalConf = signalConf
    # endDef
# endClass

class atlCorrectionsStruct:
    
    # Define class with designated fields
    def __init__(self, shiftedEasting, shiftedNorthing,
                 shiftedCrossTrack, shiftedAlongTrack,
                 shiftedLat, shiftedLon,
                 shiftedVertical, shiftedVerticalMsl,
                 time, deltaTime, classification,
                 signalConf, zone, hemi,
                 correctionsCrossTrack, correctionsAlongTrack, 
                 correctionsEasting, correctionsNorthing, 
                 correctionsVertical, 
                 correctionsMAE, correctionsRMSE, correctionsME, 
                 measRasterXCommonFinal, measRasterYCommonFinal, measRasterZCommonFinal,
                 truthRasterXCommonFinal, truthRasterYCommonFinal, truthRasterZCommonFinal):
        
        self.eastingArray = shiftedEasting
        self.northingArray = shiftedNorthing
        self.crossTrackArray = shiftedCrossTrack
        self.alongTrackArray = shiftedAlongTrack
        self.latArray = shiftedLat
        self.lonArray = shiftedLon
        self.zArray = shiftedVertical
        self.zMslArray = shiftedVerticalMsl
        self.time = time
        self.deltaTime = deltaTime
        self.classification = classification
        self.signalConf = signalConf
        self.zone = zone
        self.hemi= hemi
        self.crossTrack = correctionsCrossTrack
        self.alongTrack = correctionsAlongTrack
        self.easting = correctionsEasting
        self.northing = correctionsNorthing
        self.z = correctionsVertical
        self.mae = correctionsMAE
        self.rmse = correctionsRMSE
        self.me = correctionsME
        self.measX_raster = measRasterXCommonFinal
        self.measY_raster = measRasterYCommonFinal
        self.measZ_raster = measRasterZCommonFinal
        self.truthX_raster = truthRasterXCommonFinal
        self.truthY_raster = truthRasterYCommonFinal
        self.truthZ_raster = truthRasterZCommonFinal
    # endDef
# endClass
    
        
##### Function to read kmlBounds.txt    
def readTruthRegionsTxtFile(kmlBoundsTextFile):
    
    # Initialize output parameters
    regionName = []
    lonMin = []
    lonMax = []
    latMin = []
    latMax = []

    # Open file for reading
    f = open(kmlBoundsTextFile,'r')
    
    # Get all lines of file into a list
    allLines = list(f)
    
    # Close file for reading
    f.close()
        
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
            
        # Endif
    # Endfor
    
    # Call class to populate field names
    kmlInfo = kmlStruct(regionName, lonMin, lonMax, latMin, latMax)
    
    # Return output
    return kmlInfo

# endDef

###### Function to read header truth .mat file
#def readHeaderMatFile(headerFilePath):
#
#    if headerFilePath.endswith('.mat'):
#
#        # Initialize output data
#        coordType = []
#        zone = []
#        hemi = [] 
#        
#        # Get coordinate type from header .mat file (1 = UTM, 2 = Lat/Lon)
#        matData = loadmat(headerFilePath)
#        coordNum= matData['headerData'][0][0][10][0][0]
#        
#        # Convert lat/lon data to UTM coordinates
#        if(coordNum==1):
#            
#            coordType = 'UTM'
#            
#            # Get UTM zone and hemisphere
#            zone = str(matData['headerData'][0][0][7][0][0])
#            hemi = matData['headerData'][0][0][8][0][0]
#            ellipsoid = matData['headerData'][0][0][9][0][0:]
#            
#        else:
#            
#            coordType = 'Lat/Lon'
#            
#        # Endif
#        
#        # Get x/y min/max data and truth tile name for each truth tile
#        xmin = [matData['headerData'][0][i][0][0][0] for i in range(len(matData['headerData'][0]))]
#        xmax = [matData['headerData'][0][i][1][0][0] for i in range(len(matData['headerData'][0]))]
#        ymin = [matData['headerData'][0][i][2][0][0] for i in range(len(matData['headerData'][0]))]
#        ymax = [matData['headerData'][0][i][3][0][0] for i in range(len(matData['headerData'][0]))]
#        tileName = [matData['headerData'][0][i][12][0] for i in range(len(matData['headerData'][0]))]
#        
#        # Store data as object
#        headerData = headerStruct(coordType, zone, hemi, ellipsoid, xmin, xmax, ymin, ymax, tileName)
#        
#        # Return data
#        return headerData
#
#    elif headerFilePath.endswith('.pkl'):
#
#        from icesatReader import read_pickle
#        from icesatUtils import mode, transform
#        import pyproj
#
#        df = read_pickle(headerFilePath)
#
#        proj4_all = list(df['proj4'])
#        proj4_mode = mode(proj4_all)
#        zone0, hemi0, ellipsoid0, coordNum0 = get_params(proj4_mode)
#        tileName_v = list(df['name'])
#
#        coordType0 = 'UTM'
#        if coordNum0 != 1:
#            coordType0 = 'Lat/Lon'
#
#        xmin = np.array(df['xmin'])
#        xmax = np.array(df['xmax'])
#        ymin = np.array(df['ymin'])
#        ymax = np.array(df['ymax'])
#
#        crs_mode = pyproj.CRS.from_proj4(proj4_mode)
#        epsg_mode = crs_mode.to_epsg()
#
#        for k, proj4 in enumerate(proj4_all):
#            crs_las = pyproj.CRS.from_proj4(proj4)
#            epsg_las = crs_las.to_epsg()
#
#            epsg_in = int(epsg_las)
#            epsg_out = int(epsg_mode)
#            if epsg_in != epsg_out:
#                xmin0, ymin0 = transform(epsg_in, epsg_out, xmin[k], ymin[k])
#                xmax0, ymax0 = transform(epsg_in, epsg_out, xmax[k], ymax[k])
#                xmin[k], ymin[k] = xmin0, ymin0
#                xmax[k], ymax[k] = xmax0, ymax0
#
#        xmin_v = list(xmin)
#        xmax_v = list(xmax)
#        ymin_v = list(ymin)
#        ymax_v = list(ymax)
#
#        headerData = headerStruct(coordType0, zone0, hemi0, ellipsoid0, xmin_v, xmax_v, ymin_v, ymax_v, tileName_v)
#
#        return headerData
#    # endIf
## endDef

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

# endDef

#def writeHeaderMatFile(regionName, PKL_DIR, LAS_DIR):
#    writeHeaderFile(regionName, PKL_DIR, LAS_DIR)
## endDef
#
#def writeHeaderFile(regionName, PKL_DIR, LAS_DIR):
#
#    regionName = regionName + '_HeaderData.pkl'
#    file_pkl = os.path.join(PKL_DIR, regionName)
#
#    from icesatReader import write_pickle
#    from icesatIO import readLas
#    import laspy as las
#    import pandas as pd
#
#    files_las = [fn for fn in os.listdir(LAS_DIR) if fn.endswith('.las')]
#    files_las = sorted(files_las)
#    files_las = [os.path.join(LAS_DIR, fn) for fn in files_las]
#
#    df = pd.DataFrame(columns=['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax',
#                                        'pointcount', 'zone', 'hemi', 'ellipsoid',
#                                        'coordNum', 'numRegs', 'name', 'proj4'])
#
#    for f, file_las in enumerate(files_las):
#
#        file_las_sub = os.path.basename(file_las)
#        print(file_las_sub)
#
#        fp = las.file.File(file_las, mode='r')
#
#        min_bounds = fp.get_header().get_min()
#        max_bounds = fp.get_header().get_max()
#        xb_las = np.array([min_bounds[0], max_bounds[0]])
#        yb_las = np.array([min_bounds[1], max_bounds[1]])
#        zb_las = np.array([min_bounds[2], max_bounds[2]])
#
#        xmin, xmax = xb_las[0], xb_las[1]
#        ymin, ymax = yb_las[0], yb_las[1]
#        zmin, zmax = zb_las[0], zb_las[1]
#
#        pointcount = len(fp.points) #readLas(file_las, metadata='count')
#        proj4_las = readLas(file_las, metadata='proj4')
#        zone, hemi, ellipsoid, coordNum = get_params(proj4_las)
#
#        numRegs = 0
#        name = file_las
#
#        df.loc[f] = [xmin, xmax, ymin, ymax, zmin, zmax,
#                            pointcount, zone, hemi, ellipsoid, coordNum, numRegs, name, proj4_las]
#
#        fp.close()
#
#    # save header file
#    # convert_df_to_mat(df_mat, file_mat)
#    write_pickle(df, file_pkl)
#
#    # return df
#
## endDef


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
def readAtl03DataMapping(in_file03, label, return_delta_time=False):
#
# Reads the data from ATL03
#
# Open the file
#
    if not os.path.isfile(in_file03):
        print('File does not exist')
    # endIf
    
    try:
        f = h5py.File(in_file03, 'r')
    except Exception as e:
        print('Python message: %s\n' % e)
        return [], []
    # endTry
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
# delta_time
#
    if(return_delta_time):
        dsname=label+'/geolocation/delta_time'
        if dsname in f:
            delta_time=np.array(f[dsname])
        else:
            delta_time=[]
        # endif
    # endIf
    
    
#
# Close the file
#
    f.close()

    if(return_delta_time):
        return ph_index_beg, segment_id, delta_time
    else:
        return ph_index_beg, segment_id 
    # endIf
  
# endDef


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
    try:
        geoidalHeights = matData['geoid']['geoidalHeight'][0][0]
    except:
        geoidalHeights = matData['geoid']['fin2005n00'][0][0]        
    
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
            headerData = lasFile.headerdate = lasFile.header.get_date()
            
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
        
# endDef

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

### Function to write .las file (2.3 format)
def writeLas(xx,yy,zz,proj,output_file,classification,intensity,signalConf=None,hemi=None,zone=None):
    
    wkt = selectwkt(proj,hemi,zone)
    
    #Create new VLR
    new_vlr = laspy.VLR(user_id = "LASF_Projection",
                           record_id = 2112,
                           record_data = wkt,
                           description = "OGC Coordinate System WKT")
    inVLRs = []
    inVLRs.append(new_vlr)

    #Create new Header
    hdr = laspy.LasHeader(version='1.4')
    hdr.file_sig = 'LASF'
    
    #Create new las file with Header and VLR  
    new_las = laspy.LasData(hdr)
    new_las.header.vlr = inVLRs
    new_las.header.set_wkt = 1
    # outfile = laspy.file.File(output_file, mode="w", header=hdr)
    # outfile.header.vlr = inVLRs
    # outfile.header.set_wkt = 1
    
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
    
    new_las.header.offset = [xmin,ymin,zmin]
    
    #Establish scale
    new_las.header.scale = in_scale
    
    #Write x, y, z data and if available classification and intensity data
    new_las.x = xx
    new_las.y = yy
    new_las.z = zz
    if classification is not None:
        new_las.raw_classification = classification
    if intensity is not None:
        new_las.intensity = intensity
    if signalConf is not None:
        new_las.scan_angle_rank = signalConf
    
    #Close las
    new_las.write(output_file)
    #new_las.close()
    
# endDef
    
### Function to write .las file (1.4 format)
def writeTif(xx,yy,zz,epsg,outputFile,formatArray=False):
    
    # INPUTS:
    # xx - 2D numpy raster or list of x values (if list, then set formatArray=True)
    # yy - 2D numpy raster or list of y values (if list, then set formatArray=True)
    # zz - 2D numpy raster or list of z values (if list, then set formatArray=True)
    # epsg - EPSG code 
    # outputFile - path to output .tif file
    # formatArray - if xx, yy, zz are lists, then use this to format them into rasters
    
    xx = xx[0:10000]
    yy = yy[0:10000]
    zz = zz[0:10000]
    
    if(formatArray):
        
        # Get mean x resolution
        xIn = np.round(xx,0)
        xUnique = np.unique(xIn)
        xDiff = np.diff(xUnique)
        xRes = np.round(np.mean(xDiff),0)
        
        # Get mean y resolution
        yIn = np.round(yy,0)
        yUnique = np.unique(yIn)
        yDiff = np.diff(yUnique)
        yRes = np.round(np.mean(yDiff),0)
        
        # Set total resolution
        resolution = max([xRes,yRes])
        
        # Get raster
        rasterOutput = getRaster(xx, yy, zz, 
                                 resolution, 'Mean', 
                                 fillValue = np.nan, 
                                 time = [], 
                                 xAllArray = [], 
                                 yAllArray = [],
                                 origin=None)
        
        # Store output
        zz = rasterOutput.grid
            
    # Set geotransform
    xmin = np.min(np.ndarray.flatten(xx))
    ymax = np.max(np.ndarray.flatten(yy)) 

    # New Tranform
    new_transform = rasterio.transform.Affine(resolution, 0.0, 
                        xmin, 0.0, -resolution, ymax) 
    
    csr_out = rasterio.crs.CRS.from_epsg(3005)
    
    out_tif = rasterio.open(out_file, 'w', driver='GTiff',
                                        height=resolution, width=resolution,
                                        count=1, dtype=zz.dtype,
                                        crs=csr_out,
                                        transform=new_transform)    

    out_tif.write(zz, 1)
    out_tif.close()

    
##### Functions to write .kml files
def writeKml(lat, lon, time, kmlName):

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
    
    
##### Functions to write .kml files
def writeKml2(lats, lons, ptNames, kmlName):

    # Suppress warnings that may come from simple kml
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    # endif
    
    # Open Simple KML
    kml = simplekml.Kml()
    
    # Open Simple KML style editor
    style = simplekml.Style()
    style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
    style.iconstyle.color = 'ff0000ff'  # Red
    
    # Loop through all lat/lon values and make KML markers
    for i in range(0,len(ptNames)):
        
        ptName = ptNames[i]
        
        # Plot marker points
        pnt = kml.newpoint(name = ptName, coords=[(lons[i], lats[i])])
        pnt.style = style
        
    # EndFor

    # Save output KML file
    kml.save(kmlName)


### Function to write variables to CSV file
def writeArrayToCSV_new(csv_out, namelist, datalist_df):
    
    # Write to csv file
    datalist_df.to_csv(csv_out, mode='w', index=False, header=True)
    
# endDef
    
    
### Function to write variables to CSV file
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
    
# endDef

    
    
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

# endDef
    
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
    
# endDef

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
# endDef 
    
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
        viewer_output = output_folder + '\\Viewer_Online' + in_file03_name + '.html'
        blank = blanki
    else:
        viewer_output = output_folder + '\\Viewer_' + in_file03_name + '.html'
        
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
       
# endDef          
            
def getDEMArrays(in_tif):
    # Open DEM
    
    # Read DEM as NP Array
    # Suppress warnings that may come from rasterio
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    dem = rasterio.open(in_tif)
    data = dem.read(1)
    
    # Get Geotransform Information
    x_left = dem.bounds[0]
    y_top = dem.bounds[3]
    x_res = dem.res[0]
    y_res = -dem.res[1]
    # Generate X Array
    row = (np.arange(data.shape[1])).astype('float')
    for pos in range(0,data.shape[1]):
        row[pos] = x_left + (x_res * row[pos])
    xarr = np.array([row,]*data.shape[0])
    
    #Generate Y Array
    column = (np.arange(data.shape[0])).astype('float')  
    for pos in range(0,data.shape[0]):
        column[pos] = y_top + (y_res * column[pos])  
    yarr = np.array([column,]*data.shape[1]).transpose()
    
    # Return Data (zarray), xarray, and zarray
    return data, xarr, yarr
    
def readDEMepsg(in_file):
    dem = rasterio.open(in_file)
    epsg = str(dem.crs.to_epsg())
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

# endDef
        
### Function to write console and .log messages
def writeLog(message, fileID=False):
    
    # Print message to screen
    print(message)
        
    # Print message to .log file
    if(fileID):   
        print(message, file=fileID)
    # endIf
    
# endDef
    
### Function to load .las/.laz file
def loadLasFile(truthFilePath, atlMeasuredData, rotationData, logFileID=False):
    
    # Read .las file
    lasTruthData = readLas(truthFilePath)
    
    # Find EPSG Code from truth file
    truthHeader = readLasHeader(truthFilePath)
    epsg_truth = truthHeader['epsg'][0]
    
    # Find EPSG Code from input file
    epsg_atl = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
    
    # Get easting/northing
    lasTruthData_x = lasTruthData.x
    lasTruthData_y = lasTruthData.y
    
    if lasTruthData.headerData is not None:
        if len(lasTruthData.headerData) > 0:
            lasTruthData_year = np.zeros(lasTruthData_x.shape[0])
            lasTruthData_year = lasTruthData_year + lasTruthData.headerData[0][0].year
            lasTruthData_month = np.zeros(lasTruthData_x.shape[0])
            lasTruthData_month = lasTruthData_month + lasTruthData.headerData[0][0].month
            lasTruthData_day = np.zeros(lasTruthData_x.shape[0])
            lasTruthData_day = lasTruthData_day + lasTruthData.headerData[0][0].day
        else:
            lasTruthData_year = np.zeros(lasTruthData_x.shape[0])
            lasTruthData_month = np.zeros(lasTruthData_x.shape[0])
            lasTruthData_day = np.zeros(lasTruthData_x.shape[0])

    else:
        lasTruthData_year = np.zeros(lasTruthData_x.shape[0])
        lasTruthData_month = np.zeros(lasTruthData_x.shape[0])
        lasTruthData_day = np.zeros(lasTruthData_x.shape[0])       
        
    lasTruthData_year = np.c_[lasTruthData_year].astype(int)
    lasTruthData_month = np.c_[lasTruthData_month].astype(int)
    lasTruthData_day = np.c_[lasTruthData_day].astype(int)
        
    # Reproject if necessary
    if(epsg_truth == 'None'):
        
        writeLog('      *WARNING: Invalid reference EPSG code, skipping file.', logFileID)
        atlTruthData = False
        
    else:
        
        if(epsg_truth != epsg_atl):
            
            # If EPSG code does not match, reproject to input EPSG code
            writeLog('      *Reference file EPSG code does not match ICESat-2, reprojecting reference file...', logFileID)
            lasTruthData_x, lasTruthData_y = transform(epsg_truth, epsg_atl, lasTruthData.x, lasTruthData.y)
            
        # endIf
        
        # Rotate TRUTH data to CT/AT plane
        lasTruthData_x_newRot, lasTruthData_y_newRot, _, _, _, _ = \
        getCoordRotFwd(lasTruthData_x, lasTruthData_y, 
                       rotationData.R_mat, rotationData.xRotPt, 
                       rotationData.yRotPt, rotationData.desiredAngle)
        
        # Get reference lat/lon
        lasTruth_lat, lasTruth_lon = getUTM2LatLon(lasTruthData_x, lasTruthData_y,
                                                   atlMeasuredData.zone, atlMeasuredData.hemi)
        
        # Store data as object
        atlTruthData = atlTruthStruct(lasTruthData_x, lasTruthData_y, 
                                      lasTruthData_x_newRot, lasTruthData_y_newRot, 
                                      lasTruth_lon, lasTruth_lat,
                                      lasTruthData.z, 
                                      lasTruthData.classification, 
                                      lasTruthData.intensity,
                                      lasTruthData_year,
                                      lasTruthData_month,
                                      lasTruthData_day,                                      
                                      atlMeasuredData.zone, 
                                      atlMeasuredData.hemi,
                                      epsg_atl)
    # endIf
    
    return atlTruthData
    
# endDef
    
### Function to load .tif file
def loadTifFile(truthFilePath, atlMeasuredData, rotationData, outFilePath, logFileID=False):
    
    # Read Tif file
    xarr0, yarr0, zarr, intensity, classification, epsg = formatDEM(truthFilePath)
    
    # Convert ints to floats
    xarr0 = xarr0.astype(float)
    yarr0 = yarr0.astype(float)
    zarr = zarr.astype(float)
    
    # Find EPSG Code from tif
    epsg_truth = 'epsg:' + epsg
    
    # Determine if EPSG Code is the same for the ATL03 Measured
    epsg_atl = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
    
    # Store values
    xarr = xarr0
    yarr = yarr0
    
    # Reproject if necessary
    if(epsg_truth == 'None'):
        
        writeLog('      *WARNING: Invalid reference EPSG code, skipping file.', logFileID)
        atlTruthData = False
        
    else:
        
        if(epsg_truth != epsg_atl):
            
            # If EPSG code does not match, use GDAL Warp to create new Tif
            writeLog('      *Reference file EPSG code does not match ICESat-2, reprojecting reference file...', logFileID)
            xarr, yarr = transform(epsg_truth, epsg_atl, xarr0, yarr0)
            
        
        # Rotate Data for along-track/cross-track
        x_newRot, y_newRot, _, _, _, _ = getCoordRotFwd(xarr, yarr, 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
        
        # Get reference lat/lon
        lat, lon = getUTM2LatLon(xarr, yarr, atlMeasuredData.zone, atlMeasuredData.hemi)
        
        lasTruthData_year = np.zeros(xarr.shape[0])
        lasTruthData_month = np.zeros(xarr.shape[0])
        lasTruthData_day = np.zeros(xarr.shape[0]) 
    
        lasTruthData_year = np.c_[lasTruthData_year].astype(int)
        lasTruthData_month = np.c_[lasTruthData_month].astype(int)
        lasTruthData_day = np.c_[lasTruthData_day].astype(int)
        
        # Store Data as Object
        atlTruthData = atlTruthStruct(xarr, yarr, x_newRot, y_newRot, lon, lat,
                                      zarr, classification, intensity, 
                                      lasTruthData_year, lasTruthData_month,
                                      lasTruthData_day,
                                      atlMeasuredData.zone, 
                                      atlMeasuredData.hemi,
                                      epsg_atl, '')
          
    # endIf
      
    return atlTruthData
        
# endDef
   
### Function to read truth file info
def loadTruthFile(truthFilePath, atlMeasuredData, rotationData, truthFileType, outFilePath, logFileID=False):
    
    # Initialize output
    atlTruthData = []
    
    # Determine which file type to load
    if(('las' in truthFileType.lower()) or ('laz' in truthFileType.lower())):  
        # Get .las header info
        atlTruthData = loadLasFile(truthFilePath, atlMeasuredData, rotationData, logFileID)             
    elif('tif' in truthFileType):
        # Load .tif file
        atlTruthData = loadTifFile(truthFilePath, atlMeasuredData, rotationData, outFilePath, logFileID)
    # endIf
    
    return atlTruthData

# endDef
    
### Function to find and reproject min/max extents in header data
def reprojectHeaderData(truthHeaderDF, atlMeasuredData, logFileID=False):
    
    # Copy dataframe
    truthHeaderNewDF = truthHeaderDF.copy()
    
    # Input EPSG code
    epsg_atl = identifyEPSG(atlMeasuredData.hemi,atlMeasuredData.zone)
    
    # Loop through all truth file EPSG codes
    for i in range(0,len(truthHeaderDF)):
#        writeLog(str(i) + ' out of ' + str(len(truthHeaderDF)), logFileID)
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
                truthHeaderNewDF['xmin'][i] = xout[0]
                truthHeaderNewDF['xmax'][i] = xout[1]
                truthHeaderNewDF['ymin'][i] = yout[0]
                truthHeaderNewDF['ymax'][i] = yout[1]
                truthHeaderNewDF['epsg'][i] = epsg_atl
            
            except:
                
                # Store new extents into output dataframe
                writeLog('WARNING: Cannot reproject data, skipping file: %s' %truthHeaderNewDF['fileName'][i], logFileID)
                truthHeaderNewDF['xmin'][i] = 'None'
                truthHeaderNewDF['xmax'][i] = 'None'
                truthHeaderNewDF['ymin'][i] = 'None'
                truthHeaderNewDF['ymax'][i] = 'None'
                truthHeaderNewDF['epsg'][i] = 'None'
            
            # endTry
                
            
        # endIf
    # endFor
    
    return truthHeaderNewDF

# endDef
    
### Function to check file inputs
def getFilesFromInput(fileInput, truthFileType, logFileID=False):
    
    # Set input file check
    inputFileCheck = False
    
    # Determine input type
    if(isinstance(fileInput, str)):
        fileInput = os.path.normpath(fileInput)
        if('*' in fileInput):
            # Wilcard case
            fileExt = fileInput[-4:]
            if(truthFileType in fileExt):
                inpFiles = glob.glob(fileInput, recursive = True) 
                inputFileCheck = True
            else:
                writeLog('ERROR: File search must contain *%s in string.' %truthFileType, logFileID)
            # endIf
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
                    fileExt = fileInput[-4:]
                    if(truthFileType in fileExt):
                        inpFiles = [fileInput]
                        inputFileCheck = True
                    else:
                        writeLog('ERROR: Input file must be %s format' %truthFileType, logFileID)
                    # endIf
            else:
                writeLog('ERROR: File/Directory does not exist.', logFileID)
            # endIf
        # endIf
    elif(isinstance(fileInput, list)):
        # List case
        inpFiles = [os.path.normpath(x) for x in fileInput]
        inputFileCheck = True
    # endIf
    
    return inpFiles, inputFileCheck
    
# endDef

### Function to read header info for TIF files
def readTifHeader(tifFileInput, outputFilePath = False, logFileID = False):
    
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
               
        2) outputFilePath - option to write tif info outputs to .csv file
           Default is False (no output .csv file)
           Input is a file name or full path and file name (no extension)
    
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
               
        2) <output>.csv file - Optional, <output> name declared by user
            
    """
    
    # Import libraries

    # Get input .tif files
    tifFiles, inputFileCheck = getFilesFromInput(tifFileInput, '.tif', logFileID)
        
    # Set dataframe column names
    columnNames = ['fileName','version','xmin','xmax','ymin','ymax','zmin','zmax', 
                   'nPoints','epsg']
    
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
            
            # todo: test 
            # Get x/y min/max extents  
            dem = rasterio.open(fileName)
            xmin = dem.bounds[0] #left
            ymax = dem.bounds[3] #top
            ymin = dem.bounds[1] #bottom
            xmax = dem.bounds[2] #right

            # Get EPSG code
            epsg = 'epsg:' + readDEMepsg(fileName)
            
            # Set output array
            dataOut = [fileName, version, xmin, xmax, ymin, ymax, zmin, zmax,
                       Npoints, epsg]
            
            # Set dataframe
            dataOutDFsingle = pd.DataFrame(data=dataOut).T
            dataOutDFsingle.columns = columnNames
            
            # Append dataframe
            dataOutDF = dataOutDF.append(dataOutDFsingle, ignore_index=True)
            
        # endFor
        
        # Write output file if requested
        if(outputFilePath):
            try:
                if(os.path.exists(outputFilePath)):
                    writeLog('   Updating output header file: %s' %outputFilePath, logFileID)
                    dataOutDF.to_csv(outputFilePath, mode='a', header=False, index=False)
                else:
                    writeLog('   Writing output header file: %s' %outputFilePath, logFileID)
                    dataOutDF.to_csv(outputFilePath, index=False)
                # endIf
            except:
                writeLog('   WARNING: Could not write output .csv file with las headers.', logFileID)
            # endTry
        # endIf
    # endIf
    
    return dataOutDF

# endDef
    
### Function to read header info for LAS files
def readLasHeader(lasFileInput, outputFilePath = False, logFileID = False):
    
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
               
        2) outputFilePath - option to write las info outputs to .csv file
           Default is False (no output .csv file)
           Input is a file name or full path and file name (no extension)
    
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
               
        2) <output>.csv file - Optional, <output> name declared by user
            
    """
    
    # Import libraries
    import struct
    import numpy as np
    import pandas as pd
    import os

    # Get input .las files
    if isinstance(lasFileInput,str):
        if lasFileInput[-4:] == '.laz':
            lasFiles, inputFileCheck = getFilesFromInput(lasFileInput, '.laz', logFileID)
        elif lasFileInput[-4:] == '.las':
            lasFiles, inputFileCheck = getFilesFromInput(lasFileInput, '.las', logFileID)

    if isinstance(lasFileInput,list):
        if lasFileInput[0][-4:] == '.laz':
            lasFiles, inputFileCheck = getFilesFromInput(lasFileInput, '.laz', logFileID)
        elif lasFileInput[0][-4:] == '.las':
            lasFiles, inputFileCheck = getFilesFromInput(lasFileInput, '.las', logFileID)
        
    # Set dataframe column names
    columnNames = ['fileName','version','xmin','xmax','ymin','ymax','zmin','zmax', 
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
                        # endFor
                        
                        # Convert userID to ASCII char
                        userID = ''.join([chr(int(item)) for item in uID[0]])
                        
                        recordID = struct.unpack('H',file.read(2))[0] # uint16
                        recordLen = struct.unpack('H',file.read(2))[0] # uint16
                    
                        for j in range(0,32):
                            desc[0,j] = struct.unpack('B',file.read(1))[0] # uint8
                        # endFor
                    
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
                                        
                                        # endIf
                                        
                                    elif(sGeoKeys_wKeyID[k,0] == 3072): # UTM file
                                        
                                        if(sGeoKeys_wTIFFTagLocation[k,0] == 0):
                                            projCScode = sGeoKeys_wValue_Offset[k,0]
                                        else :
                                            projCScode = sGeoKeys_wTIFFTagLocation[k,0]
                                        # endIf
                                        
                                        # Get EPSG code
                                        epsg = 'epsg:' + str(int(projCScode))                                       
                             
                                    # endIf
                                # endFor
                            # endIf
                        
                        else:
                            
                            # Move to specified point in file
                            file.seek(recordLen, 1)
                            
                        # endIf        
                    # endFor
                    
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
                        # endFor
                        
                        # Convert userID to ASCII char
                        userID = ''.join([chr(int(item)) for item in uID[0]])
                        
                        recordID = struct.unpack('H',file.read(2))[0] # uint16
                        recordLen = struct.unpack('H',file.read(2))[0] # uint16
                        
                        for j in range(0,32):
                            desc[0,j] = struct.unpack('B',file.read(1))[0] # uint8
                        # endFor
                        
                        # Convert description to ASCII char
                        description = ''.join([chr(int(item)) for item in desc[0]])
                        
                        wkt = np.zeros((1,recordLen))
                        for j in range(0,recordLen):
                            wkt[0,j] = struct.unpack('B',file.read(1))[0] # uint8
                        # endFor
                        
                        # Convert wkt to ASCII char
                        wktAll = ''.join([chr(int(item)) for item in wkt[0]])
                        
                        # Get EPSG code
                        # proj = osr.SpatialReference(wkt=wktAll)
                        # proj = rasterio.crs.CSR.to_wtk()
                        try:
                            las_crs = rasterio.crs.CRS.from_wkt(wktAll)
                            las_crs = las_crs.crs.to_epsg()
                            epsg = 'epsg:' + str(las_crs)
                        except:
                            print('EPSG code not found')
                            pass
                        # endTry
    
#                        # Get ellipsoid
#                        ellipsoid = (findStr(wktAll,'PROJCS["','/ UTM zone')).strip()
#                        ellipsoid = ''.join(ellipsoid.split())
#                        
#                        # Get UTM Zone/Hemi
#                        utmInfo = (findStr(wktAll,'/ UTM zone','",GEOGCS["WGS')).strip()
#                        zone = str(int(utmInfo[0:-1])).zfill(2)
#                        hemi = utmInfo[-1:]
                        
                    # endFor
                # endIf
            # endWith
            
#            # Define EPSG code
#            if(hemi == 'N'):
#                epsg = 'epsg:326'
#            elif(hemi == 'S'):
#                epsg = 'epsg:327'
#            # endIf
#            epsg = epsg + zone
            
            # Set output array
            dataOut = [fileName, version, xmin, xmax, ymin, ymax, zmin, zmax,
                       Npoints, epsg]
            
            # Set dataframe
            dataOutDFsingle = pd.DataFrame(data=dataOut).T
            dataOutDFsingle.columns = columnNames
            
            # Append dataframe
            dataOutDF = dataOutDF.append(dataOutDFsingle, ignore_index=True)
        
        # endFor
        
        # Write output file if requested
        if(outputFilePath):
            try:
                if(os.path.exists(outputFilePath)):
                    writeLog('   Updating output header file: %s' %outputFilePath, logFileID)
                    dataOutDF.to_csv(outputFilePath, mode='a', header=False, index=False)
                else:
                    writeLog('   Writing output header file: %s' %outputFilePath, logFileID)
                    dataOutDF.to_csv(outputFilePath, index=False)
                # endIf
            except:
                writeLog('   WARNING: Could not write output .csv file with las headers.', logFileID)
            # endTry
        # endIf
    # endIf
    
    return dataOutDF

# endDef
    
### Function to read header info
def getTruthHeaders(truthFilePath, truthFileType, logFileID=False):
    
    # Initialize output
    truthHeaderDF = []
    writeNew = True
    appendNew = False
        
    # Create header file output name
    headerFileName = 'phoReal_headers.csv'
    truthFileDir = ntpath.dirname(truthFilePath[0])
    if os.name == 'nt':
        headerFilePath = os.path.normpath(truthFileDir + '\\' + headerFileName)
    else:
        headerFilePath = os.path.normpath(truthFileDir + '/' + headerFileName)

    # Check if header file exists
    if(os.path.exists(headerFilePath)):
        writeLog('   Previous header file exists', logFileID)
        stored_truthHeaderDF = pd.read_csv(headerFilePath)
        stored_truthFilePath = (stored_truthHeaderDF['fileName']).to_list()
        stored_filesTrue = np.isin(truthFilePath, stored_truthFilePath)
        stored_filesSame = all(stored_filesTrue)
        
        # Determine if files are the same
        if(stored_filesSame):
            if(len(truthFilePath)==len(stored_truthFilePath)):
                writeLog('   Stored header file is up-to-date, using this file', logFileID)
                truthHeaderDF = stored_truthHeaderDF
            elif(len(truthFilePath)<len(stored_truthFilePath)):
                writeLog('   Stored header file is ahead of reference directory, using only reference files selected', logFileID)
                indsTrue = np.isin(stored_truthFilePath, truthFilePath)
                truthHeaderDF = stored_truthHeaderDF[indsTrue]
                truthHeaderDF.reset_index(inplace=True)
            # endIf
            writeNew = False
        else:
            writeLog('   Stored header file is not up-to-date, updating...', logFileID)
            newInds = stored_filesTrue==False
            truthFilePath = (np.array(truthFilePath)[newInds]).tolist()
            appendNew = True
        # endIf
    # endIf
    
    if(writeNew):
        
        # Message to user
        numFiles = len(truthFilePath)
        writeLog('   Reading header info for all %d reference files in: %s' %(numFiles, truthFileDir), logFileID)
    
        if(('las' in truthFileType.lower()) or ('laz' in truthFileType.lower())):  
            # Get .las header info
            truthHeaderNewDF = readLasHeader(truthFilePath, headerFilePath, logFileID)              
        elif('tif' in truthFileType):
            # Load .tif file
            truthHeaderNewDF = readTifHeader(truthFilePath, headerFilePath, logFileID)  
        # endIf
        
        # Append dataframes together if necessary
        if(appendNew):
            frames = [stored_truthHeaderDF, truthHeaderNewDF]
            truthHeaderDF = pd.concat(frames, ignore_index=True)
        else:
            truthHeaderDF = truthHeaderNewDF
        # endIf
        
    # endIf
    
    # Write a blank line for readability
    writeLog('', logFileID)
    
    return truthHeaderDF

# endDef
    
### Function to apply geoid model to truth buffer
def applyGeoidCorrection(atlTruthData, truthHeaderNewDF, logFileID=False):

    # Apply geoid model Z corrections if necessary
    regionName = truthHeaderNewDF['fileName'][0]
    
    # Initialize geoidFile variable
    geoidFile = False
    
    # Determine geoid model to use
    if('sonoma' in regionName.lower() or 'indiana' in regionName.lower()):     
        geoidFile = 'geoid12b_latlon.mat'
    elif('finland' in regionName.lower()):
        geoidFile = 'geoidFin2005N00_latlon.mat'
    elif('rangiora' in regionName.lower() or 'wellington' in regionName.lower()):
        geoidFile = 'geoidNZ_latlon.mat'
    elif('fthood' in regionName.lower()):
        geoidFile = 'geoid_egm2008.mat'
    # endif
    
    if(geoidFile):
    
        # Print status message
        writeLog('', logFileID)
        writeLog('   STATUS: Applying Truth Z Correction (Geoid File = %s).' %geoidFile, logFileID)
            
        # Load Geoid file
        geoidData = readGeoidFile(geoidFile)

        # Get geoidal heights and add to orthometric heights
        atlTruthData = getGeoidHeight(geoidData,atlTruthData)
        
    # EndIf
    
    return atlTruthData

# endDef
    
### Function to find text between matching strings
def findStr(textIn,textLeft,textRight):

    textOut = textIn[textIn.index(textLeft)+len(textLeft):textIn.index(textRight)]
    
    return textOut

# endDef

### Function to do quick filter operation
def quickFilter(atlMeasuredData, xarr, yarr, zarr, intensity, classification):

    # Get buffered data
    buffer = 100
    xmin = np.min(atlMeasuredData.easting.flatten()) - buffer
    xmax = np.max(atlMeasuredData.easting.flatten()) + buffer
    ymin = np.min(atlMeasuredData.northing.flatten()) - buffer
    ymax = np.max(atlMeasuredData.northing.flatten()) + buffer
    
    # Get indices outside of buffer
    quickfilter = np.where((xarr > xmin) & (xarr < xmax) & (yarr > ymin) & (yarr < ymax))
    
    # Filter data outside of buffer
    xarr = xarr[quickfilter]
    yarr = yarr[quickfilter]
    zarr = zarr[quickfilter]
    intensity = intensity[quickfilter]
    classification = classification[quickfilter]
    
    # Filter out bad Z values in DTM
    zGoodInds = zarr>0
    xarr = xarr[zGoodInds]
    yarr = yarr[zGoodInds]
    zarr = zarr[zGoodInds]
    classification = classification[zGoodInds]
    intensity = intensity[zGoodInds]
    
    return xarr, yarr, zarr, intensity, classification

# endDef
    
### Fuction to create truth buffer around ground track
def makeBuffer(atlTruthData, atlMeasuredData, rotationData, buffer):
    
    # Run superfilter to get buffer points
    atlTruthDataFiltered, _ = superFilter(atlMeasuredData, atlTruthData, xBuf = buffer, classCode = [])
                    
    # Store data as object
    atlTruthDataBuffer = atlTruthStruct(atlTruthDataFiltered.easting, 
                                        atlTruthDataFiltered.northing, 
                                        atlTruthDataFiltered.crossTrack,
                                        atlTruthDataFiltered.alongTrack,
                                        atlTruthDataFiltered.lon, 
                                        atlTruthDataFiltered.lat, 
                                        atlTruthDataFiltered.z, 
                                        atlTruthDataFiltered.classification, 
                                        atlTruthDataFiltered.intensity, 
                                        atlTruthDataFiltered.year, 
                                        atlTruthDataFiltered.month, 
                                        atlTruthDataFiltered.day, 
                                        atlTruthDataFiltered.zone, 
                                        atlTruthDataFiltered.hemi,
                                        atlTruthDataFiltered.epsg)
    
    return atlTruthDataBuffer
            
# endDef
    
### Fuction to create truth buffer around ground track
def makeBuffer_legacy(atlTruthData, atlMeasuredData, rotationData, buffer):
    
    # Get MEASURED rotated buffer bounds data
    xRotL = atlMeasuredData.crossTrack - buffer
    xRotR = atlMeasuredData.crossTrack + buffer
    yRot = atlMeasuredData.alongTrack

    # Rotate TRUTH lasData to CT/AT plane
    xRotLasData, yRotLasData,  _, _, _, _ = getCoordRotFwd(atlTruthData.easting, 
                                                           atlTruthData.northing, 
                                                           rotationData.R_mat, 
                                                           rotationData.xRotPt, 
                                                           rotationData.yRotPt, 
                                                           rotationData.desiredAngle)

    # Find TRUTH lasData points inside TRUTH buffer
    yRotLocalInds = (yRot >= yRotLasData.min()) & (yRot <= yRotLasData.max()) 
    xRotLocalL = xRotL[yRotLocalInds]
    xRotLocalR = xRotR[yRotLocalInds]
    xyBufferInds = (xRotLasData >= xRotLocalL.min()) & (xRotLasData <= xRotLocalR.max())
    
    # Extract TRUTH lasData points inside TRUTH buffer
    subData_easting = atlTruthData.easting[xyBufferInds]
    subData_northing = atlTruthData.northing[xyBufferInds]
    subData_crossTrack = atlTruthData.crossTrack[xyBufferInds]
    subData_alongTrack = atlTruthData.alongTrack[xyBufferInds]
    subData_lat = atlTruthData.lat[xyBufferInds]
    subData_lon = atlTruthData.lon[xyBufferInds]
    subData_z = atlTruthData.z[xyBufferInds]
    subData_classification = atlTruthData.classification[xyBufferInds]
    subData_intensity = atlTruthData.intensity[xyBufferInds]
    subData_year = atlTruthData.year[xyBufferInds]
    subData_month = atlTruthData.month[xyBufferInds]
    subData_day = atlTruthData.day[xyBufferInds]



    # Store data as object
    atlTruthDataBuffer = atlTruthStruct(subData_easting, subData_northing, 
                                        subData_crossTrack, subData_alongTrack,
                                        subData_lon, subData_lat,
                                        subData_z, 
                                        subData_classification, 
                                        subData_intensity,
                                        subData_year,
                                        subData_month,
                                        subData_day,
                                        atlTruthData.zone, 
                                        atlTruthData.hemi,
                                        atlTruthData.epsg)
    
    return atlTruthDataBuffer
            
# endDef
    
### Function to find which truth tiles ICESat-2 crosses over
def findMatchingTruthFiles(truthHeaderNewDF, atlMeasuredData, rotationData, buffer):
       
    # Get MEASURED rotated buffer bounds data
    xRotL = atlMeasuredData.crossTrack - buffer
    xRotR = atlMeasuredData.crossTrack + buffer
    yRot  = atlMeasuredData.alongTrack
    
    # Get MEASURED buffer bounds data in easting/northing plane
    xL, yL,  _, _, _ = getCoordRotRev(xRotL, yRot, rotationData.R_mat, rotationData.xRotPt, rotationData.yRotPt)
    xR, yR,  _, _, _ = getCoordRotRev(xRotR, yRot, rotationData.R_mat, rotationData.xRotPt, rotationData.yRotPt)
    
    # Rotate truth header file min/max x,y points to CT/AT plane
    xMinyMinHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(((truthHeaderNewDF['xmin']).to_numpy()).astype('float'), 
                                                        ((truthHeaderNewDF['ymin']).to_numpy()).astype('float'), 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
    
    xMinyMaxHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(((truthHeaderNewDF['xmin']).to_numpy()).astype('float'), 
                                                        ((truthHeaderNewDF['ymax']).to_numpy()).astype('float'), 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
    xMaxyMinHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(((truthHeaderNewDF['xmax']).to_numpy()).astype('float'), 
                                                        ((truthHeaderNewDF['ymin']).to_numpy()).astype('float'), 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
    xMaxyMaxHeaderRotX, _,  _, _, _, _ = getCoordRotFwd(((truthHeaderNewDF['xmax']).to_numpy()).astype('float'), 
                                                        ((truthHeaderNewDF['ymax']).to_numpy()).astype('float'), 
                                                        rotationData.R_mat, 
                                                        rotationData.xRotPt, 
                                                        rotationData.yRotPt, 
                                                        rotationData.desiredAngle)
    
    # Find min/max x/y header points inside min/max x buffer points
    xMinyMinXPtsInBuffer = (xMinyMinHeaderRotX >= xRotL.min()) & \
    (xMinyMinHeaderRotX <= xRotR.max())
    xMinyMaxXPtsInBuffer = (xMinyMaxHeaderRotX >= xRotL.min()) & \
    (xMinyMaxHeaderRotX <= xRotR.max())
    xMaxyMinXPtsInBuffer = (xMaxyMinHeaderRotX >= xRotL.min()) & \
    (xMaxyMinHeaderRotX <= xRotR.max())
    xMaxyMaxXPtsInBuffer = (xMaxyMaxHeaderRotX >= xRotL.min()) & \
    (xMaxyMaxHeaderRotX <= xRotR.max())
    
    # Get header points inside MEASURED buffer points
    xHeaderPtsInBuffer = np.c_[xMinyMinXPtsInBuffer | xMinyMaxXPtsInBuffer| xMaxyMinXPtsInBuffer | xMaxyMaxXPtsInBuffer]
    
    # Find min/max x buffer points inside min/max x/y header points
    xyMinMaxHeaderRot = np.column_stack((xMinyMinHeaderRotX,xMinyMaxHeaderRotX,xMaxyMinHeaderRotX,xMaxyMaxHeaderRotX))
    xyMinHeaderRot = np.c_[np.amin(xyMinMaxHeaderRot,axis=1)]
    xyMaxHeaderRot = np.c_[np.amax(xyMinMaxHeaderRot,axis=1)]
    xMinBufferPtsInFile = (xRotL.min() >= xyMinHeaderRot) & (xRotL.min() <= xyMaxHeaderRot)
    xMaxBufferPtsInFile = (xRotR.max() >= xyMinHeaderRot) & (xRotR.max() <= xyMaxHeaderRot)
    
    # Get MEASURED buffer points inside header points
    xBufferPtsInHeader = xMinBufferPtsInFile | xMaxBufferPtsInFile
    
    # Get any points where buffer points are inside header and vice versa
    xPtsInFile = np.ravel(np.logical_or(xHeaderPtsInBuffer, xBufferPtsInHeader))
        
    # Get matching truth file names and x/y min/max points
    matchingFilesPre = ((truthHeaderNewDF['fileName'][xPtsInFile]).to_numpy()).astype('str')
    matchingHeaderXmin = ((truthHeaderNewDF['xmin'][xPtsInFile]).to_numpy()).astype('float')
    matchingHeaderXmax = ((truthHeaderNewDF['xmax'][xPtsInFile]).to_numpy()).astype('float')
    matchingHeaderYmin = ((truthHeaderNewDF['ymin'][xPtsInFile]).to_numpy()).astype('float')
    matchingHeaderYmax = ((truthHeaderNewDF['ymax'][xPtsInFile]).to_numpy()).astype('float')
    
    fileNumsAll = np.arange(0,len(truthHeaderNewDF))
    allFileInds = fileNumsAll[xPtsInFile]
    
    # Get all MEASURED x,y buffer points
    xAll = np.concatenate((xL, xR))
    yAll = np.concatenate((yL, yR))
    
    matchTF = []
    if(len(matchingFilesPre)>0):
        for i in range(0,len(matchingFilesPre)):
            
            # Determine TRUTH files where MEASURED data actually crosses over
            xPtsInFile = (xAll >= matchingHeaderXmin[i]) & (xAll <= matchingHeaderXmax[i])
            yPtsInFile = (yAll >= matchingHeaderYmin[i]) & (yAll <= matchingHeaderYmax[i])
            anyPtsInFile = any(xPtsInFile & yPtsInFile)
    
            # If a TRUTH file is a match, use it
            if(anyPtsInFile):
                matchTF.append(True)
            else:
                matchTF.append(False)
            # endIf
            
        # endFor
    # endIf
    
    if(len(matchTF)>0):
        matchTF = np.array(matchTF)
        matchingFiles = matchingFilesPre[matchTF]
    else:
        matchingFiles = []
    # endIf
    
    matchingFileInds = allFileInds[matchTF]
    
    return matchingFiles, matchingFileInds

# endDef
    
### Function to get truth file paths from GUI input
def getTruthFilePaths(userInput, fileExt, logFileID=False):
    
    # INPUTS:
    # userInput - string, list, tuple of path(s)
    # fileExt - '.las', '.tif', '.h5' etc.
    # logFileID - a file identifier to write to a file if desired
    
    # Initialize output
    filePaths = []
    
    # Convert tuple to list if necessary
    if(isinstance(userInput, tuple)):
        userInput = list(userInput)
    # endIf
    
    elif(isinstance(userInput, list)):
        if(len(userInput)==1):
            userInput = userInput[0]
        # endIf
    # endIf
    
    # Get truth file paths based on file, files, or directory
    if(isinstance(userInput, str)):
        if(os.path.exists(userInput)):
            if(os.path.isfile(userInput)):
                # Get truth file name and extension
                filePaths = [os.path.normpath(userInput)]
            else:
                # Get full paths to input files
                if os.name == 'nt':
                    strSearch = os.path.normpath(userInput + '\\*' + fileExt)
                else:
                    strSearch = os.path.normpath(userInput + '/*' + fileExt)
                filePaths = glob.glob(strSearch, recursive = True)   
            # endIf     
        # endIf
    elif(isinstance(userInput, list)):
        filePathsList = [os.path.normpath(x) for x in userInput]
        filePathsTrue = all([os.path.isfile(x) for x in filePathsList])
        if(filePathsTrue):
            filePaths = filePathsList
        # endIf
    # endIf

    return filePaths

# endDef 
    
# UNIT TEST


'''
beamNumToGT will convert a spot beam number into a groundtrack
Spot beams will be either strings in range 1-6, e.g., "1","2","3","4","5,"6"
Input: ATL03 file, beam num string
Output: a GT that beamnumber is assoicated with as a string
'''
def beamNumToGT(atlfilepath, beamNum):
    gt_out = []
    gt_list = ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']
    fp = h5py.File(atlfilepath,'r')
    for gt in gt_list:
        try:
            fp_a = fp[gt].attrs
            spot_number = (fp_a['atlas_spot_number']).decode()
        except:
            spot_number = 'none'
        # endTry
        
        if(any(np.isin(beamNum,spot_number))):
            gt_out.append(gt.upper())
        # endIf
        
    return gt_out
# endDef
    
def GtToBeamNum(atlfilepath, GT):
    fp = h5py.File(atlfilepath,'r')
    gt = GT.lower()
    fp_a = fp[gt].attrs
    beamNum = (fp_a['atlas_spot_number']).decode()
    
    return beamNum
# endDef
    

'''
swbeamToGT will ID which groundtrack by using a beam group (1,2,3) and a 
strongbeam or weakbeam.
Input: ATL03 file, group num string ('1','2','3') and beam type ('strong','weak')
Output: a GT that beamnumber is assoicated with as a string
'''
def swbeamToGT(atlfilepath, groupNum, sw):
    gt_out = []
    gt_list = ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']
    fp = h5py.File(atlfilepath,'r')
    for gt in gt_list:
        try:
            fp_a = fp[gt].attrs
            beam_type = (fp_a['atlas_beam_type']).decode()
        except:
            beam_type = 'none'
        # endTry
        
        if(any(np.isin(sw,beam_type) & np.isin(groupNum, gt[2]))):
            gt_out.append(gt.upper())
        # endIf
    # endFor
    
    return gt_out
# endDef
    
def GtToBeamSW(atlfilepath, GT):
    fp = h5py.File(atlfilepath,'r')
    gt = GT.lower()
    fp_a = fp[gt].attrs
    beamSW = (fp_a['atlas_beam_type']).decode()
    return beamSW
# endDef

if __name__ == "__main__":
    print("Test")
    
# endIf
