# -*- coding: utf-8 -*-
"""
Script that contains utliity functions for PhoREAL

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
import numpy as np
import sys
import pyproj as proj
import warnings
from scipy import interpolate
import pandas as pd
import h5py
import csv


# Object for getNameParts function
class fileStruct:
        
    # Define class with designated fields
    def __init__(self, atlVersion, year, month, day, hour, minute, second, trackNum, unknown, releaseNum, incrementNum):
            
        self.atlVersion = atlVersion
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        self.trackNum = trackNum
        self.unknown = unknown
        self.releaseNum = releaseNum
        self.incrementNum = incrementNum
            

# Object for gridMetricNew function
class GridStruct:
    def __init__(self, x, y, grid, time):
        self.x = x
        self.y = y
        self.grid = grid
        self.t = time
        
        
##### Function for reading parts of an .h5 file name
def getNameParts(h5FileName):
        
    # Split file name by underscores
    nameSplit = h5FileName.split('_')
    
    # Get ATL version
    atlVersion = nameSplit[0]
    
    # Get Year, month, day
    year = nameSplit[1][0:4]
    month = nameSplit[1][4:6]
    day = nameSplit[1][6:8]
    
    # Get Hour, minute, second
    hour = nameSplit[1][8:10]
    minute = nameSplit[1][10:12]
    second = nameSplit[1][12:14]
    
    # Get other details
    trackNum = nameSplit[2][0:4]
    unknown = nameSplit[2][4:9]
    releaseNum = nameSplit[3]
    incrementNum = nameSplit[4]
    
    # Get data into class structure
    fileInfo = fileStruct(atlVersion, year, month, day, hour, minute, second, trackNum, unknown, releaseNum, incrementNum)
    
    # Return class
    return fileInfo


##### Function to represent 2 numbers as 1 unique number
def cantorPairing(arrayIn):
    
    # Do Cantor pairing algorithm to express two values as one unique value
    vectorOut = 0.5 * (arrayIn[:,0] + arrayIn[:,1]) * (arrayIn[:,0] + arrayIn[:,1] + 1) + arrayIn[:,1]
    
    return vectorOut


##### Function to determine members of one array in another
def ismember(a_vec, b_vec, methodType = 'normal'):
    
    """ MATLAB equivalent ismember function """
    
    # Do Cantor pairing if necessary
    if(methodType.lower() == 'rows'):
        a_vec = cantorPairing(a_vec)
        b_vec = cantorPairing(b_vec)
    # EndIf
    
    # Find which values in set array belong to subset
    matchingTF = np.isin(a_vec,b_vec)
    common = a_vec[matchingTF]
    common_unique, common_inv  = np.unique(common, return_inverse=True)     # common = common_unique[common_inv]
    b_unique, b_ind = np.unique(b_vec, return_index=True)  # b_unique = b_vec[b_ind]
    common_ind = b_ind[np.isin(b_unique, common_unique, assume_unique=True)]
    matchingInds = common_ind[common_inv]
    
    return matchingTF, matchingInds


##### Function to determine intersection of two arrays
def getIntersection(a_vec, b_vec):
    
    # Get set intersection (common values) of two arrays
    res_set = set(map(tuple, a_vec)) & set(map(tuple, b_vec)) 
    commonVals = np.array(list(map(list, res_set))) 
    
    if(commonVals.any()):
    
        # Get indices of common values for each array
        _, a_inds = ismember(commonVals, a_vec, 'rows')
        _, b_inds = ismember(commonVals, b_vec, 'rows')
    
    else:
        
        a_inds = []
        b_inds = []
        
    # EndIf
    
    return commonVals, a_inds, b_inds


##### Function to map ATL08 to ATL03 class photons
def getAtl08Mapping(atl03_ph_index_beg, atl03_segment_id, atl08_classed_pc_indx, atl08_classed_pc_flag, atl08_segment_id):
      
    # Get ATL03 data
    indsNotZero = atl03_ph_index_beg != 0
    atl03_ph_index_beg = atl03_ph_index_beg[indsNotZero];
    atl03_segment_id = atl03_segment_id[indsNotZero];
    
    # Find ATL08 segments that have ATL03 segments
    atl03SegsIn08TF, atl03SegsIn08Inds = ismember(atl08_segment_id,atl03_segment_id)
    
    # Get ATL08 classed indices and values
    atl08classed_inds = atl08_classed_pc_indx[atl03SegsIn08TF]
    atl08classed_vals = atl08_classed_pc_flag[atl03SegsIn08TF]

    # Determine new mapping into ATL03 data
    atl03_ph_beg_inds = atl03SegsIn08Inds;
    atl03_ph_beg_val = atl03_ph_index_beg[atl03_ph_beg_inds];
    newMapping = atl08classed_inds + atl03_ph_beg_val - 2;
    
    # Get max size of output array
    sizeOutput = newMapping[-1]
    
    # Pre-populate all photon classed array with zeroes
    allph_classed = np.zeros(sizeOutput + 1).astype(int)
    
    # Populate all photon classed array from ATL08 classifications
    allph_classed[newMapping] = atl08classed_vals;
    
    # Return all photon classed array
    return allph_classed



##### Functions to convert Lat/Lon to UTM and vice versa
'''ATL Geographic Coordinate System Converter

ATL GCS Converter will translate an ATL groundtrack to different coordinate systems.
Mainly this is to convert Lonitute and latitude in WGS84 coordinates, which are
radial, to cartesian coordinates (e.g., UTM) for easy processing, and if required,
back to WGS84.  Primary package required is pyproj.


TODO: Future functionality - Exception UTM grids

Each UTM grid is standardized across the globe but there a few exceptions. Zones
    31V, 32V, 31X, 33X, 35X, and 37X are irregularly sized.  Current code does 
    not make these exceptions.  This should not make too big of a difference for
    now but will become an issue if code needs to work with other GIS.
    For more information:
    https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system
    https://upload.wikimedia.org/wikipedia/commons/e/ed/Utm-zones.jpg
    https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#/ media/File:Modified_UTM_Zones.png

TODO: Add checks on find_utm_zone_point and find_utm_zone_arr to ensure that
the points being entered into the function are in lon/lats by ensuring that the
points are do not exceed typical lon/lat ranges.

Notes:
    
UTM Zone basics:
    https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system
    https://www.e-education.psu.edu/natureofgeoinfo/c2_p23.html

Antartic coordinate systems:
    https://epsg.io/3412 (Datum not specified)
    https://epsg.io/3976 (Official NSIDC coordinate system)
    https://epsg.io/3031 (Recommended by Mike)
    
Artic coordinate systems:
    https://epsg.io/3995
    https://epsg.io/3413 (Official NSIDC coordinate system: 
    https://nsidc.org/data/oib/epsg_3413.html)
        

 '''

# Find UTM zone for individual lon/lat.
def find_utm_zone_point(lon, lat):
    utm_band = str(int((np.floor((lon + 180) / 6 ) % 60) + 1))
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
        print(utm_band)
    if lat >= 0 and lat < 84:
        epsg_code = '326' + utm_band
    elif lat >= 84:
        epsg_code = '3413'
    elif lat <= -80:
        epsg_code = '3976'
    else:
        epsg_code = '327' + utm_band
    return epsg_code


# Find UTM zone for numpy array of lon/lat.
def find_utm_zone_arr(lon, lat):
    arctic = len(lat[lat >= 84])
    nhem = len(lat[(lat < 84) & (lat >= 0)])
    shem = len(lat[(lat > -80) & (lat < 0)])
    antarctic = len(lat[lat < -80])
    if arctic > nhem: #Arctic case
        epsg_code = '3413'
    elif antarctic > shem: #Antarctic case
        epsg_code = '3976'
    else: #If not Arctic or Antarctic it is within UTM Zone
        tz = np.floor((((lon + 180) / 6 ) % 60) + 1).astype(int)
        zone = np.unique(tz)
        if len(zone) == 1:
            zone = zone[0]
        elif len(zone) == 2:
            z1 = len(tz[tz == zone[0]])
            z2 = len(tz[tz == zone[1]])
            if z1 > z2:
                zone = zone[0]
            else:
                zone = zone[1]
        elif len(zone) == 3:
            zone = zone[1]
        else:
            sys.exit("Error: Input ground track present in more than 3 UTM zones. \
                     \nRecommend manually selecting GCS.")
        if nhem >= shem:
            epsg_code = '326' + str(zone)
        else:
            epsg_code = '327' + str(zone)
    return epsg_code


# Transform GCS/PCS based on EPSG and x/y. 
def transform(epsg_in, epsg_out, x, y):
	crs_in = proj.Proj(init = epsg_in)
	crs_out = proj.Proj(init = epsg_out)
	xx, yy= proj.transform(crs_in, crs_out, x, y)
	return xx,yy


# Transform from lon/lat to given EPSG. 
def wgs84_to_epsg_transform(epsg_out, lon, lat):
    epsg_in = 'epsg:4326'
    epsg_out = ('epsg:{0}'.format(str(epsg_out)))
    xx, yy= transform(epsg_in, epsg_out, lon, lat)
    return xx,yy


# Calls functions to find EPSG code and perform GCS transform automatically.
def wgs84_to_utm_find_and_transform(lon, lat):
    epsg_out = find_utm_zone_arr(lon, lat)
    xx, yy = wgs84_to_epsg_transform(epsg_out, lon, lat)
    return xx,yy,epsg_out


# Inputs: lon, lat, (UTM zone), (UTM hemisphere)
def getLatLon2UTM(*args):
    
    # Set EPSG code for lat/lon coords
    epsg_in = 'epsg:4326'
    
    # Get lats/lons
    lon = args[0]
    lat = args[1]
    
    # Call function based on number of input args
    if(len(args) > 2):
        
        # Get zone/hemi
        zone = args[2]
        hemi = args[3]
        
        # Get EPSG out code for UTM coords
        if(hemi=='N'):
            epsg_out = 'epsg:326' + zone
        else:
            epsg_out = 'epsg:327' + zone
        # endif
        
        # Call transform function
        xx, yy = transform(epsg_in, epsg_out, lon, lat)
        
    else:
        
        # Get UTM coords
        xx, yy, epsg_out = wgs84_to_utm_find_and_transform(lon, lat)
        
        # Store zone
        zone = epsg_out[3:]
        
        # Store hemisphere
        if(epsg_out[0:3]==326):   
            hemi = 'N'  
        else:
            hemi = 'S'
        # endif
        
    # endif
    
    # Return output
    return xx, yy, zone, hemi


# Function to convert UTM to lat/lon
def getUTM2LatLon(x,y,zone,hemi):
    
    # Set EPSG code for lat/lon coords
    epsg_out = 'epsg:4326'
    
    # Get EPSG code for UTM coords
    if(hemi=='N'):
        epsg_in = 'epsg:326' + zone
    else:
        epsg_in = 'epsg:327' + zone
    # endif
    
    # Call transform function
    lon, lat = transform(epsg_in, epsg_out, x, y)
        
    return lat, lon


# Identifies midpoint for a given array.
def getMidpoint(arr):
    n = arr.shape[0] / 2.0
    n_int = int(n)
    if n % 2 == 0:
        return (arr[n_int] + arr[n_int - 1]) / 2
    else:
        return arr[n_int]
    
    
##### Functions to convert from Easting/Northing frame to Cross-Track/Along-Track frame and vice versa
def getCoordRotFwd(xIn,yIn,R_mat,xRotPt,yRotPt,desiredAngle):
   
    # Get shape of input X,Y data
    xInShape = np.shape(xIn)
    yInShape = np.shape(yIn)
    
    # If shape of arrays are (N,1), then make them (N,)
    xIn = xIn.ravel()
    yIn = yIn.ravel()
    
    # Suppress warnings that may come from np.polyfit
    if not sys.warnoptions:
        warnings.simplefilter("ignore")
    # endif
    
    # If Rmatrix, xRotPt, and yRotPt are empty, then compute them
    if(len(R_mat)==0 and len(xRotPt)==0 and len(yRotPt)==0):
        
        # Get current angle of linear fit data
        x1 = xIn[0]
        x2 = xIn[-1]
        y1 = yIn[0]
        y2 = yIn[-1]
        # endif
        deltaX = x2 - x1
        deltaY = y2 - y1
        theta = np.arctan2(deltaY,deltaX)
        
        # Get angle to rotate through
        phi = np.radians(desiredAngle) - theta
        
        # Get rotation matrix
        R_mat = np.matrix(np.array([[np.cos(phi), -np.sin(phi)],[np.sin(phi), np.cos(phi)]]))
        
        # Get X,Y rotation points
        xRotPt = x1
        yRotPt = y1
    
    else:
        
        # Get angle to rotate through
        phi = np.arccos(R_mat[0,0])
    
    # endif
    
    # Translate data to X,Y rotation point
    xTranslated = xIn - xRotPt
    yTranslated = yIn - yRotPt
    
    # Convert np array to np matrix
    xTranslated_mat = np.matrix(xTranslated)
    yTranslated_mat = np.matrix(yTranslated)
    
    # Get shape of np X,Y matrices
    (xTranslated_matRows,xTranslated_matCols) = xTranslated_mat.shape
    (yTranslated_matRows,yTranslated_matCols) = yTranslated_mat.shape
    
    # Make X input a row vector
    if(xTranslated_matRows > 1):
        xTranslated_mat = np.transpose(xTranslated_mat)
    #endif
    
    # Make Y input a row vector
    if(yTranslated_matRows > 1):
        yTranslated_mat = np.transpose(yTranslated_mat)
    #endif
    
    # Put X,Y data into separate rows of matrix
    xyTranslated_mat = np.concatenate((xTranslated_mat,yTranslated_mat))
    
    # Compute matrix multiplication to get rotated frame
    measRot_mat = np.matmul(R_mat,xyTranslated_mat)
                            
    # Pull out X,Y rotated data
    xRot_mat = measRot_mat[0,:]
    yRot_mat = measRot_mat[1,:]
    
    # Convert X,Y matrices back to np arrays for output
    xRot = np.array(xRot_mat)
    yRot = np.array(yRot_mat)
    
    # Make X,Y rotated output the same shape as X,Y input
    xRot = np.reshape(xRot,xInShape)
    yRot = np.reshape(yRot,yInShape)
    
    # Reset warnings 
    warnings.resetwarnings()
                   
    # Return outputs
    return xRot, yRot, R_mat, xRotPt, yRotPt, phi


def getCoordRotRev(xRot,yRot,R_mat,xRotPt,yRotPt):
    
    # Get shape of input X,Y data
    xRotShape = np.shape(xRot)
    yRotShape = np.shape(yRot)
    
    # Convert data to columns
    xRot_mat = np.c_[xRot]
    yRot_mat = np.c_[yRot]
    
    # Get shape of matrices
    (xRot_matRows,xRot_matCols) = xRot_mat.shape
    (yRot_matRows,yRot_matCols) = yRot_mat.shape
    
    # Make X input a row vector
    if(xRot_matRows > 1):
        xRot_mat = np.transpose(xRot_mat)
    #endif
    
    # Make Y input a row vector
    if(yRot_matRows > 1):
        yRot_mat = np.transpose(yRot_mat)
    #endif
    
    # Put X,Y data into 2 x N matrix
    xyRot_mat = np.concatenate((xRot_mat,yRot_mat))

    # Rotate data back to original frame
    measUnrot_mat = np.matmul(np.linalg.inv(R_mat),xyRot_mat)
    
    # Pull out X,Y unrotated data
    xUnrot_mat = measUnrot_mat[0,:]
    yUnrot_mat = measUnrot_mat[1,:]
    
    # Translate data back to original point
    xOut_mat = xUnrot_mat + xRotPt
    yOut_mat = yUnrot_mat + yRotPt
    
    # Convert matrices to numpy arrays for output
    xOut = np.squeeze(np.asarray(xOut_mat))
    yOut = np.squeeze(np.asarray(yOut_mat))
    
    # Make X,Y output the same shape as X,Y input
    xOut = np.reshape(xOut,xRotShape)
    yOut = np.reshape(yOut,yRotShape)
    
    # Return output variables
    return xOut, yOut, R_mat, xRotPt, yRotPt


##### Function to add in geoid model
def getGeoidHeight(geoidData,atlTruthData):
    
    # Convert truth data from UTM to Lat/Lon
    x = atlTruthData.x
    y = atlTruthData.y
    zone = atlTruthData.zone
    hemi = atlTruthData.hemi
    latsIn, lonsIn = getUTM2LatLon(x,y,zone,hemi)
        
    # Interpolate to find geoidal heights
    f = interpolate.interp2d(geoidData.lons, geoidData.lats, geoidData.geoidalHeights, kind='linear')
    geoidalHeights = interpolate.dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], lonsIn, latsIn)[0]
    geoidalHeights = np.c_[geoidalHeights]
        
    # Add geoidal heights to find new ellipsoidal heights (HAE)
    atlTruthData.z = atlTruthData.z + geoidalHeights
    
    return atlTruthData


##### Function to grid point cloud data
def getRaster(x, y, z, resolution, method, fillValue = -999, time = []):

    if(method.lower() == 'min'):
        npOperation = np.min
    elif(method.lower() == 'max'):
        npOperation = np.max
    elif(method.lower() == 'mean'):
        npOperation = np.mean
    elif(method.lower() == 'numel'):
        npOperation = np.ndarray.size
    else:
        npOperation = np.mean
    # EndIf
            
    # Round all incoming X,Y data
    xRnd = (np.round(x/resolution)*resolution).astype(int)
    yRnd = (np.round(y/resolution)*resolution).astype(int)
    
    # Get min,max of rounded X,Y data
    xRndMin = xRnd.min()
    xRndMax = xRnd.max()
    yRndMin = yRnd.min()
    yRndMax = yRnd.max()
    
    # Get all possible grid combinations
    xAll = np.arange(xRndMin, xRndMax + resolution, resolution)
    yAll = np.arange(yRndMax, yRndMin - resolution, -resolution)
    
    # Get X,Y array of all pts
    xAllArray, yAllArray = np.meshgrid(xAll,yAll)
    xyAll = np.column_stack((xAllArray.flatten(), yAllArray.flatten()))
    
    # Populate X,Y raster data
    numRows = len(yAll);
    numCols = len(xAll);
    rasterDataX = xAllArray;
    rasterDataY = yAllArray;
    
    # Get unique incoming X,Y data
    uniqueCombos, uniqueGroups = np.unique((xRnd,yRnd), return_inverse = True, axis = 1)
    uniqueCombos = np.transpose(uniqueCombos)
    
    # Find index locations of unique incoming X,Y data in X,Y array of all pts
    _, indsToFill = ismember(uniqueCombos,xyAll,'rows')
    indsToFill = np.c_[indsToFill]
    
    # Grid Z data, populate into raster array, and reshape 
    df = pd.DataFrame(np.column_stack([z, uniqueGroups]), columns=['z', 'unique_groups'])
    zGroups = df.groupby('unique_groups')
    zOut = zGroups.aggregate(npOperation)
    zOut = np.array(zOut)
    zSplit = np.c_[zOut[:,0]] 
    zRaster = fillValue*np.ones((numRows*numCols,1))
    zRaster[indsToFill,0] = zSplit
    rasterDataZ = np.reshape(zRaster,(numRows,numCols))
    
    if(any(time)):
        
        df = pd.DataFrame(np.column_stack([time, uniqueGroups]), columns=['time', 'unique_groups'])
        tGroups = df.groupby('unique_groups')
        tOut = tGroups.aggregate(npOperation)
        tOut = np.array(tOut)
        tSplit = np.c_[tOut[:,0]] 
        tRaster = fillValue*np.ones((numRows*numCols,1))
        tRaster[indsToFill,0] = tSplit
        rasterDataT = np.reshape(tRaster,(numRows,numCols))
        
    else:
        
        rasterDataT = []
        
    # EndIf

    return GridStruct(rasterDataX, rasterDataY, rasterDataZ, rasterDataT)


##### Function to find closest points in an array
def getClosest(inputArray, closestPts):
    
    # Initialize outputs
    minInd = np.zeros(np.shape(closestPts), dtype = int)
    minVal = np.zeros(np.shape(closestPts))
    
    # Loop through closest points array and find closest point to input array
    for i in range(0,len(closestPts)):
        closestPt = closestPts[i]
        arrayDif = np.abs(inputArray - closestPt)
        minInd[i] = np.argmin(arrayDif) 
        minVal[i] = inputArray[minInd[i]]
    # EndFor
        
    # Return outputs
    return minVal, minInd

def __appendGlobalList(name):
    if name:
        global_list.append(name)

def getH5Keys(h5_file,group = None, out_txt = None, verbose = False):
    global global_list
    global_list = []
    try:
        h = h5py.File(h5_file, 'r')
    except:
        print("Could not find file or file was not proper H5 file")
        sys.exit
    if group:
        group = str(group)
        h[group].visit(__appendGlobalList)
    else:
        h.visit(__appendGlobalList)
    if verbose:
        print(*global_list, sep = "\n")
    if out_txt:
        with open(out_txt, 'w', newline = '') as csvFile:
            with open(out_txt, 'w') as f:
                for item in global_list:
                    f.write("%s\n" % item)
        csvFile.close
    return global_list