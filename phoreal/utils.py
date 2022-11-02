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

# Filter Runtime warnings
import warnings

# Import Python modules
import numpy as np
import sys
import os
from fiona._transform import _transform
from scipy import interpolate
import pandas as pd
import h5py
import ctypes
from numpy.ctypeslib import ndpointer 
import copy
import sliderule
# Import ICESat-2 modules
from phoreal.gui_addins import (superFilterFile_windows, superFilterFile_linux)

EPSG_ARCTIC = '3413'
EPSG_ANTARCTIC = '3976'
EPSG_ERR = '-1'
DEG_ARCTIC = 84.0
DEG_ANTARCTIC = -80.0

EPSG_CEA = '6933' # Lambert Cylindrical Equal Area
EPSG_LAEA_N = '6931' # Lambert Azimuthal Equal Area
EPSG_LAEA_S = '6932'

GT_NUMS = ['gt1r', 'gt2r', 'gt3r', 'gt1l', 'gt2l', 'gt3l']

# Object for getNameParts function
class fileStruct:
        
    # Define class with designated fields
    def __init__(self, atlVersion, year, month, day, hour, minute, second, 
                 trackNum, unknown, releaseNum, incrementNum):
            
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
    startPos = 0
    error = False
    for i in range(0,len(nameSplit)):
        if 'ATL' in nameSplit[i]:
            startPos = startPos + i
            error = False
            break
        else:
            error = True
    if error == False:
        # Get ATL version
        atlVersion = nameSplit[startPos + 0]
        
        # Get Year, month, day
        year = nameSplit[startPos + 1][0:4]
        month = nameSplit[startPos + 1][4:6]
        day = nameSplit[startPos + 1][6:8]
        
        # Get Hour, minute, second
        hour = nameSplit[startPos + 1][8:10]
        minute = nameSplit[startPos + 1][10:12]
        second = nameSplit[startPos + 1][12:14]
        
        # Get other details
        trackNum = nameSplit[startPos + 2][0:4]
        unknown = nameSplit[startPos + 2][4:9]
        releaseNum = nameSplit[startPos + 3]
        incrementNum = nameSplit[startPos + 4]
    else:
        print('     ATL Filename Not in Standard Format, Unable to record information')
        # Get ATL version
        atlVersion = "000"
        
        # Get Year, month, day
        year = "0000"
        month = "00"
        day = "00"
        
        # Get Hour, minute, second
        hour = "00"
        minute = "00"
        second = "00"
        
        # Get other details
        trackNum = "0000"
        unknown = "0000"
        releaseNum = "0000"
        incrementNum = "0000"
    
        
    # Get data into class structure
    fileInfo = fileStruct(atlVersion, year, month, day, hour, minute, second, trackNum, unknown, releaseNum, incrementNum)
    
    # Return class
    return fileInfo

# endDef

##### Function to represent 2 numbers as 1 unique number
def cantorPairing(arrayIn):
    
    # Do Cantor pairing algorithm to express two values as one unique value
    vectorOut = 0.5 * (arrayIn[:,0] + arrayIn[:,1]) * (arrayIn[:,0] + arrayIn[:,1] + 1) + arrayIn[:,1]
    
    return vectorOut

    
##### Function to determine members of one array in another
def ismember(a_vec, b_vec, method_type = 'normal'):
    
    """ MATLAB equivalent ismember function """
    
    # Combine multi column arrays into a 1-D array of strings if necessary
    # This will ensure unique rows when using np.isin below
    if(method_type.lower() == 'rows'):
        
        # Turn a_vec into an array of strings
        a_str = a_vec.astype('str')
        b_str = b_vec.astype('str')
        
        # Concatenate each column of strings with commas into a 1-D array of strings
        for i in range(0,np.shape(a_str)[1]):
            a_char = np.char.array(a_str[:,i])
            b_char = np.char.array(b_str[:,i])
            if(i==0):
                a_vec = a_char
                b_vec = b_char
            else:
                a_vec = a_vec + ',' + a_char
                b_vec = b_vec + ',' + b_char

    matchingTF = np.isin(a_vec,b_vec)
    common = a_vec[matchingTF]
    common_unique, common_inv  = np.unique(common, return_inverse=True)     # common = common_unique[common_inv]
    b_unique, b_ind = np.unique(b_vec, return_index=True)  # b_unique = b_vec[b_ind]
    common_ind = b_ind[np.isin(b_unique, common_unique, assume_unique=True)]
    matchingInds = common_ind[common_inv]
    
    return matchingTF, matchingInds


##### Function to determine intersection of two arrays
def getIntersection2d(a_vec, b_vec, assume_unique=False):
    
    # Get max total number of rows in a_vec and b_vec
    a_vec_maxRows = np.max(a_vec[:,1]) - np.min(a_vec[:,1]) + 1
    b_vec_maxRows = np.max(b_vec[:,1]) - np.min(b_vec[:,1]) + 1
    maxRows = np.max([a_vec_maxRows, b_vec_maxRows])
    
    # Convert x,y locations to index values
    a_vec_IDs = a_vec[:,0]*maxRows + a_vec[:,1]
    b_vec_IDs = b_vec[:,0]*maxRows + b_vec[:,1]
    
    # Get common index values
    commonIDs, a_inds, b_inds = np.intersect1d(a_vec_IDs, b_vec_IDs, assume_unique, return_indices = True)
    
    # Get common values
    commonVals = a_vec[a_inds]
    
    # Return output
    return commonVals, a_inds, b_inds


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
    atl03_ph_beg_val = atl03_ph_index_beg[atl03_ph_beg_inds]
    newMapping = atl08classed_inds + atl03_ph_beg_val - 2
    
    # Get max size of output array
    sizeOutput = newMapping[-1]
    
    # Pre-populate all photon classed array with zeroes
    allph_classed = (np.zeros(sizeOutput + 1)) - 1
    
    # Populate all photon classed array from ATL08 classifications
    allph_classed[newMapping] = atl08classed_vals
    
    # Return all photon classed array
    return allph_classed



##### Functions to convert Lat/Lon to UTM and vice versa
'''ATL Geographic Coordinate System Converter

ATL GCS Converter will translate an ATL groundtrack to different coordinate systems.
Mainly this is to convert Lonitute and latitude in WGS84 coordinates, which are
radial, to cartesian coordinates (e.g., UTM) for easy processing, and if required,
back to WGS84.  Primary package required is pyproj.

Each UTM grid is standardized across the globe but there a few exceptions. Zones
    31V, 32V, 31X, 33X, 35X, and 37X are irregularly sized.  Current code does 
    not make these exceptions.  This should not make too big of a difference for
    now but will become an issue if code needs to work with other GIS.
    For more information:
    https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system
    https://upload.wikimedia.org/wikipedia/commons/e/ed/Utm-zones.jpg
    https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#/ media/File:Modified_UTM_Zones.png

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
def find_utm_zone_point(lon, lat, module=None, m=None):

    """
    Input:
    longitude, latitude, {opt} module, {opt} m


    Output:
    epsg_code string


    Desc:
        longitude, latitude - in degrees
            No checks on lon/lat bounds

        module - options: [None, 'mgrs', 'pygeodesy']
            optional
            If nothing is given, it defaults to None
            If a wrong module is given, it defaults to None condition
            If mgrs or pygeodesy are given, it uses the 
            respective module

        m - mgrs object, i.e. m = mgrs.MGRS()
            optional input to allow user to predefine
            MGRS object since it doesn't have to be
            defined for every loop

            If it defined as anything other than
            mgrs.MGRS(), code will throw an error


    Latitude Bounds
    UTM exists between -80 and 84 deg latitude. UPS or
    NSIDC Sea Ice Polar Stereographic exists either
    south of -80 deg or north of 84 deg. The exact
    latitude bounds differ by method, as shown below
    in interval set notation:

    if module == None:
        antarctic:  [-90, -80]
        UTM:        (-80, 84)
        arctic:     [84,90]
        Does not take Norway/Svalbard UTM zone changes
        into account

    if module == 'mgrs' or module == 'pygeodesy':
        antarctic:  [-90, -80)
        UTM:        [-80, 84)
        arctic:     [84,90]
        Includes all one-offs in UTM zones
        

    Precision
        None has a precision of round-off
        mgrs has a precision of 1m
        pygeodesy has a precision of round-off

    """


    if not (module == 'mgrs' or module == 'pygeodesy'): # None condition

        utm_band = str(int((np.floor((lon + 180) / 6 ) % 60) + 1))
        if len(utm_band) == 1:
            utm_band = '0'+utm_band
            # print(utm_band)
        if lat >= 0 and lat < 84:
            epsg_code = '326' + utm_band
        elif lat >= 84:
            # epsg_code = '3411'
            epsg_code = EPSG_ARCTIC
        elif lat <= -80:
            # epsg_code = '3412'
            epsg_code = EPSG_ANTARCTIC
        else:
            epsg_code = '327' + utm_band
        return epsg_code


    elif module == 'mgrs':
        import mgrs
        if m == None:
            m = mgrs.MGRS()

        # UTM defined in [-80, 84)
        if DEG_ANTARCTIC <= lat < DEG_ARCTIC:
            # [-80,84)
            mgrs_code = m.toMGRS(lat, lon, MGRSPrecision=5) # 1m precision
            UTM_code = m.MGRSToUTM(mgrs_code)
            utm_band = str(UTM_code[0]).zfill(2)
            hemi = UTM_code[1].decode() # 'N' or 'S'

            if not (hemi == 'N' or hemi == 'S'):
                # debug check since code is new
                print('warning ({}): hemi={} is something other than N or S'.format(module, hemi))

            epsg_code = '326' + utm_band # N
            if hemi == 'S':
                epsg_code = '327' + utm_band

        elif lat < DEG_ANTARCTIC:
            # [-90, -80), polar antarctic
            epsg_code = EPSG_ANTARCTIC

        elif lat >= DEG_ARCTIC:
            # [84, 90], polar arctic
            epsg_code = EPSG_ARCTIC

        return epsg_code


    elif module == 'pygeodesy':
        import pygeodesy as pg

        if DEG_ANTARCTIC <= lat < DEG_ARCTIC:
             # [-80,84)
            z = pg.utm.utmZoneBand5(lat, lon)
            utm_band = str(z.zone).zfill(2)
            hemi = z.hemipole # 'N' or 'S'

            if not (hemi == 'N' or hemi == 'S'):
                # debug check since code is new
                print('warning ({}): hemi={} is something other than N or S'.format(module, hemi))

            epsg_code = '326' + utm_band # N
            if hemi == 'S':
                epsg_code = '327' + utm_band


        elif lat < DEG_ANTARCTIC:
            # [-90, -80), polar antarctic
            epsg_code = EPSG_ANTARCTIC

        elif lat >= DEG_ARCTIC:
            # [84, 90], polar arctic
            epsg_code = EPSG_ARCTIC

        return epsg_code




# Find UTM zone for numpy array of lon/lat.
def find_utm_zone_arr(lat, lon, mode=True, module=None, m=None):

    """
    Input:
    longitude, latitude, {opt} mode, {opt} module, {opt} m
    

    Output:
    if mode:
        return a single epsg_code
    else:
        return a list of epsg_codes relative to lon/lat


    Desc:
        longitude, latitude - in degrees
            No checks on lon/lat bounds
            These must be the same length

        mode - whether to return the mode of zones or not
            default is True

        module - options: [None, 'mgrs', 'pygeodesy']
            optional
            If nothing is given, it defaults to None
            If a wrong module is given, it defaults to None condition
            If mgrs or pygeodesy are given, it uses the 
            respective module

        m - mgrs object, i.e. m = mgrs.MGRS()
            optional input to allow user to predefine
            MGRS object since it doesn't have to be
            defined for every loop

            If it defined as anything other than
            mgrs.MGRS(), code will throw an error

    See find_utm_zone_point() for additional details.

    """
    default_mod_bool = not (module == 'mgrs' or module == 'pygeodesy')

    if type(lon) != np.ndarray:
        lon = np.array(lon)
    if type(lat) != np.ndarray:
        lat = np.array(lat)

    if lon.size == 1:
        lon = np.array([lon])
    if lat.size == 1:
        lat = np.array([lat])

    if len(lon) != len(lat):
        print('error: len(lon) != len(lat)')
        return EPSG_ERR
        
    if default_mod_bool and mode:
        # OG way to find zone via the mode of all points

        arctic = len(lat[lat >= 84.0])
        nhem = len(lat[(lat < 84.0) & (lat >= 0)])
        shem = len(lat[(lat > -80.0) & (lat < 0)])
        antarctic = len(lat[lat < -80.0])
        if arctic > nhem: #Arctic case
            epsg_code = EPSG_ARCTIC # '3413'
        elif antarctic > shem: #Antarctic case
            epsg_code = EPSG_ANTARCTIC #'3976'
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
            elif len(zone) > 3:
                from scipy.stats import mode as mode_func
                zone = mode_func(zone)[0][0]
                print("Warning: Input ground track present in more than 3 UTM zones. \
                         \nRecommend manually selecting GCS.")
            else:
                # len(zone) == 0
                print("Warning: zone == [], lon/lat may not have values")
                sys.exit()

            if nhem >= shem:
                if len(str(zone)) == 1:
                    zone = "0" + str(zone)
                epsg_code = '326' + str(zone)
            else:
                if len(str(zone)) == 1:
                    zone = "0" + str(zone)
                epsg_code = '327' + str(zone)
        return epsg_code



    else:
        # OG method and mode or
        # OG method and not mode or
        # different method and mode or
        # different method and not mode

        epsg_code_arr = []
        for i in range(len(lon)):
            out = find_utm_zone_point(lon[i], lat[i], module=module, m=m)
            epsg_code_arr.append(out)

        if mode:
            from scipy.stats import mode as mode_func
            epsg_code = mode_func(epsg_code_arr)[0][0]
            return epsg_code

        else:
            return epsg_code_arr

# Transform GCS/PCS based on EPSG and x/y. 
def transform(epsg_in, epsg_out, x, y, use_old_version=False):
    
    # Using fiona
    xx_list, yy_list = _transform(epsg_in, epsg_out, x, y)
    xx = np.array(xx_list)
    yy = np.array(yy_list)
        
    return xx,yy

# Transform from lon/lat to given EPSG. 
def wgs84_to_epsg_transform(epsg_out, lon, lat):
    epsg_in = 'epsg:4326'
    epsg_out = ('epsg:{0}'.format(str(epsg_out)))
    xx, yy = transform(epsg_in, epsg_out, lon, lat)
    return xx,yy

def identifyEPSG(hemi,zone):
    if hemi == 'N':
        outstring = 'epsg:326'
    elif hemi == "S":
        outstring = 'epsg:327'
    else:
        outstring = ''
        print('error')
    
    outstring = outstring + (str(zone).zfill(2))
    
    return outstring

def identify_hemi_zone(epsg):
    epsg = str(epsg)
    epsg = epsg.split(':')[-1]
    
    if (epsg[0:3] == '326'):
        hemi = 'N'
        zone = epsg[3:5]
    elif (epsg[0:3] == '327'):
        hemi = 'N'
        zone = epsg[3:5]
    elif(epsg =='3413'):
        zone = epsg 
        hemi = 'arctic'
    elif(epsg =='3976'):
        zone = epsg 
        hemi = 'arctic'
    else:
        print('Could not read EPSG for hemi/zone')
        zone = epsg 
        hemi = ''
    return hemi, zone
        

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
            if len(zone) == 1:
                zone = '0' + zone
            epsg_out = 'epsg:326' + zone
        else:
            if len(zone) == 1:
                zone = '0' + zone 
            epsg_out = 'epsg:327' + zone
        # endif
        
        # Call transform function
        if len(lon) > 0 and len(lat) > 0:
            print(epsg_in)
            print(epsg_out)
            xx, yy = transform(epsg_in, epsg_out, lat, lon)
        else:
            xx, yy = np.array([]), np.array([])

    else:
        
        # Get UTM coords
        xx, yy, epsg_out = wgs84_to_utm_find_and_transform(lat, lon)
        
        if(epsg_out=='3413'):
            
            zone = epsg_out 
            hemi = 'arctic'
            
        elif(epsg_out=='3976'):
            
            zone = epsg_out 
            hemi = 'antarctic'
        
        else:
        
            # Store zone
            zone = epsg_out[3:]
            
            # Store hemisphere
            if(epsg_out[0:3]=='326'):   
                hemi = 'N'  
            else:
                hemi = 'S'
    
    # Return output
    return xx, yy, zone, hemi


# Function to convert UTM to lat/lon
def getUTM2LatLon(x,y,zone,hemi):
    
    # Set EPSG code for lat/lon coords
    epsg_out = 'epsg:4326'
    
    # Get EPSG code for UTM coords
    if(hemi=='N'):
        if len(zone) == 1:
            zone = "0" + zone
        epsg_in = 'epsg:326' + zone
    else:
        if len(zone) == 1:
            zone = "0" + zone
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
    
    # Make Y input a row vector
    if(yTranslated_matRows > 1):
        yTranslated_mat = np.transpose(yTranslated_mat)
    
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
    
    # Make Y input a row vector
    if(yRot_matRows > 1):
        yRot_mat = np.transpose(yRot_mat)
    
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
    x = atlTruthData.easting
    y = atlTruthData.northing
    zone = atlTruthData.zone
    hemi = atlTruthData.hemi
    latsIn, lonsIn = getUTM2LatLon(x,y,zone,hemi)
        
    # Interpolate to find geoidal heights
    lons = geoidData.lons[0]
    # lats = geoidData.lats[0]
    # geoidalHeights = geoidData.geoidalHeights[0]
    lons[lons > 180.] = lons[lons > 180] - 360
    geoidData.lons[0] = lons
    f = interpolate.interp2d(geoidData.lons, geoidData.lats, geoidData.geoidalHeights, kind='linear')
    geoidalHeights = interpolate.dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], lonsIn, latsIn)[0]
    geoidalHeights = np.c_[geoidalHeights]
        
    # Add geoidal heights to find new ellipsoidal heights (HAE)
    atlTruthData.z = atlTruthData.z + geoidalHeights
    
    return atlTruthData


##### Function to add in geoid model
def getGeoidHeight2(geoidData,x,y,z_msl,zone,hemi):
    
    # Convert truth data from UTM to Lat/Lon
    latsIn, lonsIn = getUTM2LatLon(x,y,zone,hemi)
        
    # Interpolate to find geoidal heights
    lons = geoidData.lons[0]
    # lats = geoidData.lats[0]
    # geoidalHeights = geoidData.geoidalHeights[0]
    lons[lons > 180.] = lons[lons > 180] - 360
    geoidData.lons[0] = lons
    f = interpolate.interp2d(geoidData.lons, geoidData.lats, geoidData.geoidalHeights, kind='linear')
    geoidalHeights = interpolate.dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], lonsIn, latsIn)[0]
    # geoidalHeights = np.c_[geoidalHeights]
        
    # Add geoidal heights to find new ellipsoidal heights (HAE)
    z_hae = z_msl + geoidalHeights
    
    return z_hae


##### Function to grid point cloud data
def getRaster(x, y, z, resolution, method, fillValue = -999, time = [], xAllArray = [], yAllArray = [],
                origin=None):
    
    # USER INPUTS
    # ---------------------------
    # x = input x array of values
    # y = input y array of values
    # z = input z array of values
    # resolution = resolution of grid cells (N = N x N, [M,N] = M x N)
    # method = operation to perform in each grid cell
    #   - min
    #   - max
    #   - mean (default)
    #   - median
    #   - mode
    #   - range
    #   - std (standard deviation)
    #   - numel (number of elements)
    # fillValue = value to fill in empty grid cells
    # time = secondary array (like z) to perform operation on in each grid cell
    # xAllArray = force output X grid cells (input like this: [min_x, max_x])
    # yAllArray = force output Y grid cells (input like this: [min_y, max_y])
    #
    # OUTPUTS
    # -------------
    # output.x = x (2D raster)
    # output.y = y (2D raster)
    # output.grid = grid (2D raster)
    # output.t = time (2D raster)
    
    
    # Get X,Y resolution
    if isinstance(resolution,np.integer):
        xResolution = float(resolution)
        yResolution = float(resolution)
    elif isinstance(resolution,int):
        xResolution = float(resolution)
        yResolution = float(resolution)
    elif isinstance(resolution,float):
        xResolution = resolution
        yResolution = resolution
    elif isinstance(resolution,np.ndarray):
        if len(resolution) == 1:
            xResolution = float(resolution)
            yResolution = float(resolution)
        elif len(resolution) == 2:
            xResolution = float(resolution[0])
            yResolution = float(resolution[1])
        else:
            print("Incorrect resolution input")
    elif isinstance(resolution,list):
        xResolution = float(resolution[0])
        yResolution = float(resolution[1])
    elif isinstance(resolution,str):
        strList = resolution.split(",")
        if len(strList) == 1:
            xResolution = float(resolution)
            yResolution = float(resolution)
        elif len(strList) == 2:
            xResolution = float(strList[0])
            yResolution = float(strList[1])
        else:
            print("Incorrect resolution input")
    else:
        print("Incorrect resolution input")

    # Get grid method
    if(method.lower() == 'min'):
        npOperation = np.nanmin
    elif(method.lower() == 'max'):
        npOperation = np.nanmax
    elif(method.lower() == 'mean'):
        npOperation = np.nanmean
    elif(method.lower() == 'median'):
        npOperation = np.nanmedian
    elif(method.lower() == 'mode'):
        npOperation = mode
    elif(method.lower() == 'range'):
        npOperation = np.range
    elif(method.lower() == 'std'):
        npOperation = np.nanstd
    elif(method.lower() == 'numel'):
        npOperation = np.size
    else:
        npOperation = np.mean
            
    # Round all incoming X,Y data
    xRnd = (np.round(x/xResolution)*xResolution)
    yRnd = (np.round(y/yResolution)*yResolution)
    
    # Get output X,Y grid cells
    if(any(xAllArray) and any(yAllArray)):
        
        manual_xy = True
        
        # Get xy limits
        xRndMin = (np.round(xAllArray[0]/xResolution)*xResolution)
        xRndMax = (np.round(xAllArray[1]/xResolution)*xResolution)
        yRndMin = (np.round(yAllArray[0]/yResolution)*yResolution)
        yRndMax = (np.round(yAllArray[1]/yResolution)*yResolution)
        
        xInside = np.logical_and(xRnd>=xRndMin, xRnd<=xRndMax)
        yInside = np.logical_and(yRnd>=yRndMin, yRnd<=yRndMax)
        indsInside = np.logical_and(xInside, yInside)
        xRnd = xRnd[indsInside]
        yRnd = yRnd[indsInside]
        z = z[indsInside]
        
    else:
        
        manual_xy = False
        
        # Get min,max of rounded X,Y data
        xRndMin = xRnd.min()
        xRndMax = xRnd.max()
        yRndMin = yRnd.min()
        yRndMax = yRnd.max()
        
    # Get all possible grid combinations
    xAll = np.arange(xRndMin, xRndMax + xResolution, xResolution)
    yAll = np.arange(yRndMax, yRndMin - yResolution, -yResolution)
        
    # Get X,Y array of all pts
    xAllArray, yAllArray = np.meshgrid(xAll,yAll)
    
    # Get raster X, Y, Z space
    rasterDataX = xAllArray.astype('float')
    rasterDataY = yAllArray.astype('float')
    rasterDataZ = fillValue*np.ones((np.shape(rasterDataX))).astype('float')
    
    # Make sure input data is of type float
    if(xRnd.dtype!='float64'):
        xRnd = xRnd.astype('float')
    if(yRnd.dtype!='float64'):
        yRnd = yRnd.astype('float')
    if(z.dtype!='float64'):
        z = z.astype('float')
    
    # Get x-rastered, y-rastered, and z data into array
    data = np.column_stack([xRnd, yRnd, z])
    
    # Put array into Pandas dataframe
    df = pd.DataFrame(data, columns=['xRnd', 'yRnd', 'z'])
    
    # Do groupby.agg to get get rastered operation for group
    groupedData = df.groupby(['xRnd', 'yRnd']).agg({'z': [npOperation]})
    groupedData.columns = ['z_agg']
    groupedData = groupedData.reset_index()
    zValsNew = np.array(groupedData['z_agg'])
    
    # Determine new row, column indices to place rastered Z data into
    df_xRnd_min = np.min(rasterDataX)
    df_yRnd_min = np.min(rasterDataY)
    colIndsNew = ((np.array(groupedData['xRnd']) - df_xRnd_min)/xResolution).astype(int)
    rowIndsNew = ((np.array(groupedData['yRnd']) - df_yRnd_min)/yResolution).astype(int)

    # Populate rastered Z data into array
    rasterDataZ[rowIndsNew, colIndsNew] = zValsNew
    rasterDataZ = np.flipud(rasterDataZ)
    
    # Grid 'time' array if necessary
    # if(any(time)):
        
    #     if(manual_xy):
    #         time = time[indsInside]
        
    #     # Get x-rastered, y-rastered, and time data into array
    #     if(time.dtype!='float64'):
    #         time = time.astype('float')
    #     dataTime = np.column_stack([xRnd, yRnd, time])
        
    #     # Put array into Pandas dataframe
    #     dfTime = pd.DataFrame(dataTime, columns=['xRnd', 'yRnd', 'time'])
        
    #     # Do groupby.agg to get get rastered operation for group
    #     groupedDataTime = dfTime.groupby(['xRnd', 'yRnd']).agg({'time': [npOperation]})
    #     groupedDataTime.columns = ['time_agg']
    #     groupedDataTime = groupedDataTime.reset_index()
    #     tValsNew = np.array(groupedDataTime['time_agg'])
        
    #     # Populate rastered Z data into array
    #     rasterDataT = fillValue*np.ones((np.shape(rasterDataX)))
    #     rasterDataT[rowIndsNew, colIndsNew] = tValsNew
    #     rasterDataT = np.flipud(rasterDataT)
        
    # else:
    rasterDataT = []
        
    return GridStruct(rasterDataX, rasterDataY, rasterDataZ, rasterDataT)

# endDef


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

    return minVal, minInd

def __appendGlobalList(name):
    if name:
        global_list.append(name)

def getH5Keys(h5_file,group = None, out_txt = None, verbose = False, matchText = None):
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
    if matchText:
        global_list = [s for s in global_list if matchText in s]
        
    return global_list

def sortAtlMeasured(atlMeasuredData, verbose=False):
    
    if(verbose):
        print("Sorting ATL Measured Data...", end = " ")
    sort_index = np.argsort(atlMeasuredData.alongTrack, axis = 0)
    atlMeasuredData.alongTrack = atlMeasuredData.alongTrack[sort_index[:,0]]
    atlMeasuredData.crossTrack = atlMeasuredData.crossTrack[sort_index[:,0]]
    atlMeasuredData.lat = atlMeasuredData.lat[sort_index[:,0]]
    atlMeasuredData.lon = atlMeasuredData.lon[sort_index[:,0]]
    atlMeasuredData.northing = atlMeasuredData.northing[sort_index[:,0]]
    atlMeasuredData.easting = atlMeasuredData.easting[sort_index[:,0]]
    atlMeasuredData.z = atlMeasuredData.z[sort_index[:,0]]
    atlMeasuredData.time = atlMeasuredData.time[sort_index[:,0]]
    atlMeasuredData.classification = \
    atlMeasuredData.classification[sort_index[:,0]]
    atlMeasuredData.signalConf = atlMeasuredData.signalConf[sort_index[:,0]]
   
    if(verbose):
        print("Complete")
    return atlMeasuredData
    
def sortAtlTruth(atlTruthData, verbose=False):
    if(verbose):
        print("Sorting ATL Truth Data...", end = " ")  
    sort_index = np.argsort(atlTruthData.alongTrack, axis = 0)
    atlTruthData.alongTrack = atlTruthData.alongTrack[sort_index[:,0]]
    atlTruthData.crossTrack = atlTruthData.crossTrack[sort_index[:,0]]
    atlTruthData.northing = atlTruthData.northing[sort_index[:,0]]
    atlTruthData.easting = atlTruthData.easting[sort_index[:,0]]
    atlTruthData.lat = atlTruthData.lat[sort_index[:,0]]
    atlTruthData.lon = atlTruthData.lon[sort_index[:,0]]
    atlTruthData.z = atlTruthData.z[sort_index[:,0]]
    atlTruthData.intensity = atlTruthData.intensity[sort_index[:,0]]
    atlTruthData.classification = atlTruthData.classification[sort_index[:,0]]

    if(verbose):
        print("Complete")
    return atlTruthData

def indexMatch(measuredArray,truthArray,verbose=False):
    if(verbose):
        print("Match corresponding indices...", end = " ")
    A = np.array(measuredArray)
    B = np.array(truthArray)
    C = np.empty((len(B)))
    
    
    if os.name == 'nt':
        lib = ctypes.cdll.LoadLibrary(os.path.abspath(superFilterFile_windows))
    else:
        lib = ctypes.cdll.LoadLibrary(os.path.abspath(superFilterFile_linux))
    fun = lib.cfun
    fun.restype = None
    fun.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
            ctypes.c_size_t,
            ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
            ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
            ctypes.c_size_t]
    fun(A, A.size, B, C, B.size)
    if(verbose):
        print("Complete")
    return np.array(C).astype(int)

    
def superFilter(atlMeasuredData_in, atlTruthData_in, xBuf = 7, classCode = [], verbose=False):
    if(verbose):
        print("Applying Superfilter")
    #Sort Measured
    atlMeasuredData = copy.deepcopy(atlMeasuredData_in)    
    atlTruthData = copy.deepcopy(atlTruthData_in)    

    atlMeasuredData = sortAtlMeasured(atlMeasuredData)
    #Sort Truth
    atlTruthData = sortAtlTruth(atlTruthData)
    
    
    #Find Matching Indexes
    if classCode:
        filter_class = np.isin(atlMeasuredData.classification,classCode)
        
        alongTrack = atlMeasuredData.alongTrack[filter_class]
        crossTrack = atlMeasuredData.crossTrack[filter_class]
        indexMatches = indexMatch(alongTrack,atlTruthData.alongTrack)
        #Generate filter
        indexMatches[indexMatches >= len(crossTrack)] = (len(crossTrack) - 1)
        x_check = crossTrack[indexMatches]
        x_diff = atlTruthData.crossTrack[:,0] - x_check
        filter_data = np.where((x_diff < xBuf) & (x_diff > -xBuf))    
        #Apply filter to truth data
        atlTruthData.alongTrack = atlTruthData.alongTrack[filter_data]
        atlTruthData.crossTrack = atlTruthData.crossTrack[filter_data]
        atlTruthData.northing = atlTruthData.northing[filter_data]
        atlTruthData.easting = atlTruthData.easting[filter_data]
        atlTruthData.lat = atlTruthData.lat[filter_data]
        atlTruthData.lon = atlTruthData.lon[filter_data]
        atlTruthData.z = atlTruthData.z[filter_data]
        atlTruthData.intensity = atlTruthData.intensity[filter_data]
        atlTruthData.classification = atlTruthData.classification[filter_data]
        atlTruthData.year = atlTruthData.year[filter_data]
        atlTruthData.month = atlTruthData.month[filter_data]
        atlTruthData.day = atlTruthData.day[filter_data]
        atlTruthData.time = atlTruthData.time[filter_data]
        atlTruthData.deltaTime = atlTruthData.deltaTime[filter_data]
    else:
        indexMatches = indexMatch(atlMeasuredData.alongTrack, \
                                  atlTruthData.alongTrack)
        #Generate filter
        indexMatches[indexMatches >= len(atlMeasuredData.crossTrack)] = \
            len(atlMeasuredData.crossTrack) - 1
        x_check = atlMeasuredData.crossTrack[indexMatches]
        x_diff = atlTruthData.crossTrack - x_check
        filter_data = np.where((x_diff < xBuf) & (x_diff > -xBuf))    
        #Apply filter to truth data
        atlTruthData.alongTrack = atlTruthData.alongTrack[filter_data]
        atlTruthData.crossTrack = atlTruthData.crossTrack[filter_data]
        atlTruthData.northing = atlTruthData.northing[filter_data]
        atlTruthData.easting = atlTruthData.easting[filter_data]
        atlTruthData.lat = atlTruthData.lat[filter_data]
        atlTruthData.lon = atlTruthData.lon[filter_data]
        atlTruthData.z = atlTruthData.z[filter_data]
        atlTruthData.intensity = atlTruthData.intensity[filter_data]
        atlTruthData.classification = atlTruthData.classification[filter_data] 
        atlTruthData.year = atlTruthData.year[filter_data]
        atlTruthData.month = atlTruthData.month[filter_data]
        atlTruthData.day = atlTruthData.day[filter_data]
        atlTruthData.time = atlTruthData.time[filter_data] 
        atlTruthData.deltaTime = atlTruthData.deltaTime[filter_data] 
    if(verbose):
        print("Superfilter complete")

    return atlTruthData, atlMeasuredData

# Moving mean function
def getMovingMean(x, nBack, nForward):
    
    # Get length of input array
    N = len(x)
    
    # Get starting/ending input array values
    firstVal = x[0]
    lastVal = x[-1]
            
    # Get window size
    windowSize = nBack + nForward + 1
    
    # Get starting/ending arrays to concatenate onto input array
    xStart = firstVal*np.ones((nBack + 1), dtype = 'int')
    xEnd = lastVal*np.ones(nForward, dtype = 'int')
    
    # Concatenate starting/ending arrays onto input array
    xNew = np.concatenate((xStart, x, xEnd))
    
    # Get cumulative sum of input array
    cumsum = np.cumsum(xNew)
    
    # Get back/forward sums
    backSum = cumsum[:N]
    forwardSum = cumsum[-N:]
    
    # Compute moving average
    movingMean = (forwardSum - backSum)/windowSize
    
    return movingMean

    
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

def pause(msg='enter to continue'):
    input(msg)

def mode(arr):
    from scipy.stats import mode as mode_func
    val = mode_func(arr)[0][0]
    return val

def sort_first(arr):
    """
    example: sort by x
        x_sort, y_sort, z_sort = cm.sort_first([x, y, z])
    """
    s = np.transpose(sorted(np.transpose(arr), key=lambda x: x[0]))
    return s

def remove_numpy_printing():
    np.set_printoptions(suppress=True)

def array(*argv):
    """
    Turns this:
      x, y, z = np.array(x), np.array(y), np.array(z)
    into
      import icesatUtils as it
      x, y, z = it.array(x, y, z)
    """
    var_out = []
    for var in argv:
        var_out.append(np.array(var))
    return var_out

def meshgrid_to_xy(xx, yy, zz=None):
    """
    Example:
    import icesatUtils
    import numpy as np
    x, y = np.linspace(0,1), np.linspace(0,1)
    xx, yy = np.meshgrid(x, y)
    xr, yr = icesatUtils.meshgrid_to_xy(xx, yy)
    """
    if type(zz) == type(None):
        return xx.flatten(), yy.flatten()
    else:
        return xx.flatten(), yy.flatten(), zz.flatten()


def get_date(*args, debug=0):
    
    """
    Written by Michael James

    Example:
    import icesatUtils
    doy = icesatUtils.get_date(year, month, day)
    # or
    month, day = icesatUtils.get_date(year, doy)

    """
    def help():
        print("""
    Example:
    import icesatUtils
    doy = icesatUtils.get_date(year, month, day)
    # or
    month, day = icesatUtils.get_date(year, doy)
            """)

    import datetime

    if len(args) == 2:
        y = int(args[0]) #int(sys.argv[1])
        d = int(args[1]) #int(sys.argv[2])
        yp1 = datetime.datetime(y+1, 1, 1)
        y_last_day = yp1 - datetime.timedelta(days=1)
        if d <= y_last_day.timetuple().tm_yday:
            date = datetime.datetime(y, 1, 1) + datetime.timedelta(days=d-1)
            return date.month, date.day
        else:
            print("error")
            help()

    elif len(args) == 3:
        y, m, d = args
        date = datetime.datetime(y, m, d)
        doy = int(date.timetuple().tm_yday)
        # print("doy = {}".format(date.timetuple().tm_yday))
        return str(doy).zfill(3)

    else:
        print("error: incorrect number of args")
        help()



def get_h5_meta(h5_file, meta='date', rtn_doy=False, rtn_hms=True, file_start='ATL'): #, debug=0):
    """
    This function gets metadata directly from the ATL filename.

    Input:
        h5_file - the ATL file, full-path or not
        meta - the type of metadata to output
            ['date', 'track', 'release', 'version', 
                'f_type', 'hms', 'cycle']
            f_type - either rapid or final
        rtn_doy - return day of year or not, if date is chosen
        rtn_hms - return hour/min/sec, or time in sec
        file_start - search for info based on file_start index;
                        useful if given an ATL file that starts
                        with "errorprocessing_..." or any other
                        prefix
        debug - for small errors

    Output:
        varies, but generally it's one value, unless 'hms' meta is chosen,
        in which case it is two values.

    Example:
        import icesatUtils
        fn = DIR + '/ATL03_20181016000635_02650109_200_01.h5'
        # or fn = 'ATL03_20181016000635_02650109_200_01.h5'
        year, day_of_year = icesatUtils.get_h5_meta(fn, meta='date', rtn_doy=True)
        version = icesatUtils.get_h5_meta(fn, meta='version')
        release = icesatUtils.get_h5_meta(fn, meta='release')
    """

    h5_file = os.path.basename(h5_file) # h5_file.split('/')[-1]

    meta = meta.lower()

    i0 = 0
    try:
        i0 = h5_file.index(file_start)
    except ValueError:
        print('warning: substring %s not found in %s' % (file_start, h5_file))
    # i0 = check_file_start(h5_file, file_start, debug)

    if meta == 'date':
        year = int(h5_file[i0+6:i0+10])
        month = int(h5_file[i0+10:i0+12])
        day = int(h5_file[i0+12:i0+14])

        if rtn_doy:
            doy0 = get_date(year, month, day)
            return str(year), str(doy0).zfill(3)

        return str(year), str(month).zfill(2), str(day).zfill(2)

    elif meta == 'track':
        return int(h5_file[i0+21:i0+25])

    elif meta == 'release':
        r = h5_file[i0+30:i0+34]
        if '_' in r:
            r = h5_file[i0+30:i0+33]
        return r

    elif meta == 'version':
        v = h5_file[i0+34:i0+36]
        if '_' in v:
          v = h5_file[i0+35:i0+37]
        return v

    elif meta == 'f_type':
        r = h5_file[i0+30:i0+34]
        f_type = 'rapid'
        if '_' in r:
            r = h5_file[i0+30:i0+33]
            f_type = 'final'
        return f_type

    elif meta == 'hms':
        hms = h5_file[i0+14:i0+20]
        h, m, s = hms[0:2], hms[2:4], hms[4:6]
        if rtn_hms:
            return h, m, s
        else:
            h, m, s = int(h), int(m), int(s)
            t0, t0_full = t2t(h,m,s)
            return t0, t0_full

    elif meta == 'cycle':
        return int(h5_file[i0+25:i0+29])

    else:
        print('error: unknown meta=%s' % meta)
        return 0

def get_index(t, t0, tol=2.0, err_msg=True):
    """
    Gets closest index of t0 in array t, i.e. the index of time t0
    in time array t

    tol is used in case t0 is not found; default is tol = 2.0 sec
    """
    if type(t) != np.ndarray:
        t = np.array(t)
    index = np.argmin(abs(t - t0))
    if abs(t[index] - t0) > tol:
        if err_msg:
            print('warning: index not found: looking for %4.2f vs. found %4.2f' % (t0, t[index]))
        # return -1
    return index


def t2t(hours, mins=None, sec=None):
    """
    This function changes hms into seconds, or 
    HH:MM:SSSSSS format, such as in GPS time-stamps.
    """
    if mins==None and sec== None:
        tot_sec = hours%86400
        choice = 1
    elif hours != None and mins != None and sec != None:
        choice = 2
        if hours > 24 or mins > 60 or sec > 60 or hours < 0 or mins < 0 or sec < 0:
            raise ValueError("Value out of range")
    if choice == 1:
        hours = int(tot_sec/3600)
        mins = int((tot_sec%3600)/(60))
        sec = ((tot_sec%3600)%(60))

    if choice == 2:
        tot_sec = hours*3600 + mins*60 + sec

    return [tot_sec, "{:>02d}:{:>02d}:{:>09.6f}".format(hours, mins, sec)]


def get_outlier_data(data):
    """
    This function returns simple Q1,Q2,Q3,IQR. It uses
    nanpercentile and nanmedian, so nans are removed/handled.
    """
    Q1,Q2,Q3 = np.nanpercentile(data,25), np.nanmedian(data), np.nanpercentile(data,75)
    IQR = Q3-Q1
    return Q1,Q2,Q3,IQR

def get_root_dir(dir_type=None, username=None, debug=0):

    """
    This function outputs directories relative to GLAM operatives.

    Input:
        dir_type - either 'data' or 'user', optional
        username - the function will attempt to detect the username
                    of the person that called it, but if you'd like
                    to manually specify, the option is there
        debug - optionally prints out the username the function detects

    Output:
        if dir_type is None:
            returns [data, user] root directories, relative to user
        elif dir_type is 'data':
            returns only data
        elif dir_type is 'user':
            returns only user

    Example:
        import icesatUtils as iu
        USER_ROOT, DATA_ROOT = iu.get_root_dir()
        USER_JSIPPS = iu.get_root_dir(dir_type='user', username='jsipps')

    """

    dirs = {'user': ['/foo/bar/']}

    # get or format username
    import pwd
    if type(username) == type(None):
        username = pwd.getpwuid(os.getuid()).pw_name
        if debug:
            print('username: %s' % username)
    username = username.lower()

    # get or format directory type
    #   user or data
    dir_types = {'user': 0, 'data': 1}
    if type(dir_type) != type(None):
        if type(dir_type) != str:
            dir_type = str(dir_type)
        dir_type = dir_type.lower()
        if dir_type in dir_types:
            i = dir_types[dir_type]
        else:
            keys = [key for key in dir_types.keys()]
            raise KeyError(dir_type + ' directory type unknown. Try one of', keys)

    # return directory
    if username in dirs:
        if type(dir_type) == type(None):
            return dirs[username]
        else:
            return dirs[username][i]

    else:
        keys = [key for key in dirs.keys()]
        raise KeyError('username ' + username + ' not found in', keys)

def reload(module):
    import imp
    imp.reload(module)


def get_sc_orient(file, delta_time=None):
    """
    Takes in any h5 file, but 03 and 08 should have
    the orbit_info group. If not, this function raises
    a KeyError.

    Outputs an array of sc_orient values, dependent on
    delta_time, unless delta_time is not specified, in
    which case it just outputs the sc_orient and sc_orient_time
    datasets, which are typically one element in length.
    """

    INT_MAX = np.iinfo(int).max
    group = 'orbit_info'
    datasets = [group + '/sc_orient', group + '/sc_orient_time']

    found = False
    with h5py.File(file, 'r') as fp:
        if datasets[0] in fp and datasets[1] in fp:
            sc_orient = np.array(fp[datasets[0]]).astype(int)
            sc_orient_time = np.array(fp[datasets[1]])
            found = True

    if found:
        if type(delta_time) != type(None):
            n = len(delta_time)
            sc_orient_arr = np.full(n, INT_MAX)
            for k, t0 in enumerate(sc_orient_time):
                b = delta_time >= t0
                sc_orient_arr[b] = sc_orient[k]

            return sc_orient_arr

        else:
            return sc_orient, sc_orient_time

    else:
        raise KeyError('could not find [%s or %s]' % (datasets[0], datasets[1]))


def get_beam_info(sc_orient_arr, gt):
    """
    Takes as input the sc_orient at every sample
    and a ground track; outputs what beam number
    and beam type should be at these samples.

    beam_number of -1 is assigned if sc_orient 
    value is unknown, such as 2 (during yaw-flip)
    and any other value.

    beam_type is set to unknown is sc_orient
    value is unknown.
    """

    n = len(sc_orient_arr)
    s0 = sc_orient_arr == 0
    s1 = sc_orient_arr == 1

    beam_type = np.full(n,'unknown')
    beam_number = np.full(n,-1)
    if gt == 'gt1r':
        beam_number[s0] = 2
        beam_number[s1] = 1
        beam_type[s0] = 'weak'
        beam_type[s1] = 'strong'

    elif gt == 'gt2r':
        beam_number[s0] = 4
        beam_number[s1] = 3
        beam_type[s0] = 'weak'
        beam_type[s1] = 'strong'

    elif gt == 'gt3r':
        beam_number[s0] = 6
        beam_number[s1] = 5
        beam_type[s0] = 'weak'
        beam_type[s1] = 'strong'

    elif gt == 'gt1l':
        beam_number[s0] = 1
        beam_number[s1] = 2
        beam_type[s0] = 'strong'
        beam_type[s1] = 'weak'

    elif gt == 'gt2l':
        beam_number[s0] = 3
        beam_number[s1] = 4
        beam_type[s0] = 'strong'
        beam_type[s1] = 'weak'

    elif gt == 'gt3l':
        beam_number[s0] = 5
        beam_number[s1] = 6
        beam_type[s0] = 'strong'
        beam_type[s1] = 'weak'

    return beam_number, beam_type


def calc_rdm_segment(t, c, segment_id_beg, segment_id_end, segment_id, ph_index_beg, segment_ph_cnt, debug=0):

    """
    Function to calculate radiometry (rdm)

    Input:
        t - time or delta_time of ATL03, for a given gt num
        c - classification of ATL03 for a given gt num
            ensure that no nans exist
        segment_id_beg - segment_id_beg from ATL08
        segment_id_end - segment_id_end from ATL08
        segment_id - segment_id from ATL03 geolocation/
        ph_index_beg - ph_index_beg from ATL03 geolocation/
        segment_ph_cnt - segment_ph_cnt from ATL03 geolocation/
        debug - val != 0 enables print statements if segments
                do not match from 03 to 08 (see caveats)

    Output:
        n_shots_unique - total number of unique ttg per ATL08 100m bin
        rdm_ground - rdm of ground photons (c==1)
        rdm_veg - rdm of veg photons (c==2)
        rdm_canopy - rdm of canopy photons (c==3)

    Example:
        n_shots_unique, rdm_ground, rdm_veg, rdm_canopy = \
            calc_rdm(t, c, segment_id_beg, segment_id_end, ph_index_beg, segment_ph_cnt, debug=0)

    Caveats:
        Ensure that no nans exist in classification c

        rdm_ground/veg/canopy and n_shots_unique are floating point
        b/c it's possible to have no overlap in 03 and 08 data, in
        which case the radiometry value is NaN; this is implemented by
        initializing rdm vectors are NaN. Thus, they are floating-point-
        valued.

        This functions can handle when 03/08 do not totally overlap,
        or when there is no overlap. That said, one should proceed with
        caution knowing 03 and 08 do not overlap at all. NaN values are
        initialized in rdm vectors based on these cases.

    """

    if np.isnan(c).sum() > 0 and debug:
        print('warning: NaN values found in c')

    rdm_ground = np.full(segment_id_beg.shape, np.nan)
    rdm_veg = np.full(segment_id_beg.shape, np.nan)
    rdm_canopy = np.full(segment_id_beg.shape, np.nan)
    n_shots_unique = np.full(segment_id_beg.shape, np.nan)
    
    n_id = len(segment_id)
    for s in range(len(segment_id_beg)):
#        _, k0 = iu.getClosest(segment_id, [segment_id_beg[s]])
#        _, k1 = iu.getClosest(segment_id, [segment_id_end[s]])
        _, k0 = getClosest(segment_id, [segment_id_beg[s]])
        _, k1 = getClosest(segment_id, [segment_id_end[s]])
        k0, k1 = int(k0), int(k1)
        
        warn = False
        b_edge = False
        if segment_id[k0] < segment_id_beg[s]:
            # left side incomplete
            # cm.pause('beg')
            k = k0
            while segment_id[k] < segment_id_beg[s]:
                k += 1
                if k >= n_id:
                    b_edge = True
                    break

        elif segment_id[k0] > segment_id_beg[s]:
            # print('warning: 03 seg id beg %d > 08 seg id beg %d' % (segment_id[k0], segment_id_beg[s]))
            warn = True

        # else:
            # equal, totally fine

        # if segment_id[k1] != segment_id_end[s]:
        if segment_id[k1] > segment_id_end[s]:
            # right side incomplete
            # cm.pause('end')
            k = k1
            while segment_id[k] > segment_id_end[s]:
                k -= 1
                if k < 0:
                    b_edge = True
                    break

        elif segment_id[k1] < segment_id_end[s]:
            # print('warning: 03 seg id beg %d < 08 seg id beg %d' % (segment_id[k0], segment_id_beg[s]))
            warn = True

        # else:
            # equal, totally fine

        if b_edge and debug:
            # 08 segment is entirely outside of 03 segment data
            print('outside')
            print('03: [%d, %d]' % (segment_id[k0], segment_id[k1]))
            print('08: [%d, %d]' % (segment_id_beg[s], segment_id_end[s]))
            # cm.pause()
            input('enter to continue')
            continue

        if warn and debug:
            print('partial')
            print('03: [%d, %d]' % (segment_id[k0], segment_id[k1]))
            print('08: [%d, %d]' % (segment_id_beg[s], segment_id_end[s]))
            # cm.pause()
            input('enter to continue')

        i0, i1 = ph_index_beg[k0], ph_index_beg[k1] + segment_ph_cnt[k1] - 1

        t_seg = t[i0:i1+1] # inclusive index
        c_seg = c[i0:i1+1]

        n_shots_total_uq = len(np.unique(t_seg))
        n_shots_ground = (c_seg == 1).sum()
        n_shots_veg = (c_seg == 2).sum()
        n_shots_canopy = (c_seg == 3).sum()

        n_shots_unique[s] = n_shots_total_uq
        rdm_ground[s] = float(n_shots_ground / n_shots_total_uq)
        rdm_veg[s] = float(n_shots_veg / n_shots_total_uq)
        rdm_canopy[s] = float(n_shots_canopy / n_shots_total_uq)

    return n_shots_unique, rdm_ground, rdm_veg, rdm_canopy


def unit(v):
    mag = np.linalg.norm(v)
    return np.array([vi/mag for vi in v])

def df_to_np(df, *args):
    rtn = []
    for var in args:
        rtn.append(np.array(df[var]))
    return rtn

def downsample(Fs_ds, t, *data):
    """
    Input:
        Fs_ds - downsampled new frequency
        t - time-series
        *data - one or multiple args to be
                downsampled along with t

    Output:
        t_ds - downsampled time-series
        n_ds - len(t_ds)
        data_out - one or multiple args

    Example:
        dt_ds = 1.0 # sec
        Fs_ds = 1.0 / dt_ds
        t_ds, n_ds, (arg1, arg2, ...) = downsample(Fs_ds, t, arg1, arg2, ...)

    """

    dt_ds = 1.0 / Fs_ds # sec
    tol = 0.01*dt_ds
    t_ds = []
    # data_ds = []
    data_new = []
    num_y = len(data)
    for y in data:
        data_new.append([])

    n = len(t)
    t_prev = t[0]
    for j in range(1,n):
        if abs(t[j] - t_prev) > dt_ds-tol:
            t_ds.append(t[j])
            for i in range(num_y):
                data_new[i].append(data[i][j])
            t_prev = t[j]

    n_ds = len(t_ds)

    # if rtn_numpy:
    if 1:
        t_ds = np.array(t_ds)
        for i in range(num_y):
            data_new[i] = np.array(data_new[i])

    data_out = tuple(data_new)

    if len(data_out) > 1:
        return t_ds, n_ds, data_out
    else:
        return t_ds, n_ds, data_out[0]


def b_filt(b, *args):
    arr_new = []
    for arr in args:
        arr_new.append(arr[b])
    return arr_new
# endDef
    
# Function to interpolate 1d
def interp_vals(input_x, input_y, interp_x, removeThresh=False):
    
    # Remove y values > 1e30
    if(removeThresh):
        indsUnderThresh = input_y <= 1e30
        input_x = input_x[indsUnderThresh]
        input_y = input_y[indsUnderThresh]
    # endIf
    
    # Remove x duplicate values for scipy interp
    indsUnique = np.unique(input_x, return_index=True)[1]
    input_x = input_x[indsUnique]
    input_y = input_y[indsUnique]
    
    # Interpolate delta_time
    f1 = interpolate.interp1d(input_x, input_y, kind='linear', fill_value='extrapolate')
    interp_y = f1(interp_x)
    
    return interp_y

# endDef

# Function to interpolate 2d
def interp_vals2d(input_x, input_y, input_z, interp_x, interp_y, removeThresh=False, removeDuplicates=False):
    
    # Remove y values > 1e30
    if(removeThresh):
        indsUnderThresh = input_y <= removeThresh
        input_x = input_x[indsUnderThresh]
        input_y = input_y[indsUnderThresh]
    # endIf
    
    # Remove x duplicate values for scipy interp
    if(removeDuplicates):
        indsUnique = np.unique(input_x, return_index=True)[1]
        input_x = input_x[indsUnique]
        input_y = input_y[indsUnique]
    # endIf
    
    # Interpolate delta_time
    f1 = interpolate.interp2d(input_x, input_y, input_z, kind='linear', fill_value='extrapolate')
    interp_z = f1(interp_x, interp_y)
    
    return interp_z

# endDef
    

# Function to build output dataframe for csv file
def build_dataframe_for_csv(atl03Data, namelist):
    
    # Gather data together
    datalist_df = pd.DataFrame()
    datalist_df[namelist[0]]  = np.ravel(atl03Data.time)

    if(len(atl03Data.deltaTime)>0):
        datalist_df[namelist[1]]  = np.ravel(atl03Data.deltaTime)
    else:
        datalist_df[namelist[1]]  = np.nan
    # endIf
    
    if(len(atl03Data.segmentID)>0):
        datalist_df[namelist[2]]  = np.ravel(atl03Data.segmentID)
    else:
        datalist_df[namelist[2]]  = np.nan
    # endIf
    
    if(len(atl03Data.gtNum)>0):
        datalist_df[namelist[3]]  = atl03Data.gtNum
    else:
        datalist_df[namelist[3]]  = np.nan
    # endIf
    
    if(len(atl03Data.beamNum)>0):
        datalist_df[namelist[4]]  = atl03Data.beamNum
    else:
        datalist_df[namelist[4]]  = np.nan
    # endIf
    
    if(len(atl03Data.beamStrength)>0):
        datalist_df[namelist[5]]  = atl03Data.beamStrength
    else:
        datalist_df[namelist[5]]  = np.nan
    # endIf
    
    if(len(atl03Data.lat)>0):
        datalist_df[namelist[6]]  = np.ravel(atl03Data.lat)
    else:
        datalist_df[namelist[6]]  = np.nan
    # endIf
    
    if(len(atl03Data.lon)>0):
        datalist_df[namelist[7]]  = np.ravel(atl03Data.lon)
    else:
        datalist_df[namelist[7]]  = np.nan
    # endIf
    
    if(len(atl03Data.easting)>0):
        datalist_df[namelist[8]]  = np.ravel(atl03Data.easting)
    else:
        datalist_df[namelist[8]]  = np.nan
    # endIf
    
    if(len(atl03Data.northing)>0):
        datalist_df[namelist[9]]  = np.ravel(atl03Data.northing)
    else:
        datalist_df[namelist[9]]  = np.nan
    # endIf
    
    if(len(atl03Data.zone)>0):
        datalist_df[namelist[10]] = atl03Data.zone
    else:
        datalist_df[namelist[10]]  = np.nan
    # endIf
    
    if(len(atl03Data.hemi)>0):
        datalist_df[namelist[11]] = atl03Data.hemi
    else:
        datalist_df[namelist[11]]  = np.nan
    # endIf
    
    if(len(atl03Data.crossTrack)>0):
        datalist_df[namelist[12]] = np.ravel(atl03Data.crossTrack)
    else:
        datalist_df[namelist[12]]  = np.nan
    # endIf
    
    if(len(atl03Data.alongTrack)>0):
        datalist_df[namelist[13]] = np.ravel(atl03Data.alongTrack)
    else:
        datalist_df[namelist[13]]  = np.nan
    # endIf
    
    if(len(atl03Data.z)>0):
        datalist_df[namelist[14]] = np.ravel(atl03Data.z)
    else:
        datalist_df[namelist[14]]  = np.nan
    # endIf
    
    if(len(atl03Data.zMsl)>0):
        datalist_df[namelist[15]] = np.ravel(atl03Data.zMsl)
    else:
        datalist_df[namelist[15]] = np.nan
    # endIf
    
    if(len(atl03Data.classification)>0):
        datalist_df[namelist[16]] = np.ravel(atl03Data.classification)
    else:
        datalist_df[namelist[16]]  = np.nan
    # endIf
    
    if(len(atl03Data.signalConf)>0):
        datalist_df[namelist[17]] = np.ravel(atl03Data.signalConf)
    else:
        datalist_df[namelist[17]]  = np.nan
    # endIf
    
    if(len(atl03Data.solar_elev)>0):
        datalist_df[namelist[18]] = np.ravel(atl03Data.solar_elev)
    else:
        datalist_df[namelist[18]]  = np.nan
    # endIf
    
    return datalist_df

# endDef

if __name__ == "__main__":
    print("Test")
# end