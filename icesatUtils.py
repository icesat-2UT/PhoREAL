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
import os
import pyproj as proj
import warnings
from scipy import interpolate
import pandas as pd
import h5py
import ctypes
from numpy.ctypeslib import ndpointer 
import copy
from osgeo import ogr
from gui_addins import (superFilterFile_windows, superFilterFile_linux)

EPSG_ARCTIC = '3413'
EPSG_ANTARCTIC = '3976'
EPSG_ERR = '-1'
DEG_ARCTIC = 84.0
DEG_ANTARCTIC = -80.0

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


##### Function to represent 2 numbers as 1 unique number
def cantorPairing(arrayIn):
    
    # Do Cantor pairing algorithm to express two values as one unique value
    vectorOut = 0.5 * (arrayIn[:,0] + arrayIn[:,1]) * (arrayIn[:,0] + arrayIn[:,1] + 1) + arrayIn[:,1]
    
    return vectorOut


##### Function to determine members of one array in another
def ismember(a_vec, b_vec, methodType = 'normal'):
    
    """ MATLAB equivalent ismember function """
    
    # Combine multi column arrays into a 1-D array of strings if necessary
    # This will ensure unique rows when using np.isin below
    if(methodType.lower() == 'rows'):
        
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
            # endIf
        # endFor
    # endIf
    
    # Find which values in a_vec are present in b_vec
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
    allph_classed = (np.zeros(sizeOutput + 1).astype(int)) - 1
    
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
def find_utm_zone_arr(lon, lat, mode=True, module=None, m=None):

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

        arctic = len(lat[lat >= 84])
        nhem = len(lat[(lat < 84) & (lat >= 0)])
        shem = len(lat[(lat > -80) & (lat < 0)])
        antarctic = len(lat[lat < -80])
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

def identifyEPSG(hemi,zone):
    if hemi == 'N':
        outstring = 'epsg:326'
    elif hemi == "S":
        outstring = 'epsg:327'
    else:
        outstring = ''
        print('error')
    
    outstring = outstring + (str(zone))
    
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
        xx, yy = transform(epsg_in, epsg_out, lon, lat)
        
    else:
        
        # Get UTM coords
        xx, yy, epsg_out = wgs84_to_utm_find_and_transform(lon, lat)
        
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
            # endIf
            
        # endIf
        
    # endIf
    
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
    x = atlTruthData.easting
    y = atlTruthData.northing
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
def getRaster(x, y, z, resolution, method, fillValue = -999, time = [], xAllArray = [], yAllArray = []):
    
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
    #   - range
    #   - std (standard deviation)
    #   - numel (number of elements)
    # fillValue = value to fill in empty grid cells
    # time = secondary array (like z) to perform operation on in each grid cell
    # xAllArray = force output X grid cells (use np.arange(start, stop + 1, step))
    # yAllArray = force output Y grid cells (use np.arange(start, stop + 1, step))
    
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
        npOperation = np.min
    elif(method.lower() == 'max'):
        npOperation = np.max
    elif(method.lower() == 'mean'):
        npOperation = np.mean
    elif(method.lower() == 'median'):
        npOperation = np.median
    elif(method.lower() == 'range'):
        npOperation = np.range
    elif(method.lower() == 'std'):
        npOperation = np.std
    elif(method.lower() == 'numel'):
        npOperation = np.size
    else:
        npOperation = np.mean
    # EndIf
            
    # Round all incoming X,Y data
    xRnd = (np.round(x/xResolution)*xResolution).astype(int)
    yRnd = (np.round(y/yResolution)*yResolution).astype(int)
    
    # Get output X,Y grid cells
    if(any(xAllArray) and any(yAllArray)):
        
        xAll = xAllArray
        yAll = yAllArray
        
    else:
        
        # Get min,max of rounded X,Y data
        xRndMin = xRnd.min()
        xRndMax = xRnd.max()
        yRndMin = yRnd.min()
        yRndMax = yRnd.max()
        
        # Get all possible grid combinations
        xAll = np.arange(xRndMin, xRndMax + xResolution, xResolution)
        yAll = np.arange(yRndMax, yRndMin - yResolution, -yResolution)
        
    # endIf
        
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
    
    uniqueCombos = uniqueCombos.astype('int')
    xyAll = xyAll.astype('int')
    
    # Find index locations of unique incoming X,Y data in X,Y array of all pts
    _, indsToFill = ismember(uniqueCombos,xyAll,'rows')
    indsToFill = np.c_[indsToFill]
    
    # Grid Z data, populate into raster array, and reshape 
    df = pd.DataFrame(np.column_stack([z, uniqueGroups]), columns=['z', 'unique_groups'])
    zGroups = df.groupby('unique_groups')
    zOut = zGroups.agg(npOperation)
    zOut = np.array(zOut)
    zSplit = np.c_[zOut[:,0]]
    zRaster = fillValue*np.ones((numRows*numCols,1))
    zRaster[indsToFill,0] = zSplit
    rasterDataZ = np.reshape(zRaster,(numRows,numCols))
    
    # Grid 'time' array if necessary
    if(any(time)):
        
        df = pd.DataFrame(np.column_stack([time, uniqueGroups]), columns=['time', 'unique_groups'])
        tGroups = df.groupby('unique_groups')
        tOut = tGroups.agg(npOperation)
        tOut = np.array(tOut)
        tSplit = np.c_[tOut[:,0]] 
        tRaster = fillValue*np.ones((numRows*numCols,1))
        tRaster[indsToFill,0] = tSplit
        rasterDataT = np.reshape(tRaster,(numRows,numCols))
        
    else:
        
        rasterDataT = []
        
    # EndIf

    # Return output
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

def sortAtlMeasured(atlMeasuredData):
    
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
   
    
    print("Complete")
    return atlMeasuredData
    
def sortAtlTruth(atlTruthData):
    print("Sorting ATL Truth Data...", end = " ")   
    sort_index = np.argsort(atlTruthData.alongTrack, axis = 0)
    atlTruthData.alongTrack = atlTruthData.alongTrack[sort_index[:,0]]
    atlTruthData.crossTrack = atlTruthData.crossTrack[sort_index[:,0]]
    atlTruthData.northing = atlTruthData.northing[sort_index[:,0]]
    atlTruthData.easting = atlTruthData.easting[sort_index[:,0]]
    atlTruthData.z = atlTruthData.z[sort_index[:,0]]
    atlTruthData.intensity = atlTruthData.intensity[sort_index[:,0]]
    atlTruthData.classification = atlTruthData.classification[sort_index[:,0]]

    print("Complete")
    return atlTruthData

def indexMatch(measuredArray,truthArray):
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
    print("Complete")
    return np.array(C).astype(int)

    
def superFilter(atlMeasuredData_in, atlTruthData_in, xBuf = 7, classCode = []):
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
        atlTruthData.z = atlTruthData.z[filter_data]
        atlTruthData.intensity = atlTruthData.intensity[filter_data]
        atlTruthData.classification = atlTruthData.classification[filter_data]
    else:
        indexMatches = indexMatch(atlMeasuredData.alongTrack \
                                    ,atlTruthData.alongTrack)
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
        atlTruthData.z = atlTruthData.z[filter_data]
        atlTruthData.intensity = atlTruthData.intensity[filter_data]
        atlTruthData.classification = atlTruthData.classification[filter_data]  
    print("Superfilter complete")
    #Return truthdata
    return atlTruthData, atlMeasuredData

def getBins(atlMeasuredData_in,atlTruthData_in, binsize, 
            measClassification = None, truthClassification = None, 
            operation = 'median', matchfilter = True, nanfilter = True):
  
    if measClassification:
        classlist = np.array(measClassification)
        classfilter = np.isin(atlMeasuredData_in.classification, classlist)
        measuredAlongTrack = atlMeasuredData_in.crossTrack[classfilter]
        measuredCrossTrack = atlMeasuredData_in.alongTrack[classfilter]
        measuredVal = atlMeasuredData_in.z[classfilter]
 
    else:
        measuredAlongTrack = atlMeasuredData_in.crossTrack
        measuredCrossTrack = atlMeasuredData_in.alongTrack
        measuredVal = atlMeasuredData_in.z 
    
    if truthClassification:
        classlistTruth = np.array(truthClassification)
        classfilterTruth = np.isin(atlTruthData_in.classification, classlistTruth)
        truthAlongTrack = atlTruthData_in.crossTrack[classfilterTruth]
        truthCrossTrack = atlTruthData_in.alongTrack[classfilterTruth]
        truthVal = atlTruthData_in.z[classfilterTruth]
    else:
        truthAlongTrack = atlTruthData_in.crossTrack
        truthCrossTrack = atlTruthData_in.alongTrack
        truthVal = atlTruthData_in.z
    
    #Rasterize Measured and Truth
    try:
        binMeasured = getRaster(measuredAlongTrack[:,0], measuredCrossTrack[:,0], 
                                measuredVal[:,0], [100000,binsize], str(operation), 
                                fillValue = -999, time = [])
    except:
        binMeasured = getRaster(measuredAlongTrack, measuredCrossTrack, 
                            measuredVal, [100000,binsize], str(operation), 
                            fillValue = -999, time = [])
    try:     
        binTruth = getRaster(truthAlongTrack[:,0], truthCrossTrack[:,0], 
                             truthVal[:,0], [100000,binsize], str(operation), 
                             fillValue = -999, time = [])
    except:
        binTruth = getRaster(truthAlongTrack, truthCrossTrack, 
                     truthVal, [100000,binsize], str(operation), 
                     fillValue = -999, time = [])
    if matchfilter == True:
        idx_atl03 = np.isin(binMeasured.y, binTruth.y)
        idx_truth = np.isin(binTruth.y, binMeasured.y)
        
        binMeasured.x = binMeasured.x[idx_atl03]
        binMeasured.y = binMeasured.y[idx_atl03]
        binMeasured.grid = binMeasured.grid[idx_atl03]
        
        binTruth.x = binTruth.x[idx_truth]
        binTruth.y = binTruth.y[idx_truth]
        binTruth.grid = binTruth.grid[idx_truth]
        
    if nanfilter == True:
        nan_filter = np.where((binMeasured.grid > -998) & (binTruth.grid > -998) & 
                              (binMeasured.grid < 9000) & (binTruth.grid < 9000))
        
        binMeasured.x = binMeasured.x[nan_filter]    
        binMeasured.y = binMeasured.y[nan_filter]
        binMeasured.grid = binMeasured.grid[nan_filter]
        
        binTruth.x = binTruth.x[nan_filter]
        binTruth.y = binTruth.y[nan_filter]
        binTruth.grid = binTruth.grid[nan_filter]
    return binTruth, binMeasured

def createScalogram(atlMeasured, altTruth, measClassification = [1],
                    truthClassification = [2], scalelist = []):
    superTruth, sortedMeasured = superFilter(atlMeasured, altTruth, xBuf = 7, 
                                             classCode = measClassification)

    if len(scalelist) == 0:
        scalelist = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
    minbin = 5
    count = 0
    for scale in scalelist:
        print(str(scale))
        binTruth, binMeasured = getBins(sortedMeasured,superTruth, scale, 
                measClassification, truthClassification = [2], 
                operation = 'median', matchfilter = True, nanfilter = False)
        

        binTruth.grid[binTruth.grid <= -999] = np.nan
        binMeasured.grid[binMeasured.grid <= -999] = np.nan
        me = binTruth.grid - binMeasured.grid
        
        
        repeat = scale/minbin
        
        me = np.repeat(me,[repeat],axis=0)
        if count == 0:
            lencap = np.int(len(me))
            scalogram = me
        else:
            me = me[0:lencap]
            print(me.shape)
            print(scalogram.shape)
            scalogram = np.vstack((scalogram,me))

        count = count + 1

    return scalogram


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

# endDef
    
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


    
if __name__ == "__main__":
    print("Test")
