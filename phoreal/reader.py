#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that contains ATL03 and ATL 08 H5 Reader functions for PhoREAL

Copyright 2019 Applied Research Laboratories, University of Texas at Austin

This package is free software; the copyright holder gives unlimited
permission to copy and/or distribute, with or without modification, as
long as this notice is preserved.

Authors:
    Eric Guenther
    Mike Alonzo
    
Date: February 27, 2019
"""

import pandas as pd
import numpy as np
import os
import h5py

from phoreal.utils import getH5Keys, ismember
from phoreal.io import readAtl03H5, readAtlH5, readAtl03DataMapping, readAtl08DataMapping
from phoreal.utils import getAtl08Mapping, wgs84_to_utm_find_and_transform,\
    wgs84_to_epsg_transform, getCoordRotFwd, getNameParts, get_h5_meta,\
        identify_hemi_zone, identify_hemi_zone
from phoreal.io import GtToBeamNum, GtToBeamSW, readTruthRegionsTxtFile, writeLas
# from icesatIO import readHeaderMatFile                           
# from getAtlMeasuredSwath_auto import atl03Struct as Atl03StructLegacy

class AtlRotationStruct:
    
    # Define class with designated fields
    def __init__(self, R_mat, xRotPt, yRotPt, desiredAngle, phi):
        
        self.R_mat = R_mat
        self.xRotPt = xRotPt
        self.yRotPt = yRotPt
        self.desiredAngle = desiredAngle
        self.phi = phi


class AtlStruct:
        
    # Define class with designated fields
    def __init__(self, df, gtNum, beamNum, beamStrength, epsg, zone, hemi, 
                 atlFilePath, atlFileName, 
                 trackDirection, atlProduct, alth5Info, dataIsMapped, 
                 rotation_data, ancillary=None, orbit_info=None):
            
        self.df = df
        self.gtNum = gtNum
        self.beamNum = beamNum
        self.beamStrength = beamStrength
        self.epsg = epsg
        self.zone = zone
        self.hemi = hemi
        self.atlFilePath = atlFilePath
        self.atlFileName = atlFileName
        self.trackDirection = trackDirection
        self.atlProduct = atlProduct
        self.atlVersion = alth5Info.atlVersion
        self.year = alth5Info.year
        self.month = alth5Info.month
        self.day = alth5Info.day
        self.hour = alth5Info.hour
        self.minute = alth5Info.minute
        self.second = alth5Info.second
        self.trackNum = alth5Info.trackNum
        self.unknown = alth5Info.unknown
        self.releaseNum = alth5Info.releaseNum
        self.incrementNum = alth5Info.incrementNum
        self.dataIsMapped = dataIsMapped
        self.rotationData = rotation_data
        self.ancillary = ancillary
        self.orbit_info = orbit_info
        
    def trim_by_lat(self, min_lat, max_lat):
        if self.atlProduct == 'ATL03':
            min_lat = np.min([min_lat, max_lat])
            max_lat = np.max([min_lat, max_lat])
            self.df = self.df = self.df[self.df.lat_ph > min_lat]
            self.df = self.df = self.df[self.df.lat_ph < max_lat]
            self.df = self.df.reset_index()
            self.df, self.rotationData = get_atl_alongtrack(self.df, self)
        elif self.atlProduct == 'ATL08':
            min_lat = np.min([min_lat, max_lat])
            max_lat = np.max([min_lat, max_lat])
            self.df = self.df = self.df[self.df.latitude > min_lat]
            self.df = self.df = self.df[self.df.latitude < max_lat]
            
    def trim_by_lon(self, min_lon, max_lon):
        if self.atlProduct == 'ATL03':
            min_lon = np.min([min_lon, max_lon])
            max_lon = np.max([min_lon, max_lon])
            self.df = self.df = self.df[self.df.lon_ph > min_lon]
            self.df = self.df = self.df[self.df.lon_ph < max_lon]
            self.df = self.df.reset_index()
            self.df, self.rotationData = get_atl_alongtrack(self.df, self)
        elif self.atlProduct == 'ATL08':
            min_lon = np.min([min_lon, max_lon])
            max_lon = np.max([min_lon, max_lon])
            self.df = self.df = self.df[self.df.longitude > min_lon]
            self.df = self.df = self.df[self.df.longitude < max_lon]   
            
    def to_csv(self, output_csv):
        self.df.to_csv(output_csv)
        
    def to_mat(self, output_mat):
        convert_df_to_mat(self.df,output_mat)
    
    def quick_plot(self):
        pass
    
# Read ATL03 Heights, put in Pandas DF
def read_atl03_heights_data(atl03filepath, gt):
    # Iterate through keys for "Heights"
    keys = getH5Keys(atl03filepath,gt + '/heights')
    
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtl03H5(atl03filepath, '/heights/' + key, gt)
        if idx == 0:
            df = pd.DataFrame(data,columns=[key])
        else:
            df = pd.concat([df,pd.DataFrame(data,columns=[key])],axis=1)
    return df

def read_atl03_geolocation(atl03filepath, gt):
    # Iterate through keys for "Heights"
    keys = getH5Keys(atl03filepath,gt + '/geolocation')
    key_info = get_H5_keys_info(atl03filepath, gt + '/geolocation')
    
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtl03H5(atl03filepath, '/geolocation/' + key, gt)
        if key_info[idx][1] != 'Group':
            if idx == 0:
                df = pd.DataFrame(data,columns=[key.split('/')[-1]])
            else:
                if len(data.shape) == 2:
                    cols = data.shape[1]
                    for idx2 in range(0,cols):
                        df = pd.concat([df,pd.DataFrame(data[:,idx2],columns=\
                                                        [key.split('/')[-1] +\
                                                         '_' + str(idx2)])],
                                       axis=1)                        
                else:
                    df = pd.concat([df,pd.DataFrame(data,columns=\
                                                    [key.split('/')[-1]])],
                                   axis=1)
    return df

    
# Read ATL03 Heights, put in Pandas DF
def read_atl08_land_segments(atl08filepath, gt):
    # Iterate through keys for "Land Segments"
    keys = getH5Keys(atl08filepath,gt + '/land_segments')
    key_info = get_H5_keys_info(atl08filepath,gt + '/land_segments')
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtl03H5(atl08filepath, '/land_segments/' + key, gt)
        if key_info[idx][1] != 'Group':
            if idx == 0:
                df = pd.DataFrame(data,columns=[key.split('/')[-1]])
            else:
                if len(data.shape) == 2:
                    cols = data.shape[1]
                    for idx2 in range(0,cols):
                        df = pd.concat([df,pd.DataFrame(data[:,idx2],columns=\
                                                        [key.split('/')[-1] +\
                                                         '_' + str(idx2)])],
                                       axis=1)                        
                else:
                    df = pd.concat([df,pd.DataFrame(data,columns=\
                                                    [key.split('/')[-1]])],
                                   axis=1)
    return df

# Read ATL03 Heights, put in Pandas DF
def read_atl09_hr_profile(atl09filepath, gt):
    # Iterate through keys for "Land Segments"
    subgroup = 'profile_' + gt[2] + '/high_rate/'
    keys = getH5Keys(atl09filepath, subgroup)
    key_info = get_H5_keys_info(atl09filepath, subgroup)
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtlH5(atl09filepath, subgroup + '/' + key + '/')
        if key == 'ds_layers':
            ds_layers = np.array(data)
        elif key == 'ds_va_bin_h':
            ds_va_bin_h = np.array(data)
        elif key == 'cab_prof':
            cab_prof = np.array(data)
        elif key == 'density_pass1':
            density_pass1 = np.array(data)
        elif key == 'density_pass2':
            density_pass2 = np.array(data)
        else:
            
            if key_info[idx][1] != 'Group':
                if idx == 0:
                    df = pd.DataFrame(data,columns=[key.split('/')[-1]])
                else:
                    if len(data.shape) == 2:
                        cols = data.shape[1]
                        for idx2 in range(0,cols):
                            df = pd.concat(
                                [df,pd.DataFrame(data[:,idx2],columns=\
                                                 [key.split('/')[-1] +'_' + 
                                                  str(idx2)])],axis=1)                        
                    else:
                        df = pd.concat(
                            [df,pd.DataFrame(data,columns=\
                                             [key.split('/')[-1]])],axis=1)
        
    return df, ds_layers, ds_va_bin_h, cab_prof, density_pass1, density_pass2

def read_atl09_ancillary_data(atl09filepath):
    # Iterate through keys for ancillary data
    subgroup = '/ancillary_data/'
    keys = getH5Keys(atl09filepath, subgroup)
    key_info = get_H5_keys_info(atl09filepath, subgroup)
    
    byte_encoded = ['control', 'data_start_utc', 'data_end_utc',
                    'granule_start_utc', 'granule_end_utc']
    
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtlH5(atl09filepath, subgroup + '/' + key + '/')
        
        if np.isin(key, ['release', 'version']):
            data = int(data)
            
        if np.isin(key, byte_encoded):
            data = data[0]             
            data = data.decode('utf-8')
        
        if key_info[idx][1] != 'Group':
            if len(key.split('/')) == 1:
                if idx == 0:
                    df = pd.Series(data,
                                   index=[key.split('/')[-1]],
                                   dtype=object)
                else:
                    df = pd.concat(
                        [df, pd.Series(data, 
                                       index=[key.split('/')[-1]],
                                       dtype=object)])
            
    return df

def read_atl09_orbit_info(atl09filepath):

    subgroup = '/orbit_info/'
    keys = getH5Keys(atl09filepath, subgroup)
    # key_info = get_H5_keys_info(atl09filepath, subgroup)

    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtlH5(atl09filepath, subgroup + '/' + key + '/')
        if idx == 0:
            df = pd.Series(data, index=[key], dtype=object)
        else:
            df_key = pd.Series(data, index=[key], dtype=object)
            df = pd.concat([df, df_key])
    
    return df

# Map classifications from ATL08, map back to ATL03 Photons
def get_atl03_classification(atl03filepath, atl08filepath, gt):
    # Read ATL03 metrics for class mapping
    f = h5py.File(atl03filepath, 'r')
    atl03_ph_index_beg  = np.array(f[gt + '/geolocation/ph_index_beg'])
    atl03_segment_id = np.array(f[gt + '/geolocation/segment_id'])
    atl03_heights = np.array(f[gt + '/heights/h_ph'])
    
    # Read ATL08 for class mapping
    f = h5py.File(atl08filepath, 'r')
    atl08_classed_pc_indx = np.array(f[gt + '/signal_photons/classed_pc_indx'])
    atl08_classed_pc_flag = np.array(f[gt + '/signal_photons/classed_pc_flag'])
    atl08_segment_id = np.array(f[gt + '/signal_photons/ph_segment_id'])
    
    # Map ATL08 classifications to ATL03 Photons
    allph_classed = getAtl08Mapping(atl03_ph_index_beg, atl03_segment_id, 
                                    atl08_classed_pc_indx, 
                                    atl08_classed_pc_flag, 
                                    atl08_segment_id)
    
    if len(allph_classed) < len(atl03_heights):
        n_zeros = len(atl03_heights) - len(allph_classed)
        zeros = np.zeros(n_zeros)
        allph_classed = np.append(allph_classed, zeros)
    
    return allph_classed

def merge_label_to_df(atl03filepath, atl08filepath, gt, df):
    allph_classed = get_atl03_classification(atl03filepath, atl08filepath, gt)
    # Add classifications to ATL03 DF
    df = pd.concat([df,pd.DataFrame(allph_classed,
                                    columns=['classification'])],axis=1)
    
    # Replace nan with -1 (unclassified)
    df.replace({'classification' : np.nan}, -1)
    return df

def get_atl03_heights_offset(atl03filepath, atl08filepath, gt):
    # Read ATL03 metrics for class mapping
    f = h5py.File(atl03filepath, 'r')
    atl03_ph_index_beg  = np.array(f[gt + '/geolocation/ph_index_beg'])
    atl03_segment_id = np.array(f[gt + '/geolocation/segment_id'])
    atl03_heights = np.array(f[gt + '/heights/h_ph'])

    # Read ATL08 for class mapping
    f = h5py.File(atl08filepath, 'r')
    atl08_classed_pc_indx = np.array(f[gt + '/signal_photons/classed_pc_indx'])
    atl08_heights = np.array(f[gt + '/signal_photons/ph_h'])
    atl08_segment_id = np.array(f[gt + '/signal_photons/ph_segment_id'])
            
    # Get ATL03 data
    indsNotZero = atl03_ph_index_beg != 0
    atl03_ph_index_beg = atl03_ph_index_beg[indsNotZero];
    atl03_segment_id = atl03_segment_id[indsNotZero];
    
    # Find ATL08 segments that have ATL03 segments
    atl03SegsIn08TF, atl03SegsIn08Inds = ismember(atl08_segment_id,atl03_segment_id)
    
    # Get ATL08 classed indices and values
    atl08classed_inds = atl08_classed_pc_indx[atl03SegsIn08TF]
    atl08classed_vals = atl08_heights[atl03SegsIn08TF]

    # Determine new mapping into ATL03 data
    atl03_ph_beg_inds = atl03SegsIn08Inds;
    atl03_ph_beg_val = atl03_ph_index_beg[atl03_ph_beg_inds]
    newMapping = atl08classed_inds + atl03_ph_beg_val - 2
    
    # Get max size of output array
    sizeOutput = newMapping[-1]
    
    # Pre-populate all photon classed array with zeroes
    allph_heights = (np.zeros(sizeOutput + 1)) 
    
    # Populate all photon classed array from ATL08 classifications
    allph_heights[newMapping] = atl08classed_vals
    allph_heights[allph_heights == 3.4028234663852886e+38] = np.nan

    if len(allph_heights) < len(atl03_heights):
        n_zeros = len(atl03_heights) - len(allph_heights)
        zeros = np.zeros(n_zeros)
        zeros = zeros * np.nan
        allph_heights = np.append(allph_heights, zeros)
    
    return allph_heights

def get_atl03_rate(atl03filepath, gt):
    f = h5py.File(atl03filepath, 'r')   
    bihr = np.asarray(f[gt + '/bckgrd_atlas/bckgrd_int_height_reduced'])
    bcr = np.asarray(f[gt + '/bckgrd_atlas/bckgrd_counts_reduced'])
    bapmc = np.asarray(f[gt + '/bckgrd_atlas/pce_mframe_cnt'])
    hpmc = np.asarray(f[gt + '/heights/pce_mframe_cnt'])
   
    # Calculate rate and assign to the bckgrd_atlas rate
    rate = bcr / bihr

    # Assign bckgrd_atlas attributes to photon level    
    tf, inds = ismember(hpmc, bapmc)
    
    ph_bihr = bihr[inds]
    ph_bcr = bcr[inds]
    ph_rate = rate[inds]

    ph_bihr = np.concatenate([np.zeros(len(tf[tf == False])),ph_bihr])
    ph_bcr = np.concatenate([np.zeros(len(tf[tf == False])),ph_bcr])
    ph_rate = np.concatenate([np.zeros(len(tf[tf == False])),ph_rate])
    
    return ph_bihr, ph_bcr, ph_rate

def merge_norm_h_to_df(atl03filepath, atl08filepath, gt, df):
    # Get allph_heights
    allph_heights = get_atl03_heights_offset(atl03filepath, atl08filepath, gt)
    # Add classifications to ATL03 DF
    df = pd.concat([df,pd.DataFrame(allph_heights,
                                    columns=['norm_h'])],axis=1)
    return df

def get_atl03_segment_id(atl03filepath, gt):
    f = h5py.File(atl03filepath, 'r')
    h_ph = np.asarray(f[gt + '/heights/h_ph'])
    segment_ph_count = np.array(f[gt + '/geolocation/segment_ph_cnt'])
    atl03_segment_id = np.array(f[gt + '/geolocation/segment_id'])
    atl03_ph_index_beg = np.array(f[gt + '/geolocation/ph_index_beg'])    
 
    # Get segment ID to photon level
    h_seg = np.zeros(len(h_ph))
    for i in range(0,len(atl03_segment_id)):
        if atl03_ph_index_beg[i] > 0:
            h_seg[atl03_ph_index_beg[i]-1:atl03_ph_index_beg[i]-1 +\
                      segment_ph_count[i]] = atl03_segment_id[i]
    h_seg = np.int32(h_seg) 
    return h_seg

def get_atl03_dist_ph_along(atl03filepath, gt):
    f = h5py.File(atl03filepath, 'r')
    h_ph = np.asarray(f[gt + '/heights/h_ph'])
    segment_ph_count = np.array(f[gt + '/geolocation/segment_ph_cnt'])
    atl03_segment_id = np.array(f[gt + '/geolocation/segment_dist_x'])
    atl03_ph_index_beg = np.array(f[gt + '/geolocation/ph_index_beg'])    
 
    # Get segment ID to photon level
    h_seg = np.zeros(len(h_ph))
    for i in range(0,len(atl03_segment_id)):
        if atl03_ph_index_beg[i] > 0:
            h_seg[atl03_ph_index_beg[i]-1:atl03_ph_index_beg[i]-1 +\
                      segment_ph_count[i]] = atl03_segment_id[i]
    # h_seg = np.int32(h_seg) 
    return h_seg

def merge_seg_id_to_df(atl03filepath, gt, df):
    # Get Segment ID per photon
    h_seg = get_atl03_segment_id(atl03filepath, gt)
    
    # Add classifications to ATL03 DF
    df = pd.concat([df,pd.DataFrame(h_seg,
                                    columns=['seg_id'])],axis=1)
    return df

# Calculate alongtrack time
def get_atl_time(df):
    delta_time = np.array(df['delta_time'])
    min_detla_time = np.min(delta_time[np.nonzero(delta_time)])
    time = delta_time - min_detla_time
    
    df = pd.concat([df,pd.DataFrame(time,
                                    columns=['time'])],axis=1)
    return df
    
# Calcualte Easting/Northing
def get_atl_coords(df, epsg = None):
    columns = list(df.columns)
    
    if 'lon_ph' in columns:
        lon = np.array(df['lon_ph'])
        lat = np.array(df['lat_ph'])
    elif 'longitude' in columns:
        lon = np.array(df['longitude'])
        lat = np.array(df['latitude'])
    elif 'reference_photon_lon' in columns:
        lon = np.array(df['reference_photon_lon'])
        lat = np.array(df['reference_photon_lat'])        
    
    # Specify EPSG Code or automatically find zone
    if epsg:
        xcoord, ycoord = wgs84_to_epsg_transform(epsg, lon, lat)
    else:
        xcoord, ycoord, epsg = wgs84_to_utm_find_and_transform(lon, lat)
        
    if 'easting' not in columns:
        df = pd.concat([df,pd.DataFrame(xcoord,
                                        columns=['easting'])],axis=1)
        df = pd.concat([df,pd.DataFrame(ycoord,
                                        columns=['northing'])],axis=1)
    else:
        print('Warning: Overwritting Existing Coordinates')
        df = df.drop(columns = ['easting'])
        df = df.drop(columns = ['northing'])

        df = pd.concat([df,pd.DataFrame(xcoord,
                                        columns=['easting'])],axis=1)
        df = pd.concat([df,pd.DataFrame(ycoord,
                                        columns=['northing'])],axis=1)        
    
    return df, epsg

def get_atl_alongtrack(df, atl03struct = None):
    easting = np.array(df['easting'])
    northing = np.array(df['northing'])
    
    if atl03struct:
        R_mat = atl03struct.rotationData.R_mat
        xRotPt = atl03struct.rotationData.xRotPt
        yRotPt = atl03struct.rotationData.yRotPt
        desiredAngle = 90
        crossTrack, alongTrack, R_mat, xRotPt, yRotPt, phi = \
            getCoordRotFwd(easting, northing, R_mat, xRotPt, yRotPt, [])    
    else:
        desiredAngle = 90
        crossTrack, alongTrack, R_mat, xRotPt, yRotPt, phi = \
            getCoordRotFwd(easting, northing, [], [], [], desiredAngle)

    if 'crosstrack' not in list(df.columns):
        df = pd.concat([df,pd.DataFrame(crossTrack,
                                    columns=['crosstrack'])],axis=1)
        df = pd.concat([df,pd.DataFrame(alongTrack,
                                        columns=['alongtrack'])],axis=1)
    else:
        print('Warning: Overwriting Existing Alongtrack/Crosstrack')
        df = df.drop(columns = ['crosstrack'])
        df = df.drop(columns = ['alongtrack'])
        df = pd.concat([df,pd.DataFrame(crossTrack,
                                        columns=['crosstrack'])],axis=1)
        df = pd.concat([df,pd.DataFrame(alongTrack,
                                        columns=['alongtrack'])],axis=1)  

    
    rotation_data = AtlRotationStruct(R_mat, xRotPt, yRotPt, desiredAngle, phi)
    
    return df, rotation_data
    
# def get_atl03_df(atl03filepath, atl08filepath, gt, epsg = None):
#     df = read_atl03_heights_data(atl03filepath, gt)
#     df = merge_label_to_df(atl03filepath, atl08filepath, gt, df)
#     df = merge_norm_h_to_df(atl03filepath, atl08filepath, gt, df)
#     df = merge_seg_id_to_df(atl03filepath, gt, df)

#     df = get_atl_time(df)
#     df, epsg = get_atl_coords(df, epsg = None)
#     df, rotationData = get_atl_alongtrack(df)
#     return df

def get_direction(lat):
    if(np.abs(lat[-1]) > np.abs(lat[0])):
        track_direction = 'Ascending'
    else:
        track_direction = 'Descending'
    return track_direction

def get_file_name(filepath):
    filepath = os.path.normpath(os.path.abspath(filepath))
    filename = os.path.splitext(os.path.basename(filepath))[0]
    return filename
  
def get_kml_region(lat,lon, kml_bounds_txt):
    # Determine if ATL03 track goes over a lidar truth region
    kmlBoundsTextFile = kml_bounds_txt
    kmlRegionName = False
    headerFilePath = False
    truthFilePath = False
    
    try:
        if kmlBoundsTextFile and (os.path.exists(kmlBoundsTextFile)):
            
            # Message to user
            print('   Finding Truth Region...')
            
            # Read kmlBounds.txt file and get contents
            kmlInfo = readTruthRegionsTxtFile(kmlBoundsTextFile)
            
            # Loop through kmlBounds.txt and find matching TRUTH area
            maxCounter = len(kmlInfo.regionName) - 1
            counter = 0
            while(not kmlRegionName):
                latInFile = (lat >= kmlInfo.latMin[counter]) & \
                    (lat <= kmlInfo.latMax[counter])
                lonInFile = (lon >= kmlInfo.lonMin[counter]) & \
                    (lon <= kmlInfo.lonMax[counter])
                trackInRegion = any(latInFile & lonInFile)
                if(trackInRegion):
                    
                    # Get truth region info
                    kmlRegionName = kmlInfo.regionName[counter]
    #                headerFilePath = \
    #                    os.path.normpath(kmlInfo.headerFilePath[counter])
    #                truthFilePath = \
    #                    os.path.normpath(kmlInfo.truthFilePath[counter])
                    
                    # Print truth region
                    print('   Truth File Region: %s' % kmlRegionName)
                
                if(counter >= maxCounter):
                    print('   No Truth File Region Found in kmlBounds.txt')
                    break
                counter += 1           
        else:
            kmlRegionName = None
            headerFilePath = None
            truthFilePath = None
    except:
            kmlRegionName = None
            headerFilePath = None
            truthFilePath = None        
        
            # Could not read kmlBounds.txt file
        
    return kmlRegionName, headerFilePath, truthFilePath
    
def write_atl03_las(atlstruct, outpath):
    xx = np.array(atlstruct.df.easting)
    yy = np.array(atlstruct.df.northing)
    zz = np.array(atlstruct.df.h_ph)
    cc = np.array(atlstruct.df.classification)
    ii = np.array(atlstruct.df.signal_conf_ph)
    sigconf = np.array(atlstruct.df.signal_conf_ph)
    hemi = atlstruct.hemi
    zone = atlstruct.zone
    
    print('   Writing ATL03 .las file...', end = " ")
    try:
        outname = atlstruct.atl03FileName + '_' + atlstruct.gtNum + '.las'
    except AttributeError:
        # occasionally atl03FileName is not in the atlstruct
        outname = atlstruct.atlFileName + '_' + atlstruct.gtNum + '.las'

    outfile = os.path.normpath(outpath + '/' + outname)
    
    if(not os.path.exists(os.path.normpath(outpath))):
        os.mkdir(os.path.normpath(outpath))
    
    # Get projection
    if(atlstruct.zone=='3413' or atlstruct.zone=='3976'):
        # 3413 = Arctic, 3976 = Antarctic
        lasProjection = atlstruct.hemi
        # Write .las file
        writeLas(xx,yy,zz,lasProjection,outfile,cc,ii,sigconf)

    else:
        # Write .las file for UTM projection case
        writeLas(xx,yy,zz,'utm',outfile,cc,ii,sigconf,hemi,zone)

    print('Complete') 

def get_H5_keys_info(atl08filepath,gt):
    keys = getH5Keys(atl08filepath, gt)
    h = h5py.File(atl08filepath, 'r')
    key_name = []
    key_type = []
    key_len = []
    for key in keys:
        try:
            data = h[gt + '/' + key]
            kname = str(key)
            ktype = str(data.dtype)
            klen = int(len(data))
            key_name.append(kname)
            key_type.append(ktype)
            key_len.append(klen)
        except:
            kname = str(key)
            ktype = 'Group'
            klen = 0
            key_name.append(kname)
            key_type.append(ktype)
            key_len.append(klen)
    key_info = [list(a) for a in zip(key_name, key_type, key_len)]
    return key_info

def match_atl_to_atl03(df, atl03struct):
    # Calculate Time
    delta_time03 = np.array(atl03struct.df.delta_time)
    time = np.array(df.delta_time) -\
        np.min(delta_time03[np.nonzero(delta_time03)])
    df = pd.concat([df,pd.DataFrame(time,
                            columns=['time'])],axis=1)
    # Calculate Projected Coordinates
    df, epsg = get_atl_coords(df, atl03struct.epsg)
    # Calculate Along track
    df, rotation_data = get_atl_alongtrack(df, atl03struct)        
    return df, rotation_data, epsg
        
def get_atl03_struct(atl03filepath, gt, atl08filepath = None, epsg = None, 
                     kml_bounds_txt = None, header_file_path = None):
    df = read_atl03_heights_data(atl03filepath, gt)
    if atl08filepath:
        try:
            # df = merge_label_to_df(atl03filepath, atl08filepath, gt, df)
            df['classification'] = get_atl03_classification(atl03filepath,\
                                                            atl08filepath, gt)
            # df = merge_norm_h_to_df(atl03filepath, atl08filepath, gt, df)
            df['norm_h'] = get_atl03_heights_offset(atl03filepath,\
                                                    atl08filepath, gt)
            dataIsMapped = True
        except:
            dataIsMapped = False
    else:
        dataIsMapped = False
    df['seg_id'] =  get_atl03_segment_id(atl03filepath, gt)
    df['ph_bihr'], df['ph_bcr'], df['ph_rate'] =\
        get_atl03_rate(atl03filepath, gt) 
    df['dist_ph_along'] = get_atl03_dist_ph_along(atl03filepath, gt)
    df['time'] = df['delta_time'] -\
        np.min(np.array(df.delta_time)[np.nonzero(np.array(df.delta_time))])
    df = get_atl_time(df)
    df, epsg = get_atl_coords(df, epsg)
    df, rotation_data = get_atl_alongtrack(df)
    track_direction = get_direction(np.array(df.lat_ph))
    atl03filename = get_file_name(atl03filepath)
    atl03_info = getNameParts(atl03filename)
    hemi, zone = identify_hemi_zone(epsg)

    if dataIsMapped:
        # for some reason, need to reset nan to -1,
        # evemn though it's done in get_atl03_classification..
        c = np.array(df.classification)
        nan_index = np.where(np.isnan(c))
        c[nan_index] = -1 # assign nan to -1
        c = c.astype(int)
        df.classification = c
        # df.replace({'classification' : np.nan}, -1)

    # Assign everything to the struct
    beamNum = GtToBeamNum(atl03filepath, gt)
    beamStrength = GtToBeamSW(atl03filepath, gt)
    atl03Struct = AtlStruct(df, gt, beamNum, beamStrength, epsg, zone, hemi,
                              atl03filepath, atl03filename, track_direction,
                              'ATL03', atl03_info, dataIsMapped,rotation_data)

        
    return atl03Struct

def get_atl08_struct(atl08filepath, gt, atl03struct = None, epsg = None, 
                     kml_bounds_txt = None):
    df = read_atl08_land_segments(atl08filepath, gt)
    
    # If ATL03 Struct is available
    if atl03struct:
        # Calculate Time
        df, rotation_data, epsg = match_atl_to_atl03(df, atl03struct) 

    # If ATL03 Struct is not available
    else:
        # Calculate Time
        df = get_atl_time(df)
        # Calculate Projected Coordinates
        df, epsg = get_atl_coords(df, epsg)    
        # Calculate Along Track    
        df, rotation_data = get_atl_alongtrack(df)
    track_direction = get_direction(np.array(df.latitude))
    atl08filename = get_file_name(atl08filepath)
    atl08_info = getNameParts(atl08filename)
    hemi, zone = identify_hemi_zone(epsg)
    dataIsMapped = True
    beamNum = GtToBeamNum(atl08filepath, gt)
    beamStrength = GtToBeamSW(atl08filepath, gt)

    # Assign everything to the struct
    atl08Struct = AtlStruct(df, gt, beamNum, beamStrength, epsg, zone, hemi,
                              atl08filepath, atl08filename, track_direction,
                              'ATL08', atl08_info, dataIsMapped,rotation_data)
    return atl08Struct

def get_geolocation_mapping(height_len, ph_index_beg, segment_ph_cnt, target):
    data = np.zeros(height_len)
    for i_id in range(0,len(target)):
        data[ph_index_beg[i_id]-1:\
             ph_index_beg[i_id]-1 + segment_ph_cnt[i_id]] =\
             np.full((segment_ph_cnt[i_id]), target[i_id]) 
    return data

def append_atl03_geolocation(heights, geolocation, fields = ['segment_id']):
    height_len = len(heights)
    ph_index_beg = np.array(geolocation.ph_index_beg)
    segment_ph_cnt = np.array(geolocation.segment_ph_cnt)
    for field in fields:
        target = np.array(geolocation[field])
        data = get_geolocation_mapping(height_len, ph_index_beg, 
                                       segment_ph_cnt, target)
        heights = pd.concat([heights,pd.DataFrame(data,
                                columns=[field])],axis=1)
    
    return heights
    
# def convert_atl03_to_legacy(atl03):
#     intensity = np.zeros(len(atl03.df))
#     atl03_zMsl = np.zeros(len(atl03.df))
#     segmentID = np.zeros(len(atl03.df))
#     atl03h5Info = getNameParts(atl03.atlFileName)
#     atl03legacy = Atl03StructLegacy(atl03.df.lat_ph, atl03.df.lon_ph, 
#                               atl03.df.easting, atl03.df.northing, 
#                               atl03.df.crosstrack, atl03.df.alongtrack, 
#                               atl03.df.h_ph, atl03_zMsl, atl03.df.time, 
#                               atl03.df.delta_time, atl03.df.signal_conf_ph, 
#                               atl03.df.signal_conf_ph, atl03.df.signal_conf_ph,
#                               atl03.df.signal_conf_ph, atl03.df.signal_conf_ph,
#                               atl03.df.classification, intensity, intensity,
#                               segmentID, 
#                               atl03.gtNum, atl03.beamNum, atl03.beamStrength,
#                               atl03.zone, atl03.hemi, atl03.atlFilePath, 
#                               atl03.atlFileName, atl03.trackDirection, 
#                               atl03h5Info, 
#                               atl03.dataIsMapped)
#     rotationData = atl03.rotationData
#     # headerData = atl03.headerData
#     return atl03legacy, rotationData

def write_pickle(data, filename):
    import pickle
    fp = open(filename, 'wb')
    pickle.dump(data, fp)
    fp.close()

def read_pickle(filename):
    import pickle
    fp = open(filename, 'rb')
    data = pickle.load(fp)
    fp.close()
    return data

def convert_df_to_mat(df,outfilename):
    from scipy import io
    comps =  outfilename.split('.')
    if comps[-1] != 'mat':
        outfilename = outfilename + ".mat"
    # scipy.io.savemat(outfilename, {'struct':df.to_dict("list")})
    io.savemat(outfilename, {'struct':df.to_dict("list")})
    
def get_attribute_info(atlfilepath, gt):
    # add year/doy, sc_orient, beam_number/type to 08 dataframe
    year, doy = get_h5_meta(atlfilepath, meta='date', rtn_doy=True)

    with h5py.File(atlfilepath, 'r') as fp:
        try:
            fp_a = fp[gt].attrs
            description = (fp_a['Description']).decode()
            beam_type = (fp_a['atlas_beam_type']).decode()
            atlas_pce = (fp_a['atlas_pce']).decode()
            spot_number = (fp_a['atlas_spot_number']).decode()
            atmosphere_profile = (fp_a['atmosphere_profile']).decode()
            groundtrack_id = (fp_a['groundtrack_id']).decode().lower()
            sc_orient = (fp_a['sc_orientation']).decode().lower()
        except:
            description = ''
            beam_type = ''
            atlas_pce = ''
            spot_number = ''
            atmosphere_profile = ''
            groundtrack_id = ''
            sc_orient = ''
    info_dict = {
        "description" : description,
        "atlas_beam_type" : beam_type,
        "atlas_pce" : atlas_pce,
        "atlas_spot_number" : spot_number,
        'atmosphere_profile' : atmosphere_profile,
        "groundtrack_id" : groundtrack_id,
        "sc_orientation" : sc_orient,
        "year" : year,
        "doy" : doy

        }
    
    return info_dict
    
# if __name__ == "__main__":    
    # if os.name == 'nt':
    #     basepath = 'Y:/USERS/eric/2_production/'
    # else:
    #     basepath = '/LIDAR/server/USERS/eric/2_production/'

    # atl03file = 'ATL03_20181021130238_03500103_002_01.h5'
    # atl08file = 'ATL08_20181021130238_03500103_002_01.h5'
    # atl09file = 'ATL09_20181016132105_02740101_002_01.h5'
    # # Inputs
    # atl03filepath = basepath + atl03file
    # atl08filepath = basepath + atl08file
    # atl09filepath = basepath + atl09file
    # gt = 'gt1r'
    

    # print('ATL03 Heights')
    # atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath)

    # print('ATL03 Geolocation')
    # geolocation = read_atl03_geolocation(atl03filepath, gt)

    # print('ATL08 Land Segments')
    # atl08 = get_atl08_struct(atl08filepath, gt, atl03)

    # print('ATL09 High Frequency Data')    
    # atl09 = get_atl09_struct(atl09filepath, gt, atl03)

    # print('Append ATL03 Heights Dataframe with Geolocation')
    # heights = atl03.df
    # heights = append_atl03_geolocation(heights, geolocation, 
    #                                    fields = ['segment_id'])
    
    # print('Convert ATL03 Struct to Legacy ATL03 Struct')
    # atl03legacy, rotationData = convert_atl03_to_legacy(atl03)