#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:26:23 2020

@author: eguenther
"""
import os
import time
import numpy as np
import pandas as pd
from icesatReader import get_atl03_struct
from icesatReader import convert_atl03_to_legacy
from icesatReader import get_atl_alongtrack
from icesatReference import estimate_segment_id_legacy
from icesatReference import legacy_get_meas_error
from icesatReader import read_atl03_geolocation
from icesatReader import match_atl_to_atl03
from icesatReader import append_atl03_geolocation
from icesatUtils import superFilter
from icesatCalVal import perfectClassifier
from icesatIO import read_geotiff
from icesatIO import find_intersecting_values
from icesatReader import get_atl09_struct
from icesatReader import get_atl08_struct
from scipy import stats
import scipy
from getAtlTruthSwath_auto import getAtlTruthSwath as getAtlTruth
import matplotlib.pyplot as plt
import pickle as pl

def convert_df_to_mat(df,outfilename):
    comps =  outfilename.split('.')
    if comps[-1] != 'mat':
        outfilename = outfilename + ".mat"
    scipy.io.savemat(outfilename, {'struct':df.to_dict("list")})

def get_max98(series):
    max98 = np.percentile(series, 98)
    return max98

def get_len(series):
    length = int(len(series))
    return length

def get_len_unique(series):
    try:
        length = int(len(np.unique(series)))
    except:
        length = np.nan
    return length

def get_mode(series):
    try:
        length = stats.mode(series)
    except:
        lemgth = np.nan
    return length

def calculate_seg_meteric(df_in, df_out, classification, operation, field, 
                          outfield, classfield = 'classification'):
    df_filter = df_in[df_in[classfield].isin(classification)]
    zgroup = df_filter.groupby('segment_id_beg')
    zout = zgroup.aggregate(operation)
    zout['segment_id_beg'] = zout.index
    zout = zout.reset_index(drop = True)
    # zout['segment_id_beg'] = zout['seg_id']
    zout[outfield] = zout[field]
    zout = zout.filter([outfield,'segment_id_beg'])
    df_out = df_out.merge(zout, on="segment_id_beg",how='left')  
    return df_out

def parse_times(min_time, max_time, seg):
    diff = max_time - min_time
    rem = diff % seg
    div = diff / seg  
    min_time_list = []
    max_time_list = []
    
    if diff > seg: 
        if max_time - min_time >= (2 * seg):
            for i in range(0,int(div) - 1):
                min_time_list.append(min_time + (seg * i))
                max_time_list.append(min_time + (seg * i) + seg)
        
        if rem < 5:
            min_time_list.append(min_time + (seg * (int(div) - 1)))
            max_time_list.append(max_time)
        if rem >= 5:

            min_time_list.append(min_time + (seg * (int(div) - 1)))
            seg =  ((max_time - min_time) - (seg * (int(div) - 1))) / 2
            max_time_list.append(min_time_list[-1] + seg)
            min_time_list.append(max_time_list[-1])
            max_time_list.append(max_time)
    else:
        min_time_list.append(min_time)
        max_time_list.append(max_time)
    
    return min_time_list, max_time_list

def generate_atl03_truth_plot(atl03, outfolder, filename, df_truth):
    uy = atl03.df.alongtrack[atl03.df.classification == -1]
    uz = atl03.df.h_ph[atl03.df.classification == -1]
    dy = atl03.df.alongtrack[atl03.df.classification == 0]
    dz = atl03.df.h_ph[atl03.df.classification == 0]
    gy = atl03.df.alongtrack[atl03.df.classification == 1]
    gz = atl03.df.h_ph[atl03.df.classification == 1]
    cy = atl03.df.alongtrack[atl03.df.classification == 2]
    cz = atl03.df.h_ph[atl03.df.classification == 2]
    hy = atl03.df.alongtrack[atl03.df.classification == 3]
    hz = atl03.df.h_ph[atl03.df.classification == 3]

    tgy = df_truth.alongtrack[df_truth.classification == 2]
    tgz = df_truth.z[df_truth.classification == 2]
    tcy = df_truth.alongtrack[df_truth.classification == 4]
    tcz = df_truth.z[df_truth.classification == 4]
    
    f = plt.figure()
    plt.plot(uy, uz, '.', color = [0.8,0.8,1],
             label = 'ATL03 Unclassified') 
    plt.plot(dy, dz, '.', color = [0,0.5,0.8],
             label = 'ATL03 DRAGANN') 
    plt.plot(tcy[::100], tcz[::100], '.', color = [0.6,0.6,0.6],
             label = 'Truth Canopy') 
    plt.plot(tgy[::100], tgz[::100], '.', color = [0.3,0.3,0.3],
             label = 'Truth Ground') 
    plt.plot(hy, hz, '.', color = [0,0.8,0],
             label = 'ATL03 High Canopy') 
    plt.plot(cy, cz, '.', color = [0,0.5,0],
             label = 'ATL03 Canopy')     
    plt.plot(gy, gz, '.', color = [0.8,0.5,0],
             label = 'ATL03 Ground')    
    plt.legend()

    outfile_png = outfolder + '/graphs/atl03/png/' + filename + '.png'
    outfile_pkl = outfolder + '/graphs/atl03/pkl/' + filename + '.pkl'

    plt.title(atl03.atlFileName)
    plt.ylabel('Height (m)')
    plt.xlabel('Along-track (m)')
    pl.dump(f,open(str(outfile_pkl), 'wb'))

    plt.savefig(outfile_png)
    plt.close(f)
    
def generate_atl03_plot(atl03, outfolder, filename):
    uy = atl03.df.alongtrack[atl03.df.classification == -1]
    uz = atl03.df.h_ph[atl03.df.classification == -1]
    dy = atl03.df.alongtrack[atl03.df.classification == 0]
    dz = atl03.df.h_ph[atl03.df.classification == 0]
    gy = atl03.df.alongtrack[atl03.df.classification == 1]
    gz = atl03.df.h_ph[atl03.df.classification == 1]
    cy = atl03.df.alongtrack[atl03.df.classification == 2]
    cz = atl03.df.h_ph[atl03.df.classification == 2]
    hy = atl03.df.alongtrack[atl03.df.classification == 3]
    hz = atl03.df.h_ph[atl03.df.classification == 3]  
    f = plt.figure()
    plt.plot(uy, uz, '.', color = [0.8,0.8,1],
             label = 'ATL03 Unclassified') 
    plt.plot(dy, dz, '.', color = [0,0.5,0.8],
             label = 'ATL03 DRAGANN') 
    plt.plot(hy, hz, '.', color = [0,0.8,0],
             label = 'ATL03 High Canopy') 
    plt.plot(cy, cz, '.', color = [0,0.5,0],
             label = 'ATL03 Canopy')     
    plt.plot(gy, gz, '.', color = [0.8,0.5,0],
             label = 'ATL03 Ground')      
    plt.legend()

    outfile_png = outfolder + '/graphs/atl03/png/' + filename + '.png'
    outfile_pkl = outfolder + '/graphs/atl03/pkl/' + filename + '.pkl'

    plt.title(atl03.atlFileName)
    plt.ylabel('Height (m)')
    plt.xlabel('Along-track (m)')
    pl.dump(f,open(str(outfile_pkl), 'wb'))

    plt.savefig(outfile_png)
    plt.close(f)
    
def generate_atl08_truth_plot(atl03, atl08, outfolder, filename, df_truth):
    gbfy = atl08.df.alongtrack + atl08.df.alongtrackoffset
    gbfz = atl08.df.h_te_best_fit + atl08.df.zoffset
    gbfz[gbfz > 10000] = np.nan
    
    gmy = atl08.df.alongtrack + atl08.df.alongtrackoffset
    gmz = atl08.df.h_te_median + atl08.df.zoffset
    gmz[gmz > 10000] = np.nan
    
    cmy = atl08.df.alongtrack + atl08.df.alongtrackoffset
    cmz = atl08.df.h_canopy_abs + atl08.df.zoffset
    cmz[cmz > 10000] = np.nan
    
    #ATL03 median ground
    agy = atl08.df.alongtrack
    agz = atl08.df.atl03_ground_median
    
    #ATL03 canopy
    acy = atl08.df.alongtrack
    acz = atl08.df.atl03_canopy_max98
    #PC median ground
    #ATL03 canopy
    pgy = atl08.df.alongtrack
    pgz = atl08.df.pc_ground_median
    #PC max canopy
    pcy = atl08.df.alongtrack
    pcz = atl08.df.pc_canopy_max98
    
    #Truth Ground
    tgy = atl08.df.alongtrack
    tgz = atl08.df.truth_ground_median
    #Truth Canopy
    tcy = atl08.df.alongtrack
    tcz = atl08.df.truth_canopy_max98   
    
    
    tpgy = df_truth.alongtrack[df_truth.classification == 2]
    tpgz = df_truth.z[df_truth.classification == 2]
    tpcy = df_truth.alongtrack[df_truth.classification == 4]
    tpcz = df_truth.z[df_truth.classification == 4]
    
    f = plt.figure()
    ratio = int(np.ceil(len(tpcy)/5000))
    plt.plot(tpcy[::ratio], tpcz[::ratio], '.', color = [0.6,0.6,0.6],
             label = 'Truth Canopy')
    plt.plot(tpgy[::ratio], tpgz[::ratio], '.', color = [0.3,0.3,0.3],
             label = 'Truth Ground')
    
    plt.plot(cmy, cmz, '.', color = [0,0.8,0],
             label = 'ATL08 Canopy')
    plt.plot(acy, acz, '.', color = [0,0.5,0],
             label = 'ATL03 Canopy')    
    plt.plot(pcy, pcz, '.', color = [0.1,0.8,0.5],
             label = 'PC Canopy')        
    plt.plot(tcy, tcz, '.', color = [0.1,0.3,0.7],
             label = 'Truth Canopy')    
    
    plt.plot(gbfy, gbfz, '.', color = [0.7,0.1,0.7],
             label = 'ATL08 Ground Best Fit')    
    plt.plot(gmy, gmz, '.', color = [0.8,0.3,0.1],
             label = 'ATL08 Ground Median')    
    plt.plot(agy, agz, '.', color = [0.7,0.5,0.1],
             label = 'ATL03 Ground Median')        
    plt.plot(pgy, pgz, '.', color = [0.1,0.4,0.6],
             label = 'Perfect Classifier Ground Median')     
    plt.plot(tgy, tgz, '.', color = [0.01,0.2,0.5],
             label = 'Truth Ground Median') 
    
    plt.legend()
    outfile_png = outfolder + '/graphs/atl08/png/' + filename + '.png'
    #\graphs\alt08\pkl\
    outfile_pkl = outfolder + '/graphs/atl08/pkl/' + filename + '.pkl'

    plt.title(atl03.atlFileName)
    plt.ylabel('Height (m)')
    plt.xlabel('Along-track (m)')
    pl.dump(f,open(str(outfile_pkl), 'wb'))

    plt.savefig(outfile_png)
    plt.close(f)
    
def generate_atl08_plot(atl03, atl08, outfolder, filename):
    gbfy = atl08.df.alongtrack + atl08.df.alongtrackoffset
    gbfz = atl08.df.h_te_best_fit + atl08.df.zoffset
    gbfz[gbfz > 10000] = np.nan
    
    gmy = atl08.df.alongtrack + atl08.df.alongtrackoffset
    gmz = atl08.df.h_te_median + atl08.df.zoffset
    gmz[gmz > 10000] = np.nan
    
    cmy = atl08.df.alongtrack + atl08.df.alongtrackoffset
    cmz = atl08.df.h_canopy_abs + atl08.df.zoffset
    cmz[cmz > 10000] = np.nan
    
    #ATL03 median ground
    agy = atl08.df.alongtrack
    agz = atl08.df.atl03_ground_median
    
    #ATL03 canopy
    acy = atl08.df.alongtrack
    acz = atl08.df.atl03_canopy_max98

    f = plt.figure()

    
    plt.plot(cmy, cmz, '.', color = [0,0.8,0],
             label = 'ATL08 Canopy')
    plt.plot(acy, acz, '.', color = [0,0.5,0],
             label = 'ATL03 Canopy')     

    
    plt.plot(gbfy, gbfz, '.', color = [0.7,0.1,0.7],
             label = 'ATL08 Ground Best Fit')    
    plt.plot(gmy, gmz, '.', color = [0.8,0.3,0.1],
             label = 'ATL08 Ground Median')    
    plt.plot(agy, agz, '.', color = [0.7,0.5,0.1],
             label = 'ATL03 Ground Median')    

    plt.legend()

    outfile_png = outfolder + '/graphs/atl08/png/' + filename + '.png'
    outfile_pkl = outfolder + '/graphs/atl08/pkl/' + filename + '.pkl'

    plt.title(atl03.atlFileName)
    plt.ylabel('Height (m)')
    plt.xlabel('Along-track (m)')
    pl.dump(f,open(str(outfile_pkl), 'wb'))

    plt.savefig(outfile_png)
    plt.close(f)
    
    
def segment_analysis(header_file_path, kml_bounds_txt1, truthSwathDir,
                     outFilePath_truth, outFilePath_corrected, epsg_str,
                     atl03filepath, atl08filepath, gt, min_time, 
                     max_time, outfolder):


    # header_file_path =\
    #     '/LIDAR/server/USERS/eric/1_experiment/Finland_HeaderData.mat'
        
    # kml_bounds_txt1 = '/LIDAR/server/USERS/eric/2_production/kmlBounds.txt'
    # epsg_string = '32635'
    # Read ATL03 Struct
    print('Generate ATL03 Struct')
    atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath, 
                             epsg = epsg_string, 
                             kml_bounds_txt = kml_bounds_txt1, 
                             header_file_path = header_file_path)    
    
    # Read ATL03 Geolocation Subgroup as DF
    print('Read Geolocation Subgroup')
    geolocation = read_atl03_geolocation(atl03filepath, gt)
    atl03.df = append_atl03_geolocation(atl03.df, geolocation, 
                                        fields = ['segment_id'])
    geolocation, gl_rotation_data, gl_epsg = match_atl_to_atl03(geolocation, 
                                                                atl03)

    # Trim Data by Time
    print('Trim Data by Time')    
    atl03.df = atl03.df[atl03.df['time'] < max_time]
    geolocation = geolocation[geolocation['time'] < max_time]

    atl03.df = atl03.df[atl03.df['time'] > min_time]
    geolocation = geolocation[geolocation['time'] > min_time]
    
    atl03.df = atl03.df.reset_index(drop = True)
    geolocation = geolocation.reset_index(drop = True)


    # Recalculate alongtrack/crosstrack for shortened granule
    atl03.df, rotation_data = get_atl_alongtrack(atl03.df)    
    atl03.rotationData = rotation_data
    
    # "Match" the geolocation df to the ATL03 struct
    print('Match Geolocation to ATL03')    

    geolocation, gl_rotation_data, gl_epsg = match_atl_to_atl03(geolocation, 
                                                                atl03)
    # Convert the ATL03 Struct to the legacy ATL03 Struct
    print('Convert Struct to Legacy')    
    atl03legacy, rotationData, headerData = convert_atl03_to_legacy(atl03)
    
    # Legacy Truth Swath Inputs
    buffer = 50   # Distance in cross-track (meters) around ATL03 track
    useExistingTruth = False # Option to use existing truth data if it exists
    # truthSwathDir = '/laserpewpew/data/validation/data/Finland/LAS_UTM'
    # outFilePath_truth = '/LIDAR/server/USERS/eric/1_experiment/finland_analysis/' +\
    #     'las/truth'
    createTruthFile = True      # Option to create output truth .las file

    
    # Call getAtlTruthSwath (with ACE)
    print('Run Legacy Truth Swath')
    try:
        timeStart = time.time()
        atlTruthData = getAtlTruth(atl03legacy, headerData, 
                                      rotationData, useExistingTruth, 
                                      truthSwathDir, buffer, outFilePath_truth, 
                                      createTruthFile)

        atlTruthData.classification[atlTruthData.classification == 3] = 4
        atlTruthData.classification[atlTruthData.classification == 5] = 4

        
    
        timeEnd = time.time()
        timeElapsedTotal = timeEnd - timeStart
        timeElapsedMin = np.floor(timeElapsedTotal / 60)
        timeElapsedSec = timeElapsedTotal % 60
        print('   Script Completed in %d min %d sec.' % (timeElapsedMin, 
                                                          timeElapsedSec))
        print('\n')
        
        # outFilePath_corrected = '/LIDAR/server/USERS/eric/1_experiment/' +\
        #     'finland_analysis/las/truth_corrected'
        
        print('Run Legacy Offset Coorection')    
        atlCorrections = legacy_get_meas_error(atl03legacy, atlTruthData, 
                                                rotationData, 
                                                outFilePath_corrected)
        
        # Apply ATLCorrections to the Geolocation
        geolocation.alongtrack = geolocation.alongtrack +\
            atlCorrections.alongTrack
        geolocation.crosstrack = geolocation.crosstrack +\
            atlCorrections.crossTrack
            
        # Apply ATLCorrectuibs to ATL03 Legacy
        atl03legacy.alongTrack = atl03legacy.alongTrack +\
            atlCorrections.alongTrack
        atl03legacy.crossTrack = atl03legacy.crossTrack +\
            atlCorrections.crossTrack
        atl03legacy.z = atl03legacy.z +\
            atlCorrections.z
    
        # Apply ATLCorrections to ATL03 DF        
        atl03.df.alongtrack = atl03.df.alongtrack +\
            atlCorrections.alongTrack
        atl03.df.crosstrack = atl03.df.crosstrack +\
            atlCorrections.crossTrack
        atl03.df.h_ph = atl03.df.h_ph +\
            atlCorrections.z
           
        # Run Superfilter Legacy
        superTruth, sortedMeasured = superFilter(atl03legacy, atlTruthData, 
                                                 xBuf = 5.5,classCode = [])
        
        # Run Perfect Classifier Legacy and assign to ATL03\
        truthgroundclass = 2
        truthcanopyclass = [3,4,5]
        unclassedlist = [6, 9, 13, 18]
        measpc, measoc = perfectClassifier(sortedMeasured, superTruth,
                                           ground = [truthgroundclass],
                                           canopy = truthcanopyclass, 
                                           unclassed = unclassedlist, 
                                           keepsize = True)
        
        # Sort ATL03 by Along Track
        
        ##Create new column for Index1
        print('Sort Alongtrack')
        atl03.df['index1'] = atl03.df.index
        
        ##Sort by Along-track
        atl03.df = atl03.df.sort_values(by=['alongtrack'])
        
        ##Reset Index
        atl03.df = atl03.df.reset_index(drop = True)
        
        ##Join PC to DF
        atl03.df = pd.concat([atl03.df,pd.DataFrame(
            measpc,columns=['perfect_class'])],axis=1)
        atl03.df = pd.concat([atl03.df,pd.DataFrame(
            measoc,columns=['generic_class'])],axis=1)
        
        ##Sort by Index1
        atl03.df = atl03.df.sort_values(by=['index1'])
        
        ##Reset Index
        atl03.df = atl03.df.reset_index(drop = True)
       
        ##Remove Index1
        atl03.df = atl03.df.drop(columns=['index1'])
    
        # Read ATL09
        # atl09 = get_atl09_struct(atl09filepath, gt, atl03)
        
        # df_seg = geolocation.merge(atl09.df, on="segment_id",how='left')
        # df_seg = df_seg.fillna(method='ffill',limit=14)
        
        # Assign Segment ID Values to Tuth Data
        seg_id_truth, include = estimate_segment_id_legacy(geolocation, gt, 
                                                           superTruth)
    
        # Calculate Error and Meterics
        
        ##Filter Truth Data by Include
        alongtrack = superTruth.alongTrack.flatten()
        crosstrack = superTruth.crossTrack.flatten()
        z = superTruth.z.flatten()
        easting = superTruth.easting.flatten()
        northing = superTruth.northing.flatten()
        classification = superTruth.classification.flatten()
        intensity = superTruth.intensity.flatten()
        
        alongtrack = alongtrack[include == 1]
        crosstrack = crosstrack[include == 1]
        z = z[include == 1]
        easting = easting[include == 1]
        northing = northing[include == 1]
        classification = classification[include == 1]
        intensity = intensity[include == 1]
        seg_id_truth = seg_id_truth[include == 1]
        truth_flag = True
    except:
        print('Truth Failed, continue with rest of code')
        atl03.df['perfect_class'] = np.nan
        atl03.df['generic_class'] = np.nan
      
        
        truth_flag = False
    
    ##Create ATLTruth DF
    if truth_flag == True:
        df_truth = pd.DataFrame(z,columns=['z'])
        df_truth = pd.concat([df_truth,pd.DataFrame(
            crosstrack,columns=['crosstrack'])],axis=1)
        df_truth = pd.concat([df_truth,pd.DataFrame(
            alongtrack,columns=['alongtrack'])],axis=1)
        df_truth = pd.concat([df_truth,pd.DataFrame(
            easting,columns=['easting'])],axis=1)
        df_truth = pd.concat([df_truth,pd.DataFrame(
            northing,columns=['northing'])],axis=1)
        df_truth = pd.concat([df_truth,pd.DataFrame(
            classification,columns=['classification'])],axis=1)
        # df_truth = pd.concat([df_truth,pd.DataFrame(
        # intensity,columns=['intensity'])],axis=1)
        df_truth = pd.concat([df_truth,pd.DataFrame(
            seg_id_truth,columns=['segment_id'])],axis=1)

    ##Find ATL08 Segment Range
    
    ###Read ATL08
    atl08 = get_atl08_struct(atl08filepath, gt, atl03)


    atl08.df = atl08.df[atl08.df['time'] <= (max_time - min_time)]
    atl08.df = atl08.df[atl08.df['time'] >= 0]    
    atl08.df = atl08.df.reset_index(drop = True)


    
    atl08.df, atl08_rotation_data, atl08_epsg = match_atl_to_atl03(atl08.df, 
                                                                atl03)
    
    atl03.df = atl03.df.reset_index(drop = True)
    geolocation = geolocation.reset_index(drop = True)

    ###Get ATL08 Keys    
    atl08_seg = np.array(atl08.df.segment_id_beg)
    seg_id = np.array(geolocation.segment_id)
    atl08_key_df = pd.DataFrame(atl08_seg,columns=['segment_id'])
    atl08_key_df = pd.concat([atl08_key_df,pd.DataFrame(
        atl08_seg,columns=['segment_id_beg'])],axis=1)
    atl08_key_df = pd.concat([atl08_key_df,pd.DataFrame(
        atl08_seg,columns=['seg_id'])],axis=1)
    
    key_df = pd.DataFrame(seg_id,columns=['segment_id'])
    key_df = key_df.merge(atl08_key_df, on="segment_id",how='left')  
    key_df = key_df.fillna(method='ffill',limit=4)
    
    max_seg = max(geolocation.segment_id)
    min_seg = min(geolocation.segment_id)

    key_df = key_df[key_df['segment_id'] <= max_seg]
    key_df = key_df[key_df['segment_id'] >= min_seg]
    
    

    ###Merge Geolocation/ATL09
    # df_atl09 = df_seg.merge(key_df, on="segment_id", how="left")
    

    ###Merge ATL08 Keys to Truth
    if truth_flag == True:
        df_truth = df_truth.merge(key_df, on="segment_id",how="left")
    
    
    ###Merge ATL08 Keys to ATL03    
    df_atl03 = atl03.df.merge(key_df, on="segment_id",how="left")

    # Calculate Meterics
    
    # Assign Geolocation/ATL09

    # zgroup = df_atl09.groupby('segment_id_beg')
    # zout = zgroup.aggregate(pd.Series.mode)
    # zout = zgroup.aggregate(np.median)
    # zout['segment_id_beg'] = zout.index
    # zout = zout.reset_index(drop = True)
    # zout['segment_id_beg'] = zout['seg_id']
    # df_out = df_out.merge(zout, on="segment_id_beg",how='left')  
    # return df_out    
    # zout = zout.reset_index().T.drop_duplicates().T
    # atl09_list = ['aclr_true','apparent_surf_reflec','backg_c',
    #               'backg_theoret','beam_azimuth','beam_elevation','segment_id_beg']
    # zout = zout.filter(['segment_id_beg'])
    
    # zout2 = zout2.filter([outfield,'segment_id_beg'])
    # df_test3 = atl08.df.merge(zout, on="segment_id_beg",how='left')  
    
    # Assign Tif
    ##Cornie
    cornie = '/LIDAR/server/USERS/eric/1_experiment/global_products/' +\
        'Corine_LandCover_europe/cornie_landcover_finland_UTM.tif'
    data, epsg, ulx, uly, resx, resy = read_geotiff(cornie)
    
    x = np.array(atl08.df.easting)
    y = np.array(atl08.df.northing)
    result = find_intersecting_values(x,y,data,ulx,uly,resx,resy)
    
    atl08.df['Corine_LC'] = result
    
    ##Forest Canopy Height
    simard = '/LIDAR/server/USERS/eric/1_experiment/global_products/' +\
        'Corine_LandCover_europe/Simard_Forest_Height_Finland_UTM_auto.tif'
    data, epsg, ulx, uly, resx, resy = read_geotiff(simard)
    
    x = np.array(atl08.df.easting)
    y = np.array(atl08.df.northing)
    result = find_intersecting_values(x,y,data,ulx,uly,resx,resy)
    
    atl08.df['Simard_Forest_Height'] = result
    
    ##Truth Meterics
    ###Truth Ground Median
    print('Apply meterics')
    if truth_flag == True:
        atl08.df = calculate_seg_meteric(df_truth, atl08.df, [2], np.median,'z',
                                         'truth_ground_median')
        
        ###Truth Canopy Max 98
        atl08.df = calculate_seg_meteric(df_truth, atl08.df, [4], get_max98,'z',
                                         'truth_canopy_max98')
        
        ##Measured Meterics
        ###ATL03 Ground Median
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [1], np.median,'h_ph',
                                         'atl03_ground_median', 
                                         classfield = 'generic_class')
        
        ###ATL03 Canopy Max 98
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [2], get_max98,'h_ph',
                                         'atl03_canopy_max98', 
                                         classfield = 'generic_class')
    else:
        atl08.df['truth_ground_median'] = np.nan
        atl08.df['truth_canopy_max98'] = np.nan
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [1], np.median,'h_ph',
                                         'atl03_ground_median', 
                                         classfield = 'classification')
        
        ###ATL03 Canopy Max 98
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [2,3], get_max98,'h_ph',
                                         'atl03_canopy_max98', 
                                         classfield = 'classification')        
        
  
    ##Perfect Metericsdemocratic primary count
    ###ATL03 Ground Median
    if truth_flag == True:
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [1], np.median,'h_ph',
                                         'pc_ground_median', 
                                         classfield = 'perfect_class')
        
        ###ATL03 Canopy Max 98
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [2], get_max98,'h_ph',
                                         'pc_canopy_max98', 
                                         classfield = 'perfect_class')
        
        ##Truth Size
        atl08.df = calculate_seg_meteric(df_truth, atl08.df, [2], np.size,'z',
                                         'truth_n_ground')
    
        atl08.df = calculate_seg_meteric(df_truth, atl08.df, [4], np.size,'z',
                                         'truth_n_canopy')
        
        atl08.df = calculate_seg_meteric(df_truth, atl08.df, [0], np.size,'z',
                                         'truth_n_unclassed')
    else:
        atl08.df['pc_ground_median'] = np.nan
        atl08.df['pc_canopy_max98'] = np.nan
        atl08.df['truth_n_ground'] = np.nan
        atl08.df['truth_n_unclassed'] = np.nan

    ##ATL03 Size
    atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [-1], np.size,'h_ph',
                                     'atl03_n_unclassified')
      
    atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [0], np.size,'h_ph',
                                     'atl03_n_draggan')

    atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [1], np.size,'h_ph',
                                     'atl03_n_ground')
    
    atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [2], np.size,'h_ph',
                                     'atl03_n_canopy')
    
    atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [3], np.size,'h_ph',
                                     'atl03_n_high_canopy')
    
    ##GC Size
    if truth_flag == True:
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [0], np.size,'h_ph',
                                         'gc_n_unclassified', 
                                         classfield = 'generic_class')
    
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [1], np.size,'h_ph',
                                         'gc_n_ground', 
                                         classfield = 'generic_class')
        
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [2], np.size,'h_ph',
                                         'gc_n_canopy', 
                                         classfield = 'generic_class')
    else:
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [-1,0], np.size,'h_ph',
                                         'gc_n_unclassified')
    
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [1], np.size,'h_ph',
                                         'gc_n_ground')
        
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [2, 3], np.size,'h_ph',
                                         'gc_n_canopy')        

    ##PC Size
    if truth_flag == True:
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [0], np.size,'h_ph',
                                         'pc_n_unclassified', 
                                         classfield = 'perfect_class')
    
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [1], np.size,'h_ph',
                                         'pc_n_ground', 
                                         classfield = 'perfect_class')
        
        atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [2], np.size,'h_ph',
                                         'pc_n_canopy', 
                                         classfield = 'perfect_class')
    else: 
        atl08.df['pc_n_unclassified'] = np.nan
        atl08.df['pc_n_ground'] = np.nan
        atl08.df['pc_n_canopy'] = np.nan
    
    
    ##GC Unique Time
    atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [0], get_len_unique,
                                     'time', 'gc_nshots_unclassed', 
                                     classfield = 'generic_class')

    atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [1], get_len_unique,
                                     'time', 'gc_nshots_ground', 
                                     classfield = 'generic_class')
    
    atl08.df = calculate_seg_meteric(df_atl03, atl08.df, [2], get_len_unique,
                                     'time', 'gc_nshots_canopy', 
                                     classfield = 'generic_class')
 
    if truth_flag == True:
        atl08.df['alongtrackoffset'] = float(atlCorrections.alongTrack)
        atl08.df['crosstrackoffset'] = float(atlCorrections.crossTrack)
        atl08.df['zoffset'] = float(atlCorrections.z)
    else:
        atl08.df['alongtrackoffset'] = np.nan
        atl08.df['crosstrackoffset'] = np.nan
        atl08.df['zoffset'] = np.nan
    
    # ATL03
    outfilename = outfolder + '/mat/atl03/' + atl03.atlFileName + '_' + gt + "_" +\
        str(min_time) + '_' + str(max_time) + '.mat'
    convert_df_to_mat(atl03.df,outfilename)    

    outfilename = outfolder + '/pkl/atl03/' + atl03.atlFileName + '_' + gt + "_" +\
        str(min_time) + '_' + str(max_time) + '.pkl'
    atl03.df.to_pickle(outfilename)

    # ATL08
    outfilename = outfolder + '/mat/atl08/' + atl08.atlFileName + '_' + gt + "_" +\
        str(min_time) + '_' + str(max_time) + '.mat'
    convert_df_to_mat(atl08.df,outfilename)    

    outfilename = outfolder + '/csv/atl08/' + atl08.atlFileName + '_' + gt + "_" +\
        str(min_time) + '_' + str(max_time) + '.csv'
    atl08.df.to_csv(outfilename)

    outfilename = outfolder + '/pkl/atl08/' + atl08.atlFileName + '_' + gt + "_" +\
        str(min_time) + '_' + str(max_time) + '.pkl'
    atl08.df.to_pickle(outfilename)
    
    # Truth
    if truth_flag == True:
        truth_file = atl03.atlFileName.split('ATL03_')[1]
        outfilename = outfolder + '/mat/truth/' + "truth_" + truth_file +\
            '_' + gt + "_" + str(min_time) + '_' + str(max_time) + '.mat'
        convert_df_to_mat(df_truth,outfilename)    

        outfilename = outfolder + '/pkl/truth/' + "truth_" + truth_file +\
            '_' + gt + "_" + str(min_time) + '_' + str(max_time) + '.pkl'
        df_truth.to_pickle(outfilename)    

        outfilename = atl03.atlFileName + '_' + gt + "_" + str(min_time) +\
            '_' + str(max_time)
        generate_atl03_truth_plot(atl03, outfolder, outfilename, df_truth)
        outfilename = atl08.atlFileName + '_' + gt + "_" + str(min_time) +\
            '_' + str(max_time)        
        generate_atl08_truth_plot(atl03, atl08, outfolder, outfilename, df_truth)
    else:
        outfilename = atl03.atlFileName + '_' + gt + "_" + str(min_time) +\
            '_' + str(max_time)
        generate_atl03_plot(atl03, outfolder, outfilename)       
        
        outfilename = atl08.atlFileName + '_' + gt + "_" + str(min_time) +\
            '_' + str(max_time)          
        generate_atl08_plot(atl03, atl08, outfolder, outfilename)


if __name__ == "__main__":
    import csv
    import os

    files_to_read = '/LIDAR/server/USERS/eric/1_experiment/finland_tiles-4.csv'
    outfolder = '/LIDAR/server/USERS/eric/1_experiment/finland_analysis5'

    basepath03 = '/laserpewpew/data/release/002/ATL03_r002/Finland/'
    basepath08 = '/laserpewpew/data/release/002/ATL08_r002/Finland/'
    
    header_file_path = '/LIDAR/server/USERS/eric/1_experiment/' +\
        'Finland_HeaderData.mat'
    
    kml_bounds_txt1 = '/LIDAR/server/USERS/eric/2_production/kmlBounds.txt'

    truthSwathDir = '/laserpewpew/data/validation/data/Finland/LAS_UTM'
    outFilePath_truth = '/LIDAR/server/USERS/eric/1_experiment/' +\
        'finland_analysis/las/truth'
    outFilePath_corrected = '/LIDAR/server/USERS/eric/1_experiment/' +\
        'finland_analysis/las/truth_corrected'
    epsg_string = '32635'
    
    gt_list = ['gt1r','gt1l','gt3r','gt3l','gt2r','gt2l']
    
    release = '002'
    atl03_base_file = []
    start_time = []
    end_time=  []
    atl03filelist = []
    atl08filelist = []

    
    with open(files_to_read) as csvfile:
        readCSV = csv.reader(csvfile,delimiter=',')
        
        for idx, row in enumerate(readCSV):
            if idx > 1:
                for idx2, item in enumerate(row[1::2]):
                    if item:
                        atl03_base_file.append(row[0])
                        start_time.append(row[(idx2*2) + 1])
                        end_time.append(row[(idx2*2) + 2])
                        
    
    for atl03fileshort in atl03_base_file:
        try:
            filecomponent = atl03fileshort.split('ATL03_')[1]
        except:
            filecomponent = atl03fileshort.split('ATL08_')[1]
        atl03fulllist = os.listdir(basepath03)
        atl08fulllist = os.listdir(basepath08)
        atl03filepath = ''
        atl08filepath = ''
        for item in atl03fulllist:
            if (item.split('_' + release + '_')[0] + '_' + release) ==\
                "ATL03_" + filecomponent:
                atl03filepath = item
        for item in atl08fulllist:
            if (item.split('_' + release + '_')[0] + '_' + release) ==\
                "ATL08_" + filecomponent:
                atl08filepath = item
        # atl03filepath = basepath03 + "ATL03_" + filecomponent + "_01.h5"
        # atl08filepath = basepath08 + "ATL08_" + filecomponent +  "_01.h5"
        
        if atl03filepath == '' or atl08filepath == '':
            print('skip')
        else:
            atl03filelist.append(basepath03 + atl03filepath)
            atl08filelist.append(basepath08 + atl08filepath)

    
    for gt in gt_list:
        for i in range(0,len(atl03filelist)):
            atl03filepath = atl03filelist[i]
            atl08filepath = atl08filelist[i]       
            min_t = float(start_time[i])
            max_t = float(end_time[i])
            print(str(min_t))
            print(str(max_t))
            min_time_list, max_time_list = parse_times(min_t, max_t, 10)
            try:
                for j in range(0,len(min_time_list)):
                    min_time = float(min_time_list[j])
                    max_time = float(max_time_list[j])
                    print('!!!!!!!!!!!!!!!!!!!!START!!!!!!!!!!!!!!!!!!!!!!')
                    
                    print(atl03filepath)
                    print(atl08filepath)
                    print(str(min_time))
                    print(str(max_time))
                    
                    outfile = atl08filepath.split(basepath08)[1].split('.h5')[0]
                    outfile = outfile + '_' + gt + '_' + str(min_time) + '_' +\
                        str(max_time) + '.csv'
                    outfile_list = os.listdir(outfolder + '/csv/atl08/')
                    if outfile not in outfile_list:
    
                        # segment_analysis(atl03filepath, atl08filepath, gt, min_time, 
                        #                   max_time, outfolder)
                        segment_analysis(header_file_path, kml_bounds_txt1, 
                                         truthSwathDir, outFilePath_truth, 
                                         outFilePath_corrected, epsg_string,
                                         atl03filepath, atl08filepath, gt, 
                                         min_time, max_time, outfolder)
                        print('!!!!!!!!!!!!!!!!!!!!SUCCESS!!!!!!!!!!!!!!!!!!!!!!')
                    else:
                        print('!!!!!!!!!!!!!!ALREADY PROCESSED!!!!!!!!!!!!!!!!')
                
            except:
                print('!!!!!!!!!!!!!!!!!!!!FAILED!!!!!!!!!!!!!!!!!!!!!!')


    os.system('python ./seg_update.py {} {} {}'.format(outfolder, basepath03, basepath08))