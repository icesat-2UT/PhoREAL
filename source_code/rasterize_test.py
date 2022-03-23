#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 08:26:42 2021

@author: eguenther
"""

import h5py
import os
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from icesatIO import GtToBeamNum, GtToBeamSW
from icesatReader import get_atl08_struct
from icesatReader import get_atl03_struct
from icesatBin import create_atl08_bin


def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i',deltas,deltas)
    return np.argmin(dist_2)

def get_len_unique(series):
    try:
        length = len(np.unique(series))
    except:
        length = np.nan
    return length

def get_len(series):
    try:
        length = len(series)
    except:
        length = np.nan
    return length

def readAtl03DataMapping(in_file03, label, return_delta_time=False):

    if not os.path.isfile(in_file03):
        print('File does not exist')
    
    try:
        f = h5py.File(in_file03, 'r')
    except Exception as e:
        print('Python message: %s\n' % e)
        return [], []

    dsname=label+'/geolocation/segment_id'
    if dsname in f:
        segment_id=np.array(f[dsname])
    else:
        segment_id=[]

    dsname=label+'/geolocation/ph_index_beg'
    if dsname in f:
        ph_index_beg=np.array(f[dsname])
    else:
        ph_index_beg=[]

    if(return_delta_time):
        dsname=label+'/geolocation/delta_time'
        if dsname in f:
            delta_time=np.array(f[dsname])
        else:
            delta_time=[]


    f.close()

    if(return_delta_time):
        return ph_index_beg, segment_id, delta_time
    else:
        return ph_index_beg, segment_id 


##### Function to read ATL08 .h5 files for mapping
def readAtl08DataMapping(in_file08, label):
  if not os.path.isfile(in_file08):
    print('File does not exist')
  try:
    f = h5py.File(in_file08, 'r')
  except Exception as e:
    print('Python message: %s\n' % e)
    return [], [], []

  dsname=label+'/signal_photons/classed_pc_indx'
  if dsname in f:
    classed_pc_indx=np.array(f[dsname])
  else:
    classed_pc_indx=[]

  dsname=label+'/signal_photons/classed_pc_flag'
  if dsname in f:
    classed_pc_flag=np.array(f[dsname])
  else:
    classed_pc_flag=[]

  dsname=label+'/signal_photons/ph_segment_id'
  if dsname in f:
    seg08_id=np.array(f[dsname])
  else:
    seg08_id=[]

  f.close()

  return classed_pc_indx, classed_pc_flag, seg08_id

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

def getAtl08Mapping_seg(atl03_ph_index_beg, atl03_segment_id, atl08_classed_pc_indx, atl08_classed_pc_flag, atl08_segment_id):
      
    # Get ATL03 data
    indsNotZero = atl03_ph_index_beg != 0
    atl03_ph_index_beg = atl03_ph_index_beg[indsNotZero];
    atl03_segment_id = atl03_segment_id[indsNotZero];
    
    # Find ATL08 segments that have ATL03 segments
    atl03SegsIn08TF, atl03SegsIn08Inds = ismember(atl08_segment_id,atl03_segment_id)
    
    # Get ATL08 classed indices and values
    atl08classed_inds = atl08_classed_pc_indx[atl03SegsIn08TF]
    atl08classed_segs = atl08_segment_id[atl03SegsIn08TF]

    # Determine new mapping into ATL03 data
    atl03_ph_beg_inds = atl03SegsIn08Inds;
    atl03_ph_beg_val = atl03_ph_index_beg[atl03_ph_beg_inds];
    newMapping = atl08classed_inds + atl03_ph_beg_val - 2;
    
    # Get max size of output array
    sizeOutput = newMapping[-1]
    
    # Pre-populate all photon classed array with zeroes
    allph_seg = (np.zeros(sizeOutput + 1).astype(int)) - 1
    
    # Populate all photon classed array from ATL08 classifications
    allph_seg[newMapping] = atl08classed_segs
    
    # Return all photon classed array
    return allph_seg

def get_max75(series):
    try:
        max75 = np.percentile(series, 75)
    except:
        max75 = np.nan
    return max75
def get_max25(series):
    try:
        max25 = np.percentile(series, 25)
    except:
        max25 = np.nan
    return max25
def get_max50(series):
    try:
        max50 = np.percentile(series, 50)
    except:
        max50 = np.nan
    return max50

def normalize_heights(df):
    t_ind = np.int32(np.floor((df.delta_time - np.min(df.delta_time)) / 0.001))
    df['t_ind'] = t_ind
    df_g = df[df.c == 1]
    zgroup = df_g.groupby('t_ind')
    zout = zgroup.aggregate(np.median)
    # zout = zout.drop(columns=['h_ind', 'c', 'along', 'delta_time',
    #        'bckgrd_int_height_reduced', 'bckgrd_counts_reduces', 'bckgrd_rate','norm_h'])
    zout = zout.reindex(list(range(0,np.max(t_ind) + 1)))
    zout = zout.interpolate(method='linear', axis=0).ffill().bfill()
    ground = zout.h_ph[df.t_ind]
    norm_height = np.array(df.h_ph) - np.array(ground)
    df['norm_h'] = norm_height
    df = df.drop(columns = ['t_ind'])
    return df


def get_atl03_atl08_seg_keys(atl03filepath, atl08filepath, gt):
    h_seg = get_atl03_segment_id(atl03filepath, gt)
    # Read ATL08 .h5 file
    f08 = h5py.File(atl08filepath,'r')
    seg_beg = np.asarray(f08[gt + '/land_segments/segment_id_beg']).astype(int)
    seg_end = np.asarray(f08[gt + '/land_segments/segment_id_end']).astype(int)
    
    # Create ATL03 indices (based on segment_id)
    h_ind = np.zeros(len(h_seg))
    h_ind[:] = -1
    
    # Link ATL03 segments to ATL08 segments
    for i in range(0,len(seg_beg)):
        h_ind[np.logical_and(h_seg >= seg_beg[i], h_seg <= seg_end[i])] = i
        
    return h_ind

def computer_vertical_density_n(atl03_df, upsampled_atl08_bin, max_h, min_h, title):
    df_canopy = atl03_df[atl03_df.classification > 1]
    df_select = df_canopy[df_canopy.norm_h < max_h]
    df_select = df_select[df_select.norm_h > min_h]
    # df_select = df_select[df_select.norm_h > df_select.perc0]
    zgroup = df_select.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    dense = zout['norm_h']    
    dense = dense.rename(title)
    dense = dense.reindex(list(range(0,np.max(atl03_df.h_ind)+1)),fill_value=0)
    upsampled_atl08_bin = upsampled_atl08_bin.merge(dense, left_index = True, right_index = True)
    return upsampled_atl08_bin

def computer_vertical_density_n(atl03_df, upsampled_atl08_bin, max_h, min_h, title):
    df_canopy = atl03_df[atl03_df.c > 1]
    df_select = df_canopy[df_canopy.norm_h < max_h]
    df_select = df_select[df_select.norm_h >= min_h]
    # df_select = df_select[df_select.norm_h > df_select.perc0]
    zgroup = df_select.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    dense = zout['norm_h']    
    dense = dense.rename(title)
    dense = dense.reindex(list(range(0,np.max(atl03_df.h_ind)+1)),fill_value=0)
    upsampled_atl08_bin = upsampled_atl08_bin.merge(dense, left_index = True, right_index = True)
    return upsampled_atl08_bin

def computer_vertical_density_p(atl03_df, upsampled_atl08_bin, min_h, max_h, title, h_base = 0):
    df_canopy = atl03_df[atl03_df.c > 1]
    df_canopy = df_canopy[df_canopy.norm_h > h_base]
    zgroup = df_canopy.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    total = zout['norm_h']    
    total = total.rename('total')

    df_select = df_canopy[df_canopy.norm_h < max_h]
    df_select = df_select[df_select.norm_h > min_h]
    # df_select = df_select[df_select.norm_h > df_select.perc0]
    zgroup = df_select.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    dense = zout['norm_h']    
    dense = dense.rename('dense')
    dense = dense.to_frame()
    dense = dense.merge(total, left_index = True, right_index = True)
    p = dense.dense/dense.total    
    p = p.reindex(list(range(0,np.max(atl03_df.h_ind)+1)),fill_value=0)
    p = p.rename(title)
    upsampled_atl08_bin = upsampled_atl08_bin.merge(p, left_index = True, right_index = True)
    return upsampled_atl08_bin

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
    
    return ph_bihr, ph_bcr, ph_rate

def get_binned_rate(atl03_df, upsampled_atl08_bin, atl03filepath, gt):
    ph_bihr, ph_bcr, ph_rate = get_atl03_rate(atl03filepath, gt)
    atl03_df['bckgrd_rate'] = ph_rate
    atl03_df['bckgrd_int_height_reduced'] = ph_bihr
    atl03_df['bckgrd_counts_reduces'] = ph_bcr
    zgroup = atl03_df.groupby('h_ind')
    zout = zgroup.aggregate(np.nanmean)
    f08_rate = np.asarray(zout['bckgrd_rate'])
    f08_bihr = np.asarray(zout['bckgrd_int_height_reduced'])
    f08_bcr = np.asarray(zout['bckgrd_counts_reduces'])
    
    upsampled_atl08_bin['avg_bckgrd_calc_rate'] = f08_rate
    upsampled_atl08_bin['avg_bckgrd_counts_reduces'] = f08_bcr
    upsampled_atl08_bin['avg_bckgrd_int_height_reduced'] = f08_bihr
    
    return upsampled_atl08_bin

def get_photon_rate_ground(atl03_df, upsampled_atl08_bin):
    # Count number of shots
    zgroup = atl03_df.groupby('h_ind')
    zout = zgroup.aggregate(get_len_unique)
    f08_shots = zout['delta_time']
    f08_shots = f08_shots.rename('shots')
    
    # Count number of photons
    atl03_df = atl03_df[atl03_df.c == 1]
    zgroup = atl03_df.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    f08_photons = zout['delta_time']
    f08_photons = f08_photons.rename('photons')
    
    # Calculate photon rate
    f08_shots = f08_shots.to_frame()
    f08_shots = f08_shots.merge(f08_photons, left_index = True, right_index = True)
    rate = f08_shots.photons/f08_shots.shots
    rate = rate.rename('atl03_radiometry_ground')
    upsampled_atl08_bin = upsampled_atl08_bin.merge(rate, left_index = True, right_index = True)
    
    return upsampled_atl08_bin

def get_photon_rate_can(atl03_df, upsampled_atl08_bin):
    # Count number of shots
    zgroup = atl03_df.groupby('h_ind')
    zout = zgroup.aggregate(get_len_unique)
    f08_shots = zout['delta_time']
    f08_shots = f08_shots.rename('shots')
    
    # Count number of photons
    atl03_df = atl03_df[atl03_df.c > 1]
    zgroup = atl03_df.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    f08_photons = zout['delta_time']
    f08_photons = f08_photons.rename('photons')
    
    # Calculate photon rate
    f08_shots = f08_shots.to_frame()
    f08_shots = f08_shots.merge(f08_photons, left_index = True, right_index = True)
    rate = f08_shots.photons/f08_shots.shots
    rate = rate.rename('atl03_radiometry_canopy')
    upsampled_atl08_bin = upsampled_atl08_bin.merge(rate, left_index = True, right_index = True)
    
    return upsampled_atl08_bin

def get_photon_rate_signal(atl03_df, upsampled_atl08_bin):
    # Count number of shots
    zgroup = atl03_df.groupby('h_ind')
    zout = zgroup.aggregate(get_len_unique)
    f08_shots = zout['delta_time']
    f08_shots = f08_shots.rename('shots')
    
    # Count number of photons
    atl03_df = atl03_df[atl03_df.c > 0]
    zgroup = atl03_df.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    f08_photons = zout['delta_time']
    f08_photons = f08_photons.rename('photons')
    
    # Calculate photon rate
    f08_shots = f08_shots.to_frame()
    f08_shots = f08_shots.merge(f08_photons, left_index = True, right_index = True)
    rate = f08_shots.photons/f08_shots.shots
    rate = rate.rename('atl03_radiometry_signal')
    upsampled_atl08_bin = upsampled_atl08_bin.merge(rate, left_index = True, right_index = True)
    
    return upsampled_atl08_bin

def get_n_photons_above_threshold(atl03_df, upsampled_atl08_bin, classes = [2,3], h_base = 0):
    df = atl03_df[atl03_df.c.isin(classes)]
    df = df[df.norm_h > h_base]
    zgroup = df.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    total = zout['norm_h']    
    total = total.rename('n_photons_above_threshold')
    upsampled_atl08_bin = upsampled_atl08_bin.merge(total, left_index = True, right_index = True)
    return upsampled_atl08_bin


def calc_atl03_binning_radiometry(atl03filepath, atl08filepath, gt, res, res_field):
    # Get ATL03 Bin ID

    print('Generate ATL03 Struct')
    atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath)
    
    print('Generate ATL08 Struct')
    atl08 = get_atl08_struct(atl08filepath, gt, atl03)

    
    if res_field in ['delta_time','lat_ph','lon_ph','alongtrack','northing']: 
        at03 = np.array(atl03.df[res_field])
        ind03 = np.int32(np.floor((at03 - np.min(at03))/res))
        atl03.df['h_ind'] = ind03
    elif res_field == 'atl03_seg':
        print('ok')
    elif res_field == 'atl08_seg':
        print('ok')

    print('Match ATL08 to ATL03 by segment')
    upsampled_atl08_bin = create_atl08_bin(atl03, atl08, res_at = res)

    at08 = np.array(upsampled_atl08_bin[res_field])
    ind08 = np.int32(np.floor((at08 - np.min(at03))/res))
    upsampled_atl08_bin['h_ind'] = ind08

    upsampled_atl08_bin.index = ind08
    upsampled_atl08_bin = upsampled_atl08_bin.reindex(list(range(0,np.max(ind03)+1)), method = 'ffill', limit = 3)

    # Normalize the height
    
    atl03.df.rename(columns={'classification':'c'}, inplace=True)
    atl03.df = normalize_heights(atl03.df)
    
    # Go through list of atl03 tings to aggregate
    

    # Average background rates
    upsampled_atl08_bin = get_binned_rate(atl03.df, upsampled_atl08_bin, atl03filepath, gt)

    # Photon rates
    upsampled_atl08_bin = get_photon_rate_ground(atl03.df, upsampled_atl08_bin)
    upsampled_atl08_bin = get_photon_rate_can(atl03.df, upsampled_atl08_bin)
    upsampled_atl08_bin = get_photon_rate_signal(atl03.df, upsampled_atl08_bin)

    # Height bins (%)
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 0, 1, 'h_bin_0')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 1, 2.5, 'h_bin_1')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 2.5, 5, 'h_bin_2')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 5, 7.5, 'h_bin_3')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 7.5, 10, 'h_bin_4')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 10, 12.5, 'h_bin_5')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 12.5, 15, 'h_bin_6')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 15, 17.5, 'h_bin_7')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 17.5, 20, 'h_bin_8')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 20, 22.5, 'h_bin_9')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 22.5, 25, 'h_bin_10')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 25, 27.5, 'h_bin_11')
    upsampled_atl08_bin = computer_vertical_density_p(atl03.df, upsampled_atl08_bin, 27.5, 30, 'h_bin_12')
    upsampled_atl08_bin = get_n_photons_above_threshold(atl03.df, upsampled_atl08_bin, classes = [2,3], h_base = 0)

    
    upsampled_atl08_bin['alongtrack'] = (np.array(upsampled_atl08_bin.index) + np.min(np.array(atl03.df.alongtrack))) * res
    
    # Beam Num
    beamNum = GtToBeamNum(atl03filepath, gt)
    upsampled_atl08_bin['beam'] = beamNum
   
    drop_columns = ['canopy_h_metrics_abs_0',
                    'canopy_h_metrics_abs_1','canopy_h_metrics_abs_2',
                    'canopy_h_metrics_abs_3','canopy_h_metrics_abs_4',
                    'canopy_h_metrics_abs_5','canopy_h_metrics_abs_6',
                    'canopy_h_metrics_abs_7','canopy_h_metrics_abs_8',
                    'canopy_h_metrics_abs_9','canopy_h_metrics_abs_10',
                    'canopy_h_metrics_abs_11','canopy_h_metrics_abs_12',
                    'canopy_h_metrics_abs_13','canopy_h_metrics_abs_14',
                    'canopy_h_metrics_abs_15','canopy_h_metrics_abs_16',
                    'canopy_h_metrics_abs_17']
    upsampled_atl08_bin = upsampled_atl08_bin.drop(columns=drop_columns)
    
    return upsampled_atl08_bin


def rasterize():
    basepath03 = '/laserpewpew/data/release/005/ATL03_r005/Sonoma_nsidc/'
    basepath08 = '/laserpewpew/data/release/005/ATL08_r005/Sonoma_nsidc/'  
    atl03_list = os.listdir(basepath03)
    i = 1
    atl03file = atl03_list[i]
    atl08file = 'ATL08' + atl03_list[i].split('ATL03')[1]
    atl03filepath =  os.path.join(basepath03, atl03file)
    atl08filepath =  os.path.join(basepath08, atl08file)

    gt_list = ['gt1r','gt1l','gt2l','gt2r','gt3l','gt3r']    
    gt = 'gt1r'
    print('Generate ATL03 Struct')
    atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath)
    
    
    


if __name__ == '__main__':
    print('ok')
    if os.name == 'nt':
        basepath03 = 'Z:/data/release/002/ATL03_r002/Finland/'
        basepath08 = 'Z:/data/release/002/ATL08_r002/Finland/'
    else:
        # basepath03 = '/laserpewpew/data/release/004/ATL03_r004/Alberta_nsidc/'
        # basepath08 = '/laserpewpew/data/release/004/ATL08_r004/Alberta_nsidc/'
        # basepath03 = '/laserpewpew/data/release/004/ATL03_r004/Germany_nsidc/'
        # basepath08 = '/laserpewpew/data/release/004/ATL08_r004/Germany_nsidc/'
        basepath03 = '/laserpewpew/data/release/005/ATL03_r005/Sonoma_nsidc/'
        basepath08 = '/laserpewpew/data/release/005/ATL08_r005/Sonoma_nsidc/'
        #basepath03 = '/laserpewpew/data/release/004/ATL03_r004/Congo_nsidc/'
        #basepath08 = '/laserpewpew/data/release/004/ATL08_r004/Alberta_nsidc/'
        # basepath03 = '/laserpewpew/data/release/004/ATL03_r004/Ft_Benning_nsidc/'
        # basepath08 = '/laserpewpew/data/release/004/ATL08_r004/Ft_Benning_nsidc/'
        # basepath03 = '/laserpewpew/data/release/004/ATL03_r004/Ft_Drum_nsidc/'
        # basepath08 = '/laserpewpew/data/release/004/ATL08_r004/Ft_Drum_nsidc/'
        # basepath03 = '/laserpewpew/data/release/004/ATL03_r004/JBLM_nsidc/'
        # basepath08 = '/laserpewpew/data/release/004/ATL08_r004/JBLM_nsidc/'
        # basepath03 = '/laserpewpew/data/release/004/ATL03_r004/LeonardWood_nsidc/'
        # basepath08 = '/laserpewpew/data/release/004/ATL08_r004/LeonardWood_nsidc/'
        # basepath03 = '/laserpewpew/data/release/004/ATL03_r004/Finland_nsidc/'
        # basepath08 = '/laserpewpew/data/release/004/ATL08_r004/Finland_nsidc/'  

    outFilePath = '/LIDAR/server/USERS/eric/for_amy/lai_005/Sonoma_nsidc'
    try:
        os.mkdir(outFilePath)
    except:
        print('Folder already exists')

    
    atl03_list = os.listdir(basepath03)

    # Inputs
    res = 30
    res_field = 'alongtrack'
    # for i in range(0,len(atl03_list)):
    #     try:
    #         print(i)
            # atl03file = atl03_list[i]
            # atl08file = 'ATL08' + atl03_list[i].split('ATL03')[1]
            # atl03filepath =  os.path.join(basepath03, atl03file)
            # atl08filepath =  os.path.join(basepath08, atl08file)

    #         gt_list = ['gt1r','gt1l','gt2l','gt2r','gt3l','gt3r']

    #         upsampled_list = []
    #         for gt in gt_list:
    #             print(gt)
    #             df_bin_out = calc_atl03_binning_radiometry(atl03filepath, atl08filepath, gt, res, res_field)
    
    #             a_dict = {col_name : df_bin_out[col_name].values for col_name \
    #                       in df_bin_out.columns.values}
    #             mat_file = os.path.join(outFilePath, atl08file.split('.')[0] + '_' + gt + '30m.mat')
    #             scipy.io.savemat(mat_file,{'struct':a_dict})
    #             print('Success')
    #     except:
    #         print('fail')