#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 08:26:42 2021

@author: eguenther
"""

import argparse
import h5py
import os
import numpy as np
import pandas as pd
import scipy
import scipy.stats
import random
import time
# import matplotlib.pyplot as plt
# from matplotlib.patches import Rectangle

# from icesatIO import GtToBeamNum, GtToBeamSW
from phoreal.reader import get_atl08_struct
from phoreal.reader import get_atl03_struct

import warnings
warnings.simplefilter("ignore")

# def closest_node(node, nodes):
#     nodes = np.asarray(nodes)
#     deltas = nodes - node
#     dist_2 = np.einsum('ij,ij->i',deltas,deltas)
#     return np.argmin(dist_2)

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

def normalize_heights(df):
    # Statis variables
    ground_res = 0.1
    ground_class = 2
    
    t_ind = np.int32(np.floor((df.alongtrack - np.min(df.alongtrack)) / ground_res))
    df['t_ind'] = t_ind
    df_g = df[df.classification == ground_class]
    zgroup = df_g.groupby('t_ind')
    zout = zgroup.aggregate(np.median)
    # zout = zout.drop(columns=['h_ind', 'c', 'along', 'delta_time',
    #        'bckgrd_int_height_reduced', 'bckgrd_counts_reduces', 'bckgrd_rate','norm_h'])
    zout = zout.reindex(list(range(0,np.max(t_ind) + 1)))
    zout = zout.interpolate(method='linear', axis=0).ffill().bfill()
    ground = zout.h_ph[df.t_ind]
    norm_height = np.array(df.h_ph) - np.array(ground)
    
    # Re assign points below the ground to the ground class
    df.classification[norm_height < 0] = 0 
    
    # 
    norm_height[norm_height < 0] = 0
    
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

# def computer_rh(atl03_df, upsampled_atl08_bin, percentile, title):
#     df_canopy = atl03_df[atl03_df.c > 1]
#     # df_select = df_select[df_select.norm_h > df_select.perc0]
#     zgroup = df_canopy.groupby('h_ind')
#     if percentile == 5:
#         zout = zgroup.aggregate(get_max5)
#     elif percentile == 10:
#         zout = zgroup.aggregate(get_max10)
#     elif percentile == 15:
#         zout = zgroup.aggregate(get_max15)
#     elif percentile == 20:
#         zout = zgroup.aggregate(get_max20)
#     elif percentile == 25:
#         zout = zgroup.aggregate(get_max25)
#     elif percentile == 30:
#         zout = zgroup.aggregate(get_max30)
#     elif percentile == 40:
#         zout = zgroup.aggregate(get_max40)
#     elif percentile == 50:
#         zout = zgroup.aggregate(get_max50)
#     elif percentile == 60:
#         zout = zgroup.aggregate(get_max60)
#     elif percentile == 70:
#         zout = zgroup.aggregate(get_max70)
#     elif percentile == 75:
#         zout = zgroup.aggregate(get_max75)
#     elif percentile == 80:
#         zout = zgroup.aggregate(get_max80)
#     elif percentile == 90:
#         zout = zgroup.aggregate(get_max90)
#     elif percentile == 95:
#         zout = zgroup.aggregate(get_max95)
#     elif percentile == 98:
#         zout = zgroup.aggregate(get_max98)
#     elif percentile == 100:
#         zout = zgroup.aggregate(get_max100)
#     dense = zout['norm_h']    
#     dense = dense.rename(title)
#     dense = dense.reindex(list(range(0,np.max(atl03_df.h_ind)+1)),fill_value=0)
#     upsampled_atl08_bin = upsampled_atl08_bin.merge(dense, left_index = True, right_index = True)
#     return upsampled_atl08_bin

def computer_vertical_density_n(atl03_df, upsampled_atl08_bin, max_h, min_h, title, veg_class = [2,3]):
    # df_canopy = atl03_df[atl03_df.classification > 1]
    df_canopy = atl03_df[atl03_df.classification.isin(veg_class)]
    df_select = df_canopy[df_canopy.norm_h < max_h]
    df_select = df_select[df_select.norm_h >= min_h]
    zgroup = df_select.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    dense = zout['norm_h']    
    dense = dense.rename(title)
    dense = dense.reindex(list(range(0,np.max(atl03_df.h_ind)+1)),fill_value=0)
    upsampled_atl08_bin = upsampled_atl08_bin.merge(dense, left_index = True, right_index = True)
    return upsampled_atl08_bin

def computer_vertical_density_p(atl03_df, upsampled_atl08_bin, min_h, max_h, title, veg_class = [2,3], h_base = 0):
    # df_canopy = atl03_df[atl03_df.classification > 1]
    df_canopy = atl03_df[atl03_df.classification.isin(veg_class)]
    df_canopy = df_canopy[df_canopy.norm_h > h_base]
    zgroup = df_canopy.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    total = zout['norm_h']    
    total = total.rename('total')

    df_select = df_canopy[df_canopy.norm_h < max_h]
    df_select = df_select[df_select.norm_h > min_h]
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

def get_atl03_height(atl03filepath, atl08filepath, gt):
    f = h5py.File(atl03filepath, 'r')
    segment_id = np.array(f[gt + '/geolocation/segment_id'])
    ph_index_beg = np.array(f[gt + '/geolocation/ph_index_beg'])
    # delta_time = np.array(f[gt + '/geolocation/delta_time'])
    f = h5py.File(atl08filepath, 'r')
    classed_pc_indx = np.array(f[gt + '/signal_photons/classed_pc_indx'])
    #classed_pc_flag = np.array(f[gt + '/signal_photons/classed_pc_flag'])
    classed_pc_flag = np.array(f[gt + '/signal_photons/ph_h'])
    seg08_id = np.array(f[gt + '/signal_photons/ph_segment_id'])
    f.close()

    indsNotZero = ph_index_beg != 0
    ph_index_beg = ph_index_beg[indsNotZero];
    segment_id = segment_id[indsNotZero];
    
    # Find ATL08 segments that have ATL03 segments
    atl03SegsIn08TF, atl03SegsIn08Inds = ismember(seg08_id,segment_id)
    
    # Get ATL08 classed indices and values
    atl08classed_inds = classed_pc_indx[atl03SegsIn08TF]
    atl08classed_vals = classed_pc_flag[atl03SegsIn08TF]

    # Determine new mapping into ATL03 data
    atl03_ph_beg_inds = atl03SegsIn08Inds;
    atl03_ph_beg_val = ph_index_beg[atl03_ph_beg_inds];
    newMapping = atl08classed_inds + atl03_ph_beg_val - 2;
    
    # Get max size of output array
    sizeOutput = newMapping[-1]
    
    # Pre-populate all photon classed array with zeroes
    #allph_classed = (np.zeros(sizeOutput + 1).astype(int)) - 1
    adjust_heights = (np.zeros(sizeOutput + 1))
    
    # Populate all photon classed array from ATL08 classifications
    adjust_heights[newMapping] = atl08classed_vals
    
    # Return all photon classed array
    return adjust_heights


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
    atl03_df['avg_bckgrd_calc_rate'] = ph_rate
    atl03_df['avg_bckgrd_counts_reduces'] = ph_bihr
    atl03_df['avg_bckgrd_int_height_reduced'] = ph_bcr
    zgroup = atl03_df.groupby('h_ind')
    zout = zgroup.aggregate(np.nanmean)
    # f08_rate = np.asarray(zout['bckgrd_rate'])
    # f08_bihr = np.asarray(zout['bckgrd_int_height_reduced'])
    # f08_bcr = np.asarray(zout['bckgrd_counts_reduces'])
    f08_rate = zout['avg_bckgrd_calc_rate']
    f08_bihr = zout['avg_bckgrd_counts_reduces']
    f08_bcr = zout['avg_bckgrd_int_height_reduced']
    f08_rate.reindex(range(upsampled_atl08_bin.index.min(), upsampled_atl08_bin.index.max()+1),fill_value=0)
    f08_bihr.reindex(range(upsampled_atl08_bin.index.min(), upsampled_atl08_bin.index.max()+1),fill_value=0)
    f08_bcr.reindex(range(upsampled_atl08_bin.index.min(), upsampled_atl08_bin.index.max()+1),fill_value=0)
    # upsampled_atl08_bin['avg_bckgrd_calc_rate'] = f08_rate
    # upsampled_atl08_bin['avg_bckgrd_counts_reduces'] = f08_bcr
    # upsampled_atl08_bin['avg_bckgrd_int_height_reduced'] = f08_bihr
    upsampled_atl08_bin = upsampled_atl08_bin.merge(f08_rate, 
                                        left_index = True, right_index = True)
    upsampled_atl08_bin = upsampled_atl08_bin.merge(f08_bihr, 
                                        left_index = True, right_index = True)
    upsampled_atl08_bin = upsampled_atl08_bin.merge(f08_bcr, 
                                        left_index = True, right_index = True)
    return upsampled_atl08_bin

def get_n_photons_above_threshold(atl03_df, upsampled_atl08_bin, classes = [2,3], h_base = 0):
    df = atl03_df[atl03_df.classification.isin(classes)]
    df = df[df.norm_h > h_base]
    zgroup = df.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    total = zout['norm_h']    
    total = total.rename('n_photons_above_threshold')
    # upsampled_atl08_bin = upsampled_atl08_bin.merge(total, left_index = True, right_index = True)
    df_bin_out = pd.merge(upsampled_atl08_bin, total, how='left',on='h_ind')
    return df_bin_out


def interpolate_domain(atl08_at, atl08_domain, key_df_at, kind_type):
    intep_func = scipy.interpolate.interp1d(atl08_at, atl08_domain, kind=kind_type, 
                              fill_value='extrapolate')
    interp_domain = intep_func(key_df_at)
    return interp_domain

def get_max(series):
    try:
        min0 = np.max(series)
    except:
        min0 = np.nan
    return min0

def get_min(series):
    try:
        min0 = np.min(series)
    except:
        min0 = np.nan
    return min0

def get_len(series):
    try:
        length = len(series)
    except:
        length = np.nan
    return length

def get_len_unique(series):
    try:
        length = len(np.unique(series))
    except:
        length = np.nan
    return length

def get_range(series):
    try:
        r = np.ptp(series)
    except:
        r = np.nan
    return r

def get_std(series):
    try:
        s = np.std(series)
    except:
        s = np.nan
    return s

def get_mode(series):
    try:
        mode = scipy.stats.mode(np.array(series))[0][0]
    except:
        mode = np.nan
    return mode
    
def get_mean(series):
    try:
        m = np.mean(series)
    except:
        m = np.nan
    return m

def get_median(series):
    try:
        m = np.median(series)
    except:
        m = np.nan
    return m

def get_quad_mean(series):
    try:
        quad = np.sqrt(np.mean(series**2))
    except:
        quad = np.nan
    return quad

def get_skew(series):
    try:
        skew = scipy.stats.skew(series)
    except:
        skew = np.nan
    return skew

def get_max25(series):
    try:
        max25 = np.percentile(series, 25)
    except:
        max25 = np.nan
    return max25

def _percentile_post(arr, percent):
    pos = (len(arr) - 1) * (percent/100)

    if pos.is_integer():
        out = arr[int(pos)]
    else:
        out = ((arr[int(np.ceil(pos))] - arr[int(np.floor(pos))]) *\
               (pos % 1)) + arr[int(np.floor(pos))]
    return out

def percentile_rh(arr):
    arr = np.array(arr)
    arr.sort()
    out_list = []    
    try:    
        out_list.append(_percentile_post(arr, 10))
        out_list.append(_percentile_post(arr, 20))
        out_list.append(_percentile_post(arr, 25))
        out_list.append(_percentile_post(arr, 30))
        out_list.append(_percentile_post(arr, 40))
        out_list.append(_percentile_post(arr, 50))
        out_list.append(_percentile_post(arr, 60))
        out_list.append(_percentile_post(arr, 70))
        out_list.append(_percentile_post(arr, 75))
        out_list.append(_percentile_post(arr, 80))
        out_list.append(_percentile_post(arr, 90))
        out_list.append(_percentile_post(arr, 98))
        out_list.append(_percentile_post(arr, 100))
    except:
        out_list = [np.nan, np.nan,np.nan,np.nan,np.nan,np.nan,
                                np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,
                                np.nan]

    return out_list


def calculate_seg_meteric(df_in, df_out, classification, operation, field, 
               outfield, key_field = 'bin_id', classfield = 'classification'):
    df_filter = df_in[df_in[classfield].isin(classification)]
    df_filter.drop(df_filter.columns.difference([key_field,field]),
                       1,inplace=True)
    zgroup = df_filter.groupby(key_field)
    zout = zgroup.aggregate(operation)
    zout[key_field] = zout.index
    zout = zout.reset_index(drop = True)
    # zout['segment_id_beg'] = zout['seg_id']
    zout[outfield] = zout[field]
    zout = zout.filter([outfield,key_field])
    df_out = df_out.merge(zout, on=key_field,how='left')  
    return df_out

def calculate_seg_percentile(df_in, df_out, classification, operation, field, 
               outfield, key_field = 'bin_id', classfield = 'classification'):
    q_list = [10,20,25,30,40,50,60,70,75,80,90,98,100]
    df_filter = df_in[df_in[classfield].isin(classification)]
    df_filter.drop(df_filter.columns.difference([key_field,field]), 
                   1,inplace=True)
    zgroup = df_filter.groupby(key_field)
    zout = zgroup.aggregate(percentile_rh)
    zout[key_field] = zout.index
    zout = zout.reset_index(drop = True)
    # field = field
    # outfield = outfield
    # zout['segment_id_beg'] = zout['seg_id']
    zout[outfield] = zout[field]
    zout = zout.filter([outfield,key_field])
    # print(zout.outfield[0][0])
    # print(time.time() - d3)
    for i in range(0,len(q_list)):
        arr = []
        for j in range(0,len(zout[outfield])):
            arr.append(zout[outfield][j][i])
        title = outfield + str(q_list[i])
        # print(title)
        zout[title] = np.array(arr)
    zout.drop(columns=[outfield], inplace=True)
    df_out = df_out.merge(zout, on=key_field,how='left')  
    return df_out

def slope_df(df, df_bin, ground_class):
    # Static variable
    field_name = 'terrain_slope'
    
    df = df[df.classification == ground_class]
    df = df.sort_values(by=['alongtrack'])
    df = df.reset_index()
    pos_df = df.groupby('h_ind')
    # slope = []
    indice_list = []
    slope_list = []
    # bin_df['terrain_slope'] = np.nan
    for ind in np.unique(df.h_ind):
        df_ind = df[df.h_ind == ind]
        if len(df_ind) > 1:
            # rise = np.array(df_ind.h_ph)[-1] - np.array(df_ind.h_ph)[0]
            # run = np.array(df_ind.alongtrack)[-1] - np.array(df_ind.alongtrack)[0]
            # slope.append(rise/run)
            slope_calc, _ = np.polyfit(np.array(df_ind.alongtrack),np.array(df_ind.h_ph),1)
            slope_list.append(slope_calc)
            indice_list.append(ind)
        else:
            # slope.append(np.nan)
            slope_list.append(np.nan)
            indice_list.append(ind)
            
    data = {'h_ind':indice_list,field_name:slope_list}
    df_slope = pd.DataFrame(data)
    df_bin_out = pd.merge(df_bin, df_slope, how='left',on='h_ind')
    # df_slope = df_slope.set_index('h_ind')
            
    return df_bin_out



def sub_bin_canopy_metrics(df, res_at = 2):
    at03 = np.array(df['alongtrack'])
    ind03 = np.int32(np.floor((at03 - np.min(at03))/res_at))
    df['c_ind'] = ind03
    df_c = df[df.classification > 1]
    # df_c = df_c.reindex(list(range(ind03.min(),ind03.max()+1)),fill_value=0)
    zgroup = df_c.groupby('c_ind')
    sub_bin = zgroup.aggregate(np.max)
    sub_bin = sub_bin.reindex(list(range(ind03.min(),ind03.max()+1)),fill_value=0)
    sub_bin['alongtrack'] = (sub_bin.index * res_at) + (res_at / 2)

    domain_interp = interpolate_domain(df.alongtrack, 
                                       df.h_ind, sub_bin.alongtrack + np.min(at03), 'nearest')
    
    norm_h = np.array(sub_bin.norm_h)
    a_h = np.array([x - norm_h[i - 1] for i, x in enumerate(norm_h)][1:])
    a_h = np.append(a_h, 0)
    # test_out = np.pad(test_out, (0,1), 'constant')
    surface_area = np.sqrt((a_h**2) + (res_at**2))
    sub_bin = sub_bin.drop(columns=['delta_time', 'dist_ph_across', 'dist_ph_along', 'h_ph', 'lat_ph',
       'lon_ph', 'pce_mframe_cnt', 'ph_id_channel', 'ph_id_count',
       'ph_id_pulse', 'quality_ph', 'signal_conf_ph', 'classification','seg_id', 'time', 'easting', 'northing', 'crosstrack',
      'alongtrack','ph_bihr','ph_bcr','ph_rate'])
    sub_bin['surface_len'] = surface_area
    sub_bin['ground_len'] = res_at
    sub_bin['h_ind'] = domain_interp.astype(int)

    # sub_bin = sub_bin[sub_bin['h_ind'] != 0]
    return sub_bin

def sub_bin_canopy_metrics_truth(df, res_at = 2):
    at = np.array(df['alongtrack'])
    ind = np.int32(np.floor((at - np.min(at))/res_at))
    df['c_ind'] = ind
    df_c = df[df.classification == 4]
    # df_c = df_c.reindex(list(range(ind03.min(),ind03.max()+1)),fill_value=0)
    zgroup = df_c.groupby('c_ind')
    sub_bin = zgroup.aggregate(np.max)
    sub_bin = sub_bin.reindex(list(range(ind.min(),ind.max()+1)),fill_value=0)
    sub_bin['alongtrack'] = (sub_bin.index * res_at) + (res_at / 2)

    domain_interp = interpolate_domain(df.alongtrack, 
                                       df.h_ind, sub_bin.alongtrack + np.min(at), 'nearest')
    
    norm_h = np.array(sub_bin.norm_h)
    a_h = np.array([x - norm_h[i - 1] for i, x in enumerate(norm_h)][1:])
    a_h = np.append(a_h, 0)
    # test_out = np.pad(test_out, (0,1), 'constant')
    surface_area = np.sqrt((a_h**2) + (res_at**2))
    if 'date' in sub_bin.columns:
        sub_bin = sub_bin.drop(columns='date')
    sub_bin = sub_bin.drop(columns=['h_ph', 'lat',
       'lon', 'x', 'y','classification','intensity', 'easting', 'northing', 'crosstrack',
      'alongtrack'])
    sub_bin['surface_len'] = surface_area
    sub_bin['ground_len'] = res_at
    # sub_bin['h_ind'] = domain_interp.astype(int)

    # sub_bin = sub_bin[sub_bin['h_ind'] != 0]
    return sub_bin

def calc_rumple_index(sub_bin):
    zgroup = sub_bin.groupby('h_ind')
    veg_df = zgroup.aggregate(np.sum)
    veg_df['rumple_index'] = veg_df.surface_len / veg_df.ground_len
    veg_df['veg_area'] = veg_df.norm_h
    veg_df = veg_df.drop(columns=['norm_h'])    
    return veg_df

def rebin_atl08(atl03, atl08, gt, res, res_field):
    # Get ATL03 Bin ID

    t1 = time.time()
    percentile_intervals = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 
                            70, 75, 80, 85, 90, 95, 98]
    
    # print('Generate ATL03 Struct')
    # atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath)
    # # atl03.df = atl03.df[atl03.df.lat_ph > 57]
    # # atl03.df = atl03.df[atl03.df.lat_ph < 65]
    # #atl03.df = atl03.df[atl03.df.lon_ph > 88]
    # #atl03.df = atl03.df[atl03.df.lon_ph < 106]
    
    # print('Generate ATL08 Struct')
    # atl08 = get_atl08_struct(atl08filepath, gt, atl03)
    # atl08.df = atl08.df[atl08.df.latitude > 57]
    # atl08.df = atl08.df[atl08.df.latitude < 65]
    #atl08.df = atl08.df[atl08.df.longitude > 88]
    #atl08.df = atl08.df[atl08.df.longitude < 106]

    
    if res_field in ['delta_time','lat_ph','lon_ph','alongtrack','northing']: 
        at03 = np.array(atl03.df[res_field])
        ind03 = np.int32(np.floor((at03 - np.min(at03))/res))
        atl03.df['h_ind'] = ind03

    #Place holder for ATL08
        at08 = np.array(atl08.df[res_field])
        ind08 = np.int32(np.floor((at08 - np.min(at03))/res))
        atl08.df['h_ind'] = ind08
        
    elif res_field == 'atl03_seg':
        print('ok')
    elif res_field == 'atl08_seg':
        print('ok')
        
    
    linear_field_list = ['alongtrack', 'crosstrack', 'time', 'easting', 
                           'northing', 'asr', 'atlas_pa', 'beam_azimuth', 
                           'beam_coelev', 
                         'delta_time', 'dem_h', 'h_canopy_uncertainty', 
                         'latitude', 'longitude', 'sigma_across', 
                         'sigma_along', 'sigma_atlas_land', 'sigma_h', 
                         'sigma_topo', 'snr', 'solar_azimuth', 
                         'solar_elevation', 'h_te_best_fit','h_te_interp']
    
    bin_df = pd.DataFrame(np.arange(np.min(atl03.df.h_ind),np.max(atl03.df.h_ind)),columns=['h_ind'])
    
    for field in linear_field_list:
        domain_interp = interpolate_domain(atl08.df.h_ind, 
                                           atl08.df[field], bin_df.h_ind, 'linear')
        bin_df = pd.concat([bin_df,pd.DataFrame(domain_interp,
                                                columns=[field])],axis=1)

    nearest_field_list = ['brightness_flag','segment_cover','subset_can_flag',
                          'cloud_flag_atm','cloud_fold_flag','dem_flag',
                          'dem_removal_flag','layer_flag','msw_flag',
                          'night_flag','ph_removal_flag','psf_flag','rgt',
                          'sat_flag','segment_landcover','segment_snowcover',
                          'segment_watermask','surf_type','subset_te_flag',
                          'terrain_flg','urban_flag', 'surf_type_0',
                          'surf_type_1','surf_type_2','surf_type_3',
                          'surf_type_4','h_te_median','h_canopy', 
                          'ph_removal_flag']

    for field in nearest_field_list:
        try:
            domain_interp = interpolate_domain(atl08.df.h_ind, 
                                               atl08.df[field], bin_df.h_ind, 'nearest')
            bin_df = pd.concat([bin_df,pd.DataFrame(domain_interp,
                                                    columns=[field])],axis=1)
        except:
            print(field)
    
    # Apply nan terrain and canopy filter from ATL08'
    ### put new code here***
    bin_df['ground_filter'] = [np.array(bin_df.h_te_median) < 3.402823e+37][0]
    bin_df['canopy_filter'] = [np.array(bin_df.h_canopy) < 3.402823e+37][0]
    bin_df = bin_df.drop(columns=['h_te_median'])
    bin_df = bin_df.drop(columns=['h_canopy'])
    
    print('The first part...')
    print(time.time() - t1)
    t1 = time.time()
    # Compute standard metrics
    # canopy_rh
    bin_df = calculate_seg_percentile(atl03.df, bin_df, [1, 2,3], percentile_rh, 'norm_h', 
                   'canopy_rh_ground_', key_field = 'h_ind', classfield = 'classification')

    bin_df = calculate_seg_percentile(atl03.df, bin_df, [2,3], percentile_rh, 'norm_h', 
                   'canopy_rh_', key_field = 'h_ind', classfield = 'classification')
    # canopy_rh_abs
    bin_df = calculate_seg_percentile(atl03.df, bin_df, [2,3], percentile_rh, 'h_ph', 
                   'canopy_rh_abs_', key_field = 'h_ind', classfield = 'classification')
    #canopy_openness
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_std, 'norm_h', 
                   'canopy_openness', key_field = 'h_ind', classfield = 'classification')
    #centroid_height
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1,2,3], get_mean, 'h_ph', 
                   'centroid_height', key_field = 'h_ind', classfield = 'classification')
    #h_canopy_quad
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_quad_mean, 'norm_h', 
                   'h_canopy_quad', key_field = 'h_ind', classfield = 'classification')
    #h_max_canopy
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_max, 'norm_h', 
                   'h_max_canopy', key_field = 'h_ind', classfield = 'classification')
    #h_max_canopy_abs
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_max, 'h_ph', 
                   'h_max_canopy_abs', key_field = 'h_ind', classfield = 'classification')
    #h_mean_canopy
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_mean, 'norm_h', 
                   'h_mean_canopy', key_field = 'h_ind', classfield = 'classification')
    #h_mean_canopy_abs
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_mean, 'h_ph', 
                   'h_canopy_quad', key_field = 'h_ind', classfield = 'classification')
    #h_median_canopy
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_median, 'norm_h', 
                   'h_median_canopy', key_field = 'h_ind', classfield = 'classification')
    #h_median_canopy_abs 
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_median, 'h_ph', 
                   'h_median_canopy_abs', key_field = 'h_ind', classfield = 'classification')
    #h_min_canopy
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_min, 'norm_h', 
                   'h_min_canopy', key_field = 'h_ind', classfield = 'classification')
    #h_min_canopy_abs 
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2,3], get_min, 'h_ph', 
                   'h_min_canopy_abs', key_field = 'h_ind', classfield = 'classification')
    #n_ca_photons
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [2], get_len, 'h_ph', 
                   'n_ca_photons', key_field = 'h_ind', classfield = 'classification')
    #n_toc_photons
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [3], get_len, 'h_ph', 
                   'n_toc_photons', key_field = 'h_ind', classfield = 'classification')
    #toc_roughness
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [3], get_std, 'norm_h', 
                   'toc_roughness', key_field = 'h_ind', classfield = 'classification')
    #h_te_skew
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1], get_skew, 'alongtrack', 
                   'h_canopy_quad', key_field = 'h_ind', classfield = 'classification')
    #h_te_std
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1], get_std, 'h_ph', 
                   'h_te_std', key_field = 'h_ind', classfield = 'classification')
    #h_te_max
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1], get_max, 'h_ph', 
                   'h_te_max', key_field = 'h_ind', classfield = 'classification')
    #h_te_min
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1], get_min, 'h_ph', 
                   'h_te_min', key_field = 'h_ind', classfield = 'classification')
    #h_te_mean
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1], get_mean, 'h_ph', 
                   'h_canopy_quad', key_field = 'h_ind', classfield = 'classification')
    #h_te_median
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1], get_median, 'h_ph', 
                   'h_te_median', key_field = 'h_ind', classfield = 'classification')
    #h_te_mode
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1], get_mode, 'h_ph', 
                   'h_te_mode', key_field = 'h_ind', classfield = 'classification')
    #n_te_photons
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1], get_len, 'h_ph', 
                   'n_te_photons', key_field = 'h_ind', classfield = 'classification')
    #n_unique_shots
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [-1,0,1,2,3], get_len_unique, 'delta_time', 
                   'n_unique_shots', key_field = 'h_ind', classfield = 'classification')
    #h_te_rh25
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [1], get_max25, 'h_ph', 
                   'h_te_rh25', key_field = 'h_ind', classfield = 'classification')
    #n_seg_ph
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [-1,0,1,2,3], get_len, 'h_ph', 
                   'n_seg_ph', key_field = 'h_ind', classfield = 'classification')
    #avg_bckgrd_calc_rate
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [-1,0,1,2,3], get_mean, 'ph_bihr', 
                   'avg_bckgrd_calc_rate', key_field = 'h_ind', classfield = 'classification')
    #avg_bckgrd_counts_reduces
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [-1,0,1,2,3], get_mean, 'ph_bcr', 
                   'avg_bckgrd_counts_reduces', key_field = 'h_ind', classfield = 'classification')
    #avg_bckgrd_int_height_reduced
    bin_df = calculate_seg_meteric(atl03.df, bin_df, [-1,0,1,2,3], get_mean, 'ph_rate', 
                   'avg_bckgrd_int_height_reduced', key_field = 'h_ind', classfield = 'classification')
    print('The big compute...')
    print(time.time() - t1)
    t1 = time.time()
    
    
    # Compute photon_rate_can
    bin_df['photon_rate_can'] = (bin_df.n_ca_photons + bin_df.n_toc_photons) / bin_df.n_unique_shots
    
    # Compute photon_rate_te
    bin_df['photon_rate_te'] = bin_df.n_te_photons / bin_df.n_unique_shots
    
    # # Unstack canopy_rh and rename (possible to drop canopy_rh_abs)
    # for i in range(0,len(percentile_intervals)):
    #     arr = []
    #     for j in range(0,len(bin_df)):
    #         if type(bin_df.iloc[j].canopy_rh) == np.ndarray:
    #             arr.append(bin_df.iloc[j].canopy_rh[i])
    #         else:
    #             arr.append(0)
            
    #     arr = np.array(arr)
    #     title = 'canopy_rh_' + str(percentile_intervals[i])
    #     bin_df[title] = arr

    # # Unstack canopy_rh and rename (possible to drop canopy_rh_abs)
    # for i in range(0,len(percentile_intervals)):
    #     arr = []
    #     for j in range(0,len(bin_df)):
    #         if type(bin_df.iloc[j].canopy_rh_abs) == np.ndarray:
    #             arr.append(bin_df.iloc[j].canopy_rh_abs[i])
    #         else:
    #             arr.append(0)
            
    #     arr = np.array(arr)
    #     title = 'canopy_rh_abs_' + str(percentile_intervals[i])
    #     bin_df[title] = arr
        
    print('The unpacking...')
    print(time.time() - t1)
    t1 = time.time()

    bin_df['canopy_h_ground'] = bin_df['canopy_rh_98'] 
    bin_df['canopy_h'] = bin_df['canopy_rh_98'] 
    bin_df['h_canopy_abs'] = bin_df['canopy_rh_abs_98'] #        
    
    # bin_df = bin_df.drop(columns=['canopy_rh'])
    # bin_df = bin_df.drop(columns=['canopy_rh_abs'])
    
    # Compute h_dif_canopy
    bin_df['photon_rate_te'] = bin_df.n_te_photons / bin_df.n_unique_shots
        
    # Compute slope
    # calc_df['terrain_slope'] = calc_df.h_te_rise / calc_df.h_te_run
    # calc_df = calc_df.drop(columns=['h_te_rise','h_te_run'])

    # Compute h_dif_can
    bin_df['h_dif_canopy'] = bin_df.canopy_h - bin_df.canopy_rh_50
    
    # Compute traditional canopy cover
    bin_df['traditional_canopy_cover'] = (bin_df['n_ca_photons'] + \
            bin_df['n_toc_photons']) / (bin_df['n_ca_photons'] + \
            bin_df['n_toc_photons'] + bin_df['n_te_photons'])
    
    print('The trad canopy cover...')
    print(time.time() - t1)
    t1 = time.time()
    # Compute Rumple Index (canopy surface area/ground surface area)
    sub_bin = sub_bin_canopy_metrics(atl03.df, res_at = 2)
    veg_bin = calc_rumple_index(sub_bin)
    bin_df = bin_df.join(veg_bin)
            
    print('The sub_bin_canopy_metrics...')
    print(time.time() - t1)
    t1 = time.time()
    # Compute Area-Under-Canopy/Open Sky Ratio
    bin_df['canopy_ratio'] = bin_df.veg_area / (bin_df.h_max_canopy * 
                                                  bin_df.ground_len)
    # Compute slope
    bin_df = slope_df(atl03.df, bin_df, 1)
    print('Slope calc...')
    print(time.time() - t1)
    t1 = time.time()
    
    # # Merge Computed columns to bin_df
    # bin_df = bin_df.join(calc_df)
    
    # Compute h_dif_ref
    bin_df['h_dif_ref'] = bin_df.h_te_median - bin_df.dem_h
    
    # Height bins (%)
    h_bin_list = [[0,1],[1,2.5],[2.5,5],[5,7.5],[7.5,10],[10,12.5],[12.5,15],
                  [15,17.5],[17.5,20],[20,22.5],[22.5,25],[25,27.5],[27.5,30],
                  [30,200]]
    h_base_val = 0
    h_bin_prefix = 'h_bin_'
    for i in range(0,len(h_bin_list)):
        title = h_bin_prefix + str(i)
        bin_df = computer_vertical_density_p(atl03.df, bin_df, 
                        h_bin_list[i][0], h_bin_list[i][1], title, 
                        veg_class = [2,3], h_base = h_base_val)

    print('Vertical Height Bins...')
    print(time.time() - t1)
    t1 = time.time()
    # Compute photons above set threshold
    bin_df = get_n_photons_above_threshold(atl03.df, bin_df, classes = [2,3], 
                                           h_base = h_base_val)

    print('Photons above threshold...')
    print(time.time() - t1)
    t1 = time.time()

    bin_df['beamNum'] = atl03.beamNum
    bin_df['beamStrength'] = atl03.beamStrength
    bin_df['year'] = atl03.year
    bin_df['month'] = atl03.month
    bin_df['day'] = atl03.day
    bin_df['epsg'] = atl03.epsg

    print('The last part...')
    print(time.time() - t1)
    t1 = time.time()

    # Return bin_df
    return bin_df

def rebin_truth(atl03, truth_swath, res, res_field):
    ground_class = 2
    veg_class = 4
    other_veg_classes = [3, 5]
    
    # Get ATL03 Bin ID

    percentile_intervals = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 
                            70, 75, 80, 85, 90, 95, 98]
    
    for other_veg_class in other_veg_classes:
        truth_swath.classification\
            [truth_swath.classification == other_veg_class] = veg_class

    
    if res_field in ['delta_time','lat_ph','lon_ph','alongtrack','northing']: 
        at03 = np.array(atl03.df[res_field])

        at_truth = np.array(truth_swath[res_field])
        indtruth = np.int32(np.floor((at_truth - np.min(at03))/res))
        truth_swath['h_ind'] = indtruth
        
    elif res_field == 'atl03_seg':
        print('ok')
    elif res_field == 'atl08_seg':
        print('ok')
        
    bin_df = pd.DataFrame(np.arange(np.min(truth_swath.h_ind),np.max(truth_swath.h_ind)),columns=['h_ind'])

    
    linear_field_list = ['alongtrack', 'crosstrack', 'easting', 
                           'northing', 'lon', 'lat']
    
    
    for field in linear_field_list:
        domain_interp = interpolate_domain(truth_swath.h_ind, 
                                           truth_swath[field], bin_df.h_ind, 'linear')
        bin_df = pd.concat([bin_df,pd.DataFrame(domain_interp,
                                                columns=[field])],axis=1)


    # Rename h_ph
    truth_swath = truth_swath.rename(columns={'z':'h_ph'})

    # Create normalized ground
    truth_swath = normalize_heights(truth_swath)

    # Calculate fields      
    # canopy_rh
    bin_df = calculate_seg_percentile(truth_swath, bin_df, [4], percentile_rh, 'norm_h', 
                   'canopy_rh_', key_field = 'h_ind', classfield = 'classification')
    # canopy_rh_abs
    bin_df = calculate_seg_percentile(truth_swath, bin_df, [4], percentile_rh, 'h_ph', 
                   'canopy_rh_abs_', key_field = 'h_ind', classfield = 'classification')
    #canopy_openness
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_std, 'norm_h', 
                   'canopy_openness', key_field = 'h_ind', classfield = 'classification')
    #centroid_height
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2,4], get_mean, 'h_ph', 
                   'centroid_height', key_field = 'h_ind', classfield = 'classification')
    #h_canopy_quad
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_quad_mean, 'norm_h', 
                   'h_canopy_quad', key_field = 'h_ind', classfield = 'classification')
    #h_max_canopy
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_max, 'norm_h', 
                   'h_max_canopy', key_field = 'h_ind', classfield = 'classification')
    #h_max_canopy_abs
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_max, 'h_ph', 
                   'h_max_canopy_abs', key_field = 'h_ind', classfield = 'classification')
    #h_mean_canopy
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_mean, 'norm_h', 
                   'h_mean_canopy', key_field = 'h_ind', classfield = 'classification')
    #h_mean_canopy_abs
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_mean, 'h_ph', 
                   'h_canopy_quad', key_field = 'h_ind', classfield = 'classification')
    #h_median_canopy
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_median, 'norm_h', 
                   'h_median_canopy', key_field = 'h_ind', classfield = 'classification')
    #h_median_canopy_abs 
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_median, 'h_ph', 
                   'h_median_canopy_abs', key_field = 'h_ind', classfield = 'classification')
    #h_min_canopy
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_min, 'norm_h', 
                   'h_min_canopy', key_field = 'h_ind', classfield = 'classification')
    #h_min_canopy_abs 
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_min, 'h_ph', 
                   'h_min_canopy_abs', key_field = 'h_ind', classfield = 'classification')
    #n_ca_photons
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [4], get_len, 'h_ph', 
                   'n_ca_photons', key_field = 'h_ind', classfield = 'classification')
    #h_te_skew
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2], get_skew, 'alongtrack', 
                   'h_canopy_quad', key_field = 'h_ind', classfield = 'classification')
    #h_te_std
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2], get_std, 'h_ph', 
                   'h_te_std', key_field = 'h_ind', classfield = 'classification')
    #h_te_max
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2], get_max, 'h_ph', 
                   'h_te_max', key_field = 'h_ind', classfield = 'classification')
    #h_te_min
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2], get_min, 'h_ph', 
                   'h_te_min', key_field = 'h_ind', classfield = 'classification')
    #h_te_mean
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2], get_mean, 'h_ph', 
                   'h_canopy_quad', key_field = 'h_ind', classfield = 'classification')
    #h_te_median
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2], get_median, 'h_ph', 
                   'h_te_median', key_field = 'h_ind', classfield = 'classification')
    #h_te_mode
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2], get_mode, 'h_ph', 
                   'h_te_mode', key_field = 'h_ind', classfield = 'classification')
    #n_te_photons
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2], get_len, 'h_ph', 
                   'n_te_photons', key_field = 'h_ind', classfield = 'classification')
    #h_te_rh25
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2], get_max25, 'h_ph', 
                   'h_te_rh25', key_field = 'h_ind', classfield = 'classification')
    #n_seg_ph
    bin_df = calculate_seg_meteric(truth_swath, bin_df, [2,4], get_len, 'h_ph', 
                   'n_seg_ph', key_field = 'h_ind', classfield = 'classification')
    
    # # Unstack canopy_rh and rename (possible to drop canopy_rh_abs)
    # for i in range(0,len(percentile_intervals)):
    #     arr = []
    #     for j in range(0,len(calc_df)):
    #         if type(calc_df.iloc[j].canopy_rh) == np.ndarray:
    #             arr.append(bin_df.iloc[j].canopy_rh[i])
    #         else:
    #             arr.append(0)
            
        # arr = np.array(arr)
        # title = 'canopy_rh_' + str(percentile_intervals[i])
        # bin_df[title] = arr

    # # Unstack canopy_rh and rename (possible to drop canopy_rh_abs)
    # for i in range(0,len(percentile_intervals)):
    #     arr = []
    #     for j in range(0,len(bin_df)):
    #         if type(bin_df.iloc[j].canopy_rh_abs) == np.ndarray:
    #             arr.append(bin_df.iloc[j].canopy_rh_abs[i])
    #         else:
    #             arr.append(0)
            
        # arr = np.array(arr)
        # title = 'canopy_rh_abs_' + str(percentile_intervals[i])
        # bin_df[title] = arr

    bin_df['canopy_h'] = bin_df['canopy_rh_98'] 
    bin_df['h_canopy_abs'] = bin_df['canopy_rh_abs_98'] #        
    
    # bin_df = bin_df.drop(columns=['canopy_rh_98'])
    # bin_df = bin_df.drop(columns=['canopy_rh_abs_98'])
        
    # # Compute slope
    # calc_df['terrain_slope'] = calc_df.h_te_rise / calc_df.h_te_run
    # calc_df = calc_df.drop(columns=['h_te_rise','h_te_run'])


    # Compute h_dif_can
    bin_df['h_dif_canopy'] = bin_df.canopy_h - bin_df.canopy_rh_50
    
    # Compute traditional canopy cover
    bin_df['traditional_canopy_cover'] = (bin_df['n_ca_photons']) /\
        (bin_df['n_ca_photons']  + bin_df['n_te_photons'])
            
    # Compute Rumple Index (canopy surface area/ground surface area)
    sub_bin = sub_bin_canopy_metrics_truth(truth_swath, res_at = 2)
    veg_bin = calc_rumple_index(sub_bin)
    bin_df = bin_df.join(veg_bin)
            
    # Compute Area-Under-Canopy/Open Sky Ratio
    bin_df['canopy_ratio'] = bin_df.veg_area / (bin_df.h_max_canopy * 
                                                  bin_df.ground_len)
    
    # Merge Computed columns to bin_df
    bin_df = bin_df.set_index('h_ind')

    bin_df = slope_df(truth_swath, bin_df, ground_class)
    
    # Compute h_dif_ref
    # bin_df['h_dif_ref'] = bin_df.h_te_median - bin_df.dem_h
    
    # Height bins (%)
    h_bin_list = [[0,1],[1,2.5],[2.5,5],[5,7.5],[7.5,10],[10,12.5],[12.5,15],
                  [15,17.5],[17.5,20],[20,22.5],[22.5,25],[25,27.5],[27.5,30],
                  [30,200]]
    h_base_val = 0
    h_bin_prefix = 'h_bin_'
    for i in range(0,len(h_bin_list)):
        title = h_bin_prefix + str(i)
        bin_df = computer_vertical_density_p(truth_swath, bin_df, 
                        h_bin_list[i][0], h_bin_list[i][1], title, veg_class = [3,4,5],
                        h_base = h_base_val)

    # Compute photons above set threshold
    bin_df = get_n_photons_above_threshold(truth_swath, bin_df, 
                                           classes = [3,4,5], 
                                           h_base = h_base_val)

    if 'date' in truth_swath.columns:
        bin_df['year'] = truth_swath.date.iloc[0].year
        bin_df['month'] = truth_swath.date.iloc[0].month
        bin_df['day'] = truth_swath.date.iloc[0].day
    return bin_df

def match_truth_fields(truth_bin, atl08_bin):
    # Set core fields to be equal to atl03_bin
    drop_fields = ['alongtrack', 'crosstrack', 'easting', 
                           'northing', 'lon', 'lat']
    
    truth_bin = truth_bin.loc[:, ~truth_bin.columns.isin(drop_fields)]
    
    include_fields = ['h_ind','alongtrack', 'crosstrack', 'easting', 
                      'northing', 'longitude', 'latitude', 'delta_time','time']
    
    atl08_bin = atl08_bin.loc[:, atl08_bin.columns.isin(include_fields)]
    
    truth_bin = pd.merge(truth_bin, atl08_bin, on = 'h_ind')
    
    return truth_bin

def match_atl_fields(truth_bin, atl08_bin):
    # Set core fields to be equal to atl03_bin
    drop_fields = ['h_ind','alongtrack', 'crosstrack', 'easting', 
                      'northing', 'longitude', 'latitude', 'delta_time','time']
    
    truth_bin = truth_bin.loc[:, ~truth_bin.columns.isin(drop_fields)]
    
    include_fields = ['h_ind','alongtrack', 'crosstrack', 'easting', 
                      'northing', 'longitude', 'latitude', 'delta_time','time']
    
    atl08_bin = atl08_bin.loc[:, atl08_bin.columns.isin(include_fields)]
    
    truth_bin = pd.merge(truth_bin, atl08_bin, on = 'h_ind')
    
    return truth_bin
    
# def main(in_atl03, in_atl08, output_dir, res, rando, v):
#     try:
#         os.mkdir(output_dir)
#     except:
#         print('Folder already exists')
    
#     atl03_list = os.listdir(in_atl03)

#     if rando is True:
#         print('Random sort ATL03 List')
#         random.shuffle(atl03_list)

#     # # Inputs
#     # res = 30
#     res_field = 'alongtrack'
#     for i in range(0,len(atl03_list)):
#         print(atl03_list[i])
#         atl03file = atl03_list[i]
#         atl08file = 'ATL08' + atl03_list[i].split('ATL03')[1]
#         atl03filepath =  os.path.join(in_atl03, atl03file)
#         atl08filepath =  os.path.join(in_atl08, atl08file)
#         gt_list = ['gt1r','gt1l','gt2l','gt2r','gt3l','gt3r']
#         for gt in gt_list:
#             try:
#                 csv_file = os.path.join(output_dir, atl08file.split('.')[0] + '_' + gt + '30m.csv')
#                 if os.path.exists(csv_file) == False:
#                     bin_df = rebin_atl08(atl03filepath, atl08filepath, gt, res, res_field)
#                     bin_df.to_csv(csv_file)
#                 else:
#                     print('File already exists')
#             except:
#                 print('Binning failed')
    
def main(in_atl03, in_atl08, output_dir, res, rando, v):
    try:
        os.mkdir(output_dir)
    except:
        print('Folder already exists')
    
    atl03_list = os.listdir(in_atl03)

    if rando is True:
        print('Random sort ATL03 List')
        random.shuffle(atl03_list)

    # # Inputs
    # res = 30
    res_field = 'alongtrack'
    for i in range(0,len(atl03_list)):
        print(atl03_list[i])
        atl03file = atl03_list[i]
        atl08file = 'ATL08' + atl03_list[i].split('ATL03')[1]
        atl03filepath =  os.path.join(in_atl03, atl03file)
        atl08filepath =  os.path.join(in_atl08, atl08file)
        gt_list = ['gt1r','gt1l','gt2l','gt2r','gt3l','gt3r']
        for gt in gt_list:
            try:
                csv_file = os.path.join(output_dir, atl08file.split('.')[0] + '_' + gt + '30m.csv')
                if os.path.exists(csv_file) == False:
                    print('Generate ATL03 Struct')
                    atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath)
                    # atl03.df = atl03.df[atl03.df.lat_ph > 57]
                    # atl03.df = atl03.df[atl03.df.lat_ph < 65]
                    #atl03.df = atl03.df[atl03.df.lon_ph > 88]
                    #atl03.df = atl03.df[atl03.df.lon_ph < 106]
                    
                    print('Generate ATL08 Struct')
                    atl08 = get_atl08_struct(atl08filepath, gt, atl03)
                    # atl08.df = atl08.df[atl08.df.latitude > 57]
                    # atl08.df = atl08.df[atl08.df.latitude < 65]
                    #atl08.df = atl08.df[atl08.df.longitude > 88]
                    #atl08.df = atl08.df[atl08.df.longitude < 106]
                    
                    # Rebin ATL08
                    bin_df = rebin_atl08(atl03, atl08, gt, res, res_field)
                    bin_df.to_csv(csv_file)
                else:
                    print('File already exists')
            except:
                print('Binning failed')
    
if __name__ == '__main__':
    """ Command line entry point """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("atl03_dir",
                        help="Input ATL03 file or ATL03 directory")

    parser.add_argument("atl08_dir",
                        help="Input ATL08 file or ATL08 directory")

    parser.add_argument("out_dir",
                        help="Output directory")

    parser.add_argument("-r", "--res", nargs='?', const=1, type=int, default=30,
                        help="Alongtrack resolution (m)")

    parser.add_argument('-rando', '--random', action='store_true',
                        help='Randomly iterate through directory')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose Output')

    args = parser.parse_args()

    args.atl03_dir = os.path.normpath(args.atl03_dir)
    args.atl08_dir = os.path.normpath(args.atl08_dir)
    args.out_dir = os.path.normpath(args.out_dir)

    main(args.atl03_dir,
        args.atl08_dir,
        args.out_dir,
        args.res,
        args.rando,
        args.verbose)
