#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 13:45:53 2020

@author: eguenther
"""
import os
import ctypes
from numpy.ctypeslib import ndpointer 
import pandas as pd
import numpy as np
from scipy import stats
from scipy import interpolate

root = os.path.dirname(__file__)
superFilterFile_windows = os.path.join(root,'closest_x64.dll')
superFilterFile_linux = os.path.join(root,'phorealc.so')


def indexMatch(measuredArray,truthArray):
    # print("Match corresponding indices...", end = " ")
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
    # print("Complete")
    return np.array(C).astype(int)

def get_max100(series):
    try:
        max98 = np.percentile(series, 100)
    except:
        max98 = np.nan
    return max98

def get_max98(series):
    try:
        max98 = np.percentile(series, 98)
    except:
        max98 = np.nan
    return max98

def get_max75(series):
    try:
        max75 = np.percentile(series, 75)
    except:
        max75 = np.nan
    return max75

def get_max50(series):
    try:
        max50 = np.percentile(series, 50)
    except:
        max50 = np.nan
    return max50

def get_max25(series):
    try:
        max25 = np.percentile(series, 25)
    except:
        max25 = np.nan
    return max25

def get_min(series):
    try:
        min0 = np.percentile(series, 0)
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
        mode = stats.mode(np.array(series))[0][0]
    except:
        mode = np.nan
    return mode

def calculate_seg_meteric(df_in, df_out, classification, operation, field, 
                          outfield, key_field = 'bin_id',
                          classfield = 'classification'):
    df_filter = df_in[df_in[classfield].isin(classification)]
    zgroup = df_filter.groupby(key_field)
    zout = zgroup.aggregate(operation)
    zout[key_field] = zout.index
    zout = zout.reset_index(drop = True)
    # zout['segment_id_beg'] = zout['seg_id']
    zout[outfield] = zout[field]
    zout = zout.filter([outfield,key_field])
    df_out = df_out.merge(zout, on=key_field,how='left')  
    return df_out

def create_key_df(df, field, res):
    start = np.min(df[field])
    stop = np.max(df[field])
    beg = np.arange(start,stop,res)
    mid = beg + (res / 2)
    end = beg + (res)
    res_field = end - beg
    key_df = pd.DataFrame(beg,columns=['beg_id'])
    key_df = pd.concat([key_df,pd.DataFrame(mid
                                        ,columns=['mid_id'])],axis=1)
    key_df = pd.concat([key_df,pd.DataFrame(
        end,columns=['end_id'])],axis=1)
    key_df = pd.concat([key_df,pd.DataFrame(
        res_field,columns=['res'])],axis=1)
    key_df = pd.concat([key_df,pd.DataFrame(
        mid,columns=['bin_id'])],axis=1)
    return key_df

def get_target_keys(key_df, df, field, res = -1):
 
    if len(np.unique(key_df.res)) == 1:
        res = key_df.res[0]
    elif res != -1:
        res = res
    else:
        res = key_df.res[0]
    # res = 1
    mid = np.array(df[field])

    target_mid = np.array(key_df.mid_id)

    minInd = indexMatch(target_mid, mid)
    
    datalen = len(target_mid)
    minInd[minInd >= datalen] = (datalen - 1)
    include = np.zeros(len(mid))
    seg_at_comp = np.array([target_mid[x] for x in minInd])
    # seg_id_truth = np.array([seg_id[x] for x in index])
    diff = np.abs(seg_at_comp - mid)
    
    include[diff <= (res/2)] = 1    
    
    target_key = np.array(key_df.bin_id[minInd])
    
    return target_key, include

def agg_keys(key_df, df, agg_list, key = 'bin_id'):
    for i in agg_list:
        agg = i.split(',')
        if 'perfect' in agg[2]:
            c = 'perfect_class'
        else:
            c = 'classification'

        class_str = agg[4]
        class_list = class_str.strip('][').split(',')
        if 'median' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                                             np.median,agg[2],agg[1], 
                                             key_field = key,
                                             classfield = c)
        elif 'count' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         np.size,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'count_unique' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_len_unique,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'range' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         np.ptp,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max100' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max100,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max98' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max98,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max75' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max75,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max50' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max50,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max25' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                        get_max25,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'min' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_min,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'mode' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                          get_mode,agg[2],agg[1], key_field = key,
                                              classfield = c)        
        elif 'std' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_std,agg[2],agg[1], key_field = key,
                                             classfield = c) 
        elif 'rh_canopy' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, [1], 
                         np.median,'h_ph','ground', key_field = key,
                                             classfield = c) 
            key_df = calculate_seg_meteric(df, key_df, [2,3], 
                         np.median,'h_ph','canopy', key_field = key,
                                             classfield = c) 
            key_df[agg[1]] = key_df.canopy - key_df.ground
            key_df[agg[1]][key_df[agg[1]] < 0] = 0
            key_df[agg[1]][key_df[agg[1]] > 130] = np.nan
            key_df = key_df.drop(columns=['ground','canopy'])
        elif 'radiometry' in agg[3]:
            
            key_df = calculate_seg_meteric(df, key_df, [-1,0,1,2,3], 
                         get_len_unique,agg[2],'unique_time', key_field = key,
                                             classfield = c) 
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         np.size,agg[2],'target_time', key_field = key,
                                             classfield = c) 
            t = np.asarray(key_df.target_time)
            a = np.asarray(key_df.unique_time)
            key_df[agg[1]] = np.divide(t,a,out=np.zeros_like(t),where=a!=0)
            key_df = key_df.drop(columns=['unique_time','target_time'])
    return key_df

def orient_df(df, field):
    df = df.sort_values(field)
    df = df.reset_index(drop = True)
    return df

def get_bin_df(df, unit, res, agg_list):
    '''
    

    Parameters
    ----------
    df : Pandas DF
        Input DF with values that need to be aggregated.
    unit : String (key field)
        Field name which key needs to be calculated (alongtrack, seg, etc.)
    res : Float (Bin size)
        Resolution of bin given a specified unit.
    agg_list : List of strings
        Define your agg_list.  It is a list of strings on how to aggregate.
            The strings should be arranged so that:
                1. The type of file (ATL03, ATL08, Truth, Skinny)
                2. The name of the output field
                3. The field to aggregate on
                4. The operation to use
                5. The classes to run the operation on
            Each item is to be separated by a comma like:
                'skinny,atl03_ground_median,h_ph,median,[1]'

    Returns
    -------
    df_bin : Pandas DF
        DF with binned values.

    '''
    df = orient_df(df,unit)
    key_df = create_key_df(df, unit, res)
    target_key, include = get_target_keys(key_df, df, unit)
    df = df.reset_index(drop = True)
    df = pd.concat([df,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)
    df = pd.concat([df,pd.DataFrame(include,columns=['include'])],axis=1)
    df = df[df.include == 1]
    df_bin = agg_keys(key_df, df, agg_list, key = 'bin_id')
    return df_bin
    

def get_ATL08_seg_mid(atl08_df):
    '''
    Description
    ----------    
    Get the keys for ATL08 upsampling and downsampling

    Parameters
    ----------
    atl08 : ATL08 Struct
        ATL08 Struct that contained ATL08 Pandas DF.
    geolocation : Pandas DF
        Geolocation Pandas DF associated with ATL08 file.

    Returns
    -------
    key_df : Pandas DF
        Pandas DF used of keys for ATL08 Segment Analysis.

    '''    

    res = 5
    beg = (np.unique(np.array(atl08_df.segment_id_beg)))
    beg = beg[~np.isnan(beg)].astype(int)
    
    # Find 'target' type and set it in the middle 

    mid = beg + ((res)/2)
    # atl08.df.reset_index()
    # atl08.df = pd.concat([df,pd.DataFrame(mid,
    #                                       columns=['segment_id_mid'])],axis=1)
    
    return mid

def get_specific_key(data, field):
    # If converting ATL03
    specific_field = field.lower()
    if data.lower() == 'atl03':
        if field.lower() == 'seg':
            specific_field = 'segment_id'
        if field.lower() == 'z':
            specific_field = 'h_ph'
        if field.lower() == 'z':
            specific_field = 'h_ph'
    if data.lower() == 'atl08':
        if field.lower() == 'seg' or field.lower() == 'segment_id':
            specific_field = 'segment_id_mid'

    return specific_field
    

# def get_multibin_df(atl03, atl08, truth, unit, res, agg_list):
#     if atl03:
#         key_df = create_key_df(atl03.df, unit, res)
#         target_key, include = get_target_keys(key_df, atl03.df, unit)
#         df = atl03.df.reset_index(drop = True)
#         df = pd.concat([df,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)
#         df = pd.concat([df,pd.DataFrame(include,columns=['include'])],axis=1)
#         df = df[df.include == 1]
#         in_agg_list = []
#         for agg in agg_list:
#             if agg.split(',')[0] == 'atl03':
#                 in_agg_list.append(agg)
#         df_bin = agg_keys(key_df, df, in_agg_list, key = 'bin_id')
#     if atl08:
#         atl08_unit = get_specific_key('atl08', unit)
#         if atl08_unit.lower() == 'segment_id_mid':
#             seg_mid = get_ATL08_seg_mid(atl08.df)
#             atl08.df.reset_index(drop = True)
#             atl08.df = pd.concat([atl08.df,pd.DataFrame(seg_mid,
#                                         columns=['segment_id_mid'])],axis=1)            
#         target_key, include = get_target_keys(key_df, atl08.df, atl08_unit)
#         df = atl08.df.reset_index(drop = True)
#         df = pd.concat([df,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)
#         df = pd.concat([df,pd.DataFrame(include,columns=['include'])],axis=1)
#         df = df[df.include == 1]
#         in_agg_list = []
#         for agg in agg_list:
#             if ((agg.split(',')[0] == 'atl08') and \
#                 (agg.split(',')[1] == 'all') and (res <= 5)):
#                 res = 5
#                 df_bin = pd.merge(key_df, df, on="segment_id_mid",how='left')
#             elif agg.split(',')[0] == 'atl08':
#                 in_agg_list.append(agg)
#                 df_bin = agg_keys(key_df, df, in_agg_list, key = 'bin_id')
#     if truth:
#         target_key, include = get_target_keys(key_df, atl08.df, unit)
#         df = atl03.df.reset_index(drop = True)
#         df = pd.concat([df,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)
#         df = pd.concat([df,pd.DataFrame(include,columns=['include'])],axis=1)
#         df = df[df.include == 1]
#         in_agg_list = []
#         for agg in agg_list:
#             if agg.split(',')[0] == 'truth':
#                 in_agg_list.append(agg)
#         df_bin = agg_keys(key_df, df, in_agg_list, key = 'bin_id')        
    
#     return df_bin

def interpolate_domain(atl08_at, atl08_domain, key_df_at):
    intep_func = interpolate.interp1d(atl08_at, atl08_domain, kind='linear', 
                              fill_value='extrapolate')
    interp_domain = intep_func(key_df_at)
    return interp_domain

def compute_relative_canopy_height(canopy, ground):
    rh = canopy - ground
    rh[rh < 0] = 0
    return rh

def create_atl08_bin(atl03, atl08, res_at = 30):
    #Create Key DF based on ATL08
    atl03_df = orient_df(atl03.df, 'alongtrack')
    atl08_df = orient_df(atl08.df, 'alongtrack')
    
    key_df = create_key_df(atl08_df, 'alongtrack', res_at)

    domain_list = ['time','easting','northing','crosstrack','alongtrack',
               'delta_time','delta_time_beg','delta_time_end','latitude',
               'longitude','segment_id_beg','segment_id_end']
    
    for domain in domain_list:
        domain_interp = interpolate_domain(atl08_df.alongtrack, 
                                           atl08_df[domain], key_df.bin_id)
        key_df = pd.concat([key_df,pd.DataFrame(domain_interp,
                                                columns=[domain])],axis=1)

    key_id = key_df.bin_id
    atl08_id = atl08_df.alongtrack
    
    minInd = indexMatch(atl08_id,key_id)
    res = 50
    datalen = len(atl08_id)
    minInd[minInd >= datalen] = (datalen - 1)
    include = np.zeros(len(key_id))
    seg_at_comp = np.array([atl08_id[x] for x in minInd])
    # seg_id_truth = np.array([seg_id[x] for x in index])
    diff = np.abs(seg_at_comp - key_id)
    
    include[diff <= (res)] = 1  
    
    atl08_df = atl08_df.loc[minInd]
    atl08_df = atl08_df.drop(columns = domain_list)


    df_bin_08 = pd.concat([key_df.reset_index(drop=True),
                     atl08_df.reset_index(drop=True)], axis=1)
    
    
    df_bin_08 = df_bin_08[include == 1]
    
    #Compute the ATL03 
    
    agg_list = ['atl03,atl03_ground_median,h_ph,median,[1]',
            'atl03,atl03_canopy_h,h_ph,get_max98,[2,3]']
    
    target_key, include = get_target_keys(key_df, atl03_df, 'alongtrack')
    atl03_df = atl03_df.reset_index(drop = True)
    atl03_df = pd.concat([atl03_df,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)
    atl03_df = pd.concat([atl03_df,pd.DataFrame(include,columns=['include'])],axis=1)
    atl03_df = atl03_df[atl03_df.include == 1]
    
    
    
    df_bin_03 = agg_keys(key_df, atl03_df, agg_list, key = 'bin_id')
    
    atl03_canopy_rh =\
        compute_relative_canopy_height(df_bin_03.atl03_canopy_h, 
                                       df_bin_03.atl03_ground_median)
        
    df_bin_03['atl03_canopy_rh'] = atl03_canopy_rh

    atl03_df = pd.concat([atl03_df,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)

    
    df_bin_03.drop(df_bin_03.columns.difference(['bin_id','atl03_ground_median',
                                         'atl03_canopy_h','atl03_canopy_rh']),
                   1,inplace=True)
    
    
    
    df_bin = pd.merge(df_bin_08, df_bin_03, on="bin_id",how='left')    
    
    df_bin.drop(columns = ['beg_id', 'mid_id', 'end_id', 'res'])
    
    return df_bin
    
    

if __name__ == "__main__":
    if os.name == 'nt':
        basepath = 'Y:/USERS/eric/2_production/'
    else:
        basepath = '/LIDAR/server/USERS/eric/2_production/'

    atl03file = 'ATL03_20181021130238_03500103_002_01.h5'
    atl08file = 'ATL08_20181021130238_03500103_002_01.h5'
    
    atl03filepath = basepath + atl03file
    atl08filepath = basepath + atl08file
    gt = 'gt1r'    

    from icesatReader import get_atl03_struct
    from icesatReader import get_atl08_struct
    
    
    print('ATL03 Heights')
    atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath)
    
    print('ATL08 Land Segments')
    atl08 = get_atl08_struct(atl08filepath, gt, atl03)
    
    
    upsampled_atl08_bin = create_atl08_bin(atl03, atl08, res_at = 30)
    
    agg_list = ['atl03,atl03_min,h_ph,min,[1]']

    downsampled_atl03_bin = get_bin_df(atl03.df, 'time', 0.01, agg_list)
