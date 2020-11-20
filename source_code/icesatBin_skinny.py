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
        length = int(len(series))
    except:
        length = np.nan
    return length

def get_len_unique(series):
    try:
        length = int(len(np.unique(series)))
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

def get_mean(series):
    try:
        m = np.mean(series)
    except:
        m = np.nan
    return m

# Can be included if 'from scipy import stats' is used
# def get_mode(series):
#     try:
#         mode = stats.mode(np.array(series))[0][0]
#     except:
#         mode = np.nan
#     return mode

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
    column = np.array(df[field])
    start = column[0]
    stop = column[len(column)-1]
    if start < stop:      
        beg = np.arange(start,stop,res)
        mid = beg + (res / 2)
        end = beg + (res)
        res_field = end - beg
    else:
        beg = np.arange(start,stop,(-1*res))
        mid = beg - (res / 2)
        end = beg - (res)
        res_field = np.abs(end - beg)  
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

def get_target_keys(key_df, df, field):
 
    if len(np.unique(key_df.res)) == 1:
        res = key_df.res[0]
    else:
        res = key_df.res[0]
    mid = np.array(df[field])

    target_mid = np.array(key_df.mid_id)
    # if 'seg' in field:
    #     beg = beg[~np.isnan(beg)].astype(int)
    #     mid = beg + ((res)/2)
    # else:
    #     beg = beg[~np.isnan(beg)].astype(float)
    #     mid = beg
 
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
        agg = i.split(';')
        if 'perfect' in agg[2]:
            c = 'perfect_class'
        elif('bathy' in agg[0]):
            c = 'Label'
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
                         np.size,agg[2],agg[1], key_field = key,
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
        # elif 'mode' in agg[3]:
        #     key_df = calculate_seg_meteric(df, key_df, class_list, 
        #                  get_mode,agg[2],agg[1], key_field = key,
        #                                      classfield = c)        
        elif 'std' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_std,agg[2],agg[1], key_field = key,
                                             classfield = c) 
        elif 'mean' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_mean,agg[2],agg[1], key_field = key,
                                             classfield = c) 
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
                'skinny,atl03_ground_median,h_ph,median,[1]

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
    

if __name__ == "__main__":
    # Stuff to load a df for an example
    if os.name == 'nt':
        example_file = 'Y:/USERS/eric/1_experiment/thin_file_gui/test.pkl'
    else:
        example_file = '/LIDAR/server/poseidon_files/USERS/eric/' +\
            '1_experiment/thin_file_gui/test.pkl'

    df = pd.read_pickle(example_file)
    
    '''
    Define your agg_list.  It is a list of strings on how to aggregate.
    The strings should be arranged so that:
        1. The type of file (ATL03, ATL08, Truth, Skinny)
        2. The name of the output field
        3. The field to aggregate on
        4. The operation to use
        5. The classes to run the operation on
    Each item is to be separated by a comma like:
        'skinny,atl03_ground_median,h_ph,median,[1]'
    '''
    agg_list = ['skinny,atl03_ground_median,h_ph,median,[1]',
                'skinny,atl03_canopy_max98,h_ph,get_max98,[2,3]']
    

    df_bin = get_bin_df(df, 'alongtrack', 10, agg_list)