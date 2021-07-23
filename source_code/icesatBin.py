#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 13:45:53 2020

@author: eguenther
"""
import os
import argparse
import ctypes
from shutil import copyfile    
from numpy.ctypeslib import ndpointer 
import pandas as pd
import numpy as np
try:
    from scipy import stats
except:
    print('scipy.stats import failed')
from scipy import interpolate
import h5py
import time


from icesatReader import get_atl03_struct
from icesatReader import get_atl08_struct


root = os.path.dirname(__file__)
superFilterFile_windows = os.path.join(root,'closest_x64.dll')
superFilterFile_linux = os.path.join(root,'phorealc.so')

def indexMatch(measuredArray,truthArray):
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
        max98 = np.nanmax(series)
    except:
        max98 = np.nan
    return max98

def get_max98(series):
    try:
        max98 = np.percentile(series, 98)
    except:
        max98 = np.nan
    return max98

def get_max95(series):
    try:
        max95 = np.percentile(series, 95)
    except:
        max95 = np.nan
    return max95

def get_max90(series):
    try:
        max90 = np.percentile(series, 90)
    except:
        max90 = np.nan
    return max90

def get_max85(series):
    try:
        max85 = np.percentile(series, 85)
    except:
        max85 = np.nan
    return max85

def get_max80(series):
    try:
        max80 = np.percentile(series, 80)
    except:
        max80 = np.nan
    return max80

def get_max75(series):
    try:
        max75 = np.percentile(series, 75)
    except:
        max75 = np.nan
    return max75

def get_max70(series):
    try:
        max70 = np.percentile(series, 70)
    except:
        max70 = np.nan
    return max70

def get_max65(series):
    try:
        max65 = np.percentile(series, 65)
    except:
        max65 = np.nan
    return max65

def get_max60(series):
    try:
        max60 = np.percentile(series, 60)
    except:
        max60 = np.nan
    return max60

def get_max55(series):
    try:
        max55 = np.percentile(series, 55)
    except:
        max55 = np.nan
    return max55

def get_max50(series):
    try:
        max50 = np.percentile(series, 50)
    except:
        max50 = np.nan
    return max50

def get_max45(series):
    try:
        max45 = np.percentile(series, 45)
    except:
        max45 = np.nan
    return max45

def get_max40(series):
    try:
        max40 = np.percentile(series, 40)
    except:
        max40 = np.nan
    return max40

def get_max35(series):
    try:
        max35 = np.percentile(series, 35)
    except:
        max35 = np.nan
    return max35

def get_max30(series):
    try:
        max30 = np.percentile(series, 30)
    except:
        max30 = np.nan
    return max30

def get_max25(series):
    try:
        max25 = np.percentile(series, 25)
    except:
        max25 = np.nan
    return max25

def get_max20(series):
    try:
        max20 = np.percentile(series, 20)
    except:
        max20 = np.nan
    return max20

def get_max15(series):
    try:
        max15 = np.percentile(series, 15)
    except:
        max15 = np.nan
    return max15

def get_max10(series):
    try:
        max10 = np.percentile(series, 10)
    except:
        max10 = np.nan
    return max10

def get_max5(series):
    try:
        max5 = np.percentile(series, 5)
    except:
        max5 = np.nan
    return max5

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
    
def get_mean(series):
    try:
        m = np.mean(series)
    except:
        m = np.nan
    return m

def above_op(series):
    q_list = [0,10,20,25,30,40,50,60,70,75,80,90,98,100]
    try:
        all_percent = [np.percentile(series, q_list)]
    except:
        # all_percent = np.nan
        all_percent = np.array([np.nan, np.nan,np.nan,np.nan,np.nan,np.nan,
                                np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,
                                np.nan,np.nan])
    return all_percent

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
        # print(i)
        dt2 = time.time()
       
        agg = i.split(';')
        if 'perfect' in agg[2]:
            c = 'perfect_class'
        elif('bathy' in agg[0]):
            c = 'Label'
        else:
            c = 'classification'

        class_str = agg[4]
        class_list_str = class_str.strip('][').split(',')
        class_list = []
        if len(class_list_str) > 0:
            for j in range(0,len(class_list_str)):
                class_list.append(int(class_list_str[j]))
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
        elif 'max95' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max95,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max90' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max90,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max85' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max85,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max80' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max95,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max75' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max75,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max70' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max70,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max65' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max65,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max60' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max60,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max55' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max55,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max50' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max50,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max45' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max45,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max40' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max40,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max35' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max35,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max30' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max30,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max25' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                        get_max25,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max20' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max20,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max15' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                        get_max15,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max10' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_max10,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'max5' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                        get_max5,agg[2],agg[1], key_field = key,
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
        elif 'mean' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_mean,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'get_len' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, class_list, 
                         get_len,agg[2],agg[1], key_field = key,
                                             classfield = c)
        elif 'rh_canopy' in agg[3]:
            key_df = calculate_seg_meteric(df, key_df, [1], 
                         np.median,agg[2],'ground', key_field = key,
                                             classfield = c) 
            key_df = calculate_seg_meteric(df, key_df, [2,3], 
                         np.median,agg[2],'canopy', key_field = key,
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
        elif 'above_percentile' in agg[3]:
            d3 = time.time()
            q_list = [10,20,25,30,40,50,60,70,75,80,90,98,100]
            df_filter = df[df[c].isin(class_list)]
            zgroup = df_filter.groupby(key)
            zout = zgroup.aggregate(above_op)
            zout[key] = zout.index
            zout = zout.reset_index(drop = True)
            field = agg[2]
            outfield = 'outfield'
            # zout['segment_id_beg'] = zout['seg_id']
            zout[outfield] = zout[field]
            zout = zout.filter([outfield,key])
            print(zout.outfield[0][0])
            print(time.time() - d3)
            for i in range(0,len(q_list)):
                arr = []
                for j in range(0,len(zout.outfield)):
                    arr.append(zout.outfield[j][0][i])
                title = agg[1] + str(q_list[i])
                print(title)
                zout[title] = np.array(arr)
            zout.drop(columns=['outfield'], inplace=True)
            key_df = key_df.merge(zout, on=key,how='left')  
        # print(time.time() - dt2)
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
    
def interpolate_domain(atl08_at, atl08_domain, key_df_at):
    intep_func = interpolate.interp1d(atl08_at, atl08_domain, kind='linear', 
                              fill_value='extrapolate')
    interp_domain = intep_func(key_df_at)
    return interp_domain

def compute_rh(canopy, ground):
    rh = canopy - ground
    rh[canopy.isnull()] = 0
    rh[canopy.isnull() & ground.isnull()] = np.nan
    return rh

def create_atl08_bin(atl03, atl08, res_at = 30):
    import time
    #Create Key DF based on ATL08
    dt = time.time()
    atl03_df = orient_df(atl03.df, 'alongtrack')
    atl08_df = orient_df(atl08.df, 'alongtrack')
    
    key_df = create_key_df(atl08_df, 'alongtrack', res_at)

    domain_list = [
               'delta_time','delta_time_beg','delta_time_end','latitude',
               'longitude','segment_id_beg','segment_id_end']
    
    for domain in domain_list:
        domain_interp = interpolate_domain(atl08_df.alongtrack, 
                                           atl08_df[domain], key_df.bin_id)
        key_df = pd.concat([key_df,pd.DataFrame(domain_interp,
                                                columns=[domain])],axis=1)

    key_id = key_df.bin_id
    atl08_id = atl08_df.alongtrack
    print('Part1')
    print(time.time() - dt)
    dt = time.time()

    
    minInd = indexMatch(atl08_id,key_id)
    res = res_at
    datalen = len(atl08_id)
    minInd[minInd >= datalen] = (datalen - 1)
    include = np.zeros(len(key_id))
    seg_at_comp = np.array([atl08_id[x] for x in minInd])
    diff = np.abs(seg_at_comp - key_id)
    
    include[diff <= (res)] = 1  
    
    atl08_df = atl08_df.loc[minInd]
    atl08_df = atl08_df.drop(columns = domain_list)

    print('Part2')
    print(time.time() - dt)
    dt = time.time()

    df_bin_08 = pd.concat([key_df.reset_index(drop=True),
                     atl08_df.reset_index(drop=True)], axis=1)
    
    
    df_bin_08 = df_bin_08[include == 1]
    
    #Compute the ATL03 
    
    print(time.time() - dt)
    dt = time.time()
    print('Part3')
    
    agg_list = ['atl03;atl03_ground_median;h_ph;median;[1]',
            'atl03;atl03_canopy_h;h_ph;get_max98;[2,3]',                
            'atl03;gedi_rh_;h_ph;above_percentile;[1,2,3]',
            'atl03;atl03_rh_;h_ph;above_percentile;[2,3]',
            'atl03;atl03_n_unclass;h_ph;get_len;[0]',
            'atl03;atl03_n_ground;h_ph;get_len;[1]',
            'atl03;atl03_n_canopy;h_ph;get_len;[2]',
            'atl03;atl03_n_hi_canopy;h_ph;get_len;[3]',
            'atl03;atl03_radiometry_ground;h_ph;radiometry;[1]',
            'atl03;atl03_radiometry_canopy;h_ph;radiometry;[2,3]',
            'atl03;atl03_radiometry_all;h_ph;radiometry;[1,2,3]']
    
    target_key, include = get_target_keys(key_df, atl03_df, 'alongtrack')
    atl03_df = atl03_df.reset_index(drop = True)
    atl03_df = pd.concat([atl03_df,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)
    atl03_df = pd.concat([atl03_df,pd.DataFrame(include,columns=['include'])],axis=1)
    atl03_df = atl03_df[atl03_df.include == 1]
    
    print(time.time() - dt)
    dt = time.time()
    print('Part4a')
        
    df_bin_03 = agg_keys(key_df, atl03_df, agg_list, key = 'bin_id')

    print(time.time() - dt)
    dt = time.time()
    print('Part4b')

    atl03_df = pd.concat([atl03_df,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)

    print(time.time() - dt)
    dt = time.time()
    print('Part5')
    
    df_bin_03.drop(df_bin_03.columns.difference(['bin_id','atl03_ground_median',
                                                 'atl03_canopy_h',
                                          'gedi_rh_100','gedi_rh_98',
                                           'gedi_rh_90','gedi_rh_80',
                                           'gedi_rh_75','gedi_rh_70',
                                           'gedi_rh_60','gedi_rh_50',
                                           'gedi_rh_40','gedi_rh_30',
                                           'gedi_rh_25','gedi_rh_20',
                                           'gedi_rh_10','atl03_rh_100',
                                           'atl03_rh_98','atl03_rh_90',
                                           'atl03_rh_80','atl03_rh_75',
                                           'atl03_rh_70','atl03_rh_60',
                                           'atl03_rh_50','atl03_rh_40',
                                           'atl03_rh_30','atl03_rh_25',
                                           'atl03_rh_20','atl03_rh_10',
                                         'atl03_n_unclass','atl03_n_ground',
                                         'atl03_n_canopy','atl03_n_hi_canopy'
                                         ]),
                   1,inplace=True)
        
    print(time.time() - dt)
    dt = time.time()
    print('Part6')
    
    df_bin = pd.merge(df_bin_08, df_bin_03, on="bin_id",how='left')       
    
    df_bin['gedi_rh_100'] = compute_rh(df_bin['gedi_rh_100'], 
                                        df_bin['atl03_ground_median'])
    df_bin['gedi_rh_98'] = compute_rh(df_bin['gedi_rh_98'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_90'] = compute_rh(df_bin['gedi_rh_90'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_80'] = compute_rh(df_bin['gedi_rh_80'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_75'] = compute_rh(df_bin['gedi_rh_75'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_70'] = compute_rh(df_bin['gedi_rh_70'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_60'] = compute_rh(df_bin['gedi_rh_60'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_50'] = compute_rh(df_bin['gedi_rh_50'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_40'] = compute_rh(df_bin['gedi_rh_40'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_30'] = compute_rh(df_bin['gedi_rh_30'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_25'] = compute_rh(df_bin['gedi_rh_25'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_20'] = compute_rh(df_bin['gedi_rh_20'],
                                      df_bin['atl03_ground_median'])
    df_bin['gedi_rh_10'] = compute_rh(df_bin['gedi_rh_10'],
                                      df_bin['atl03_ground_median'])
    df_bin['atl03_rh_100'] = compute_rh(df_bin['atl03_rh_100'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_98'] = compute_rh(df_bin['atl03_rh_98'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_90'] = compute_rh(df_bin['atl03_rh_90'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_80'] = compute_rh(df_bin['atl03_rh_80'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_75'] = compute_rh(df_bin['atl03_rh_75'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_70'] = compute_rh(df_bin['atl03_rh_70'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_60'] = compute_rh(df_bin['atl03_rh_60'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_50'] = compute_rh(df_bin['atl03_rh_50'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_40'] = compute_rh(df_bin['atl03_rh_40'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_30'] = compute_rh(df_bin['atl03_rh_30'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_25'] = compute_rh(df_bin['atl03_rh_25'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_20'] = compute_rh(df_bin['atl03_rh_20'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_rh_10'] = compute_rh(df_bin['atl03_rh_10'],
                                        df_bin['atl03_ground_median'])
    df_bin['atl03_trad_cc'] = (df_bin['atl03_n_canopy'] +\
        df_bin['atl03_n_hi_canopy']) / (df_bin['atl03_n_canopy'] +\
                        df_bin['atl03_n_hi_canopy'] + df_bin['atl03_n_ground'])

    print('Part7')
    print(time.time() - dt)
    dt = time.time()
                                        
    # df_bin.drop(columns = ['beg_id', 'mid_id', 'end_id', 'res','bin_id',
    #                        'time','easting','northing','crosstrack',
    #                        'alongtrack'], axis = 1, inplace = True)
    
    df_bin.drop(columns = ['beg_id', 'mid_id', 'end_id'],
                axis = 1, inplace = True)

    print('Part8')
    print(time.time() - dt)
    dt = time.time()
    
    return df_bin

def create_truth_bin(df_truth, atl08, res_at = 30):
    import time
    #Create Key DF based on ATL08
    dt = time.time()
    df_truth = orient_df(df_truth, 'alongtrack')
    atl08_df = orient_df(atl08.df, 'alongtrack')
    
    key_df = create_key_df(atl08_df, 'alongtrack', res_at)

    domain_list = [
               'delta_time','delta_time_beg','delta_time_end','latitude',
               'longitude','segment_id_beg','segment_id_end','crosstrack']
    
    for domain in domain_list:
        domain_interp = interpolate_domain(atl08_df.alongtrack, 
                                           atl08_df[domain], key_df.bin_id)
        key_df = pd.concat([key_df,pd.DataFrame(domain_interp,
                                                columns=[domain])],axis=1)

    key_id = key_df.bin_id
    atl08_id = atl08_df.alongtrack

    print(dt - time.time())
    dt = time.time()

    
    minInd = indexMatch(atl08_id,key_id)
    res = res_at
    datalen = len(atl08_id)
    minInd[minInd >= datalen] = (datalen - 1)
    include = np.zeros(len(key_id))
    seg_at_comp = np.array([atl08_id[x] for x in minInd])
    diff = np.abs(seg_at_comp - key_id)
    
    include[diff <= (res)] = 1  
    
    atl08_df = atl08_df.loc[minInd]
    atl08_df = atl08_df.drop(columns = domain_list)

    print(dt - time.time())
    dt = time.time()

    df_bin_08 = pd.concat([key_df.reset_index(drop=True),
                     atl08_df.reset_index(drop=True)], axis=1)
    
    
    df_bin_08 = df_bin_08[include == 1]
    
    #Compute the ATL03 
    
    print(dt - time.time())
    dt = time.time()
    print('part4')
    
    agg_list = ['truth;truth_ground_median;z;median;[2]',
            'truth;truth_canopy_max98;z;get_max98;[4]',
            'atl03;gedi_rh_;z;above_percentile;[2,4]',
            'atl03;truth_rh_;z;above_percentile;[4]',
            'truth;truth_n_ground;z;get_len;[2]',
            'truth;truth_n_canopy;z;get_len;[4]']
    
    target_key, include = get_target_keys(key_df, df_truth, 'alongtrack')
    df_truth = df_truth.reset_index(drop = True)
    df_truth = pd.concat([df_truth,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)
    df_truth = pd.concat([df_truth,pd.DataFrame(include,columns=['include'])],axis=1)
    df_truth = df_truth[df_truth.include == 1]
    
    print(dt - time.time())
    dt = time.time()
    print('part5')
        
    df_bin_truth = agg_keys(key_df, df_truth, agg_list, key = 'bin_id')
    
    df_truth = pd.concat([df_truth,pd.DataFrame(target_key,columns=['bin_id'])],axis=1)

    print(dt - time.time())
    dt = time.time()
    print('part6')
    
    df_bin_truth.drop(df_bin_truth.columns.difference(['bin_id','truth_ground_median',
                                         'truth_canopy_max98',
                                          'gedi_rh_100','gedi_rh_98',
                                          'gedi_rh_90','gedi_rh_80',
                                          'gedi_rh_75','gedi_rh_70',
                                          'gedi_rh_60','gedi_rh_50',
                                          'gedi_rh_40','gedi_rh_30',
                                          'gedi_rh_25','gedi_rh_20',
                                          'gedi_rh_10','truth_rh_100',
                                          'truth_rh_98','truth_rh_90',
                                          'truth_rh_80','truth_rh_75',
                                          'truth_rh_70','truth_rh_60',
                                          'truth_rh_50','truth_rh_40',
                                          'truth_rh_30','truth_rh_25',
                                          'truth_rh_20','truth_rh_10',
                                         'truth_n_ground','truth_n_canopy'
                                         ]),
                   1,inplace=True)

    df_bin_truth['gedi_rh_100'] = compute_rh(df_bin_truth['gedi_rh_100'], 
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_98'] = compute_rh(df_bin_truth['gedi_rh_98'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_90'] = compute_rh(df_bin_truth['gedi_rh_90'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_80'] = compute_rh(df_bin_truth['gedi_rh_80'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_75'] = compute_rh(df_bin_truth['gedi_rh_75'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_70'] = compute_rh(df_bin_truth['gedi_rh_70'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_60'] = compute_rh(df_bin_truth['gedi_rh_60'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_50'] = compute_rh(df_bin_truth['gedi_rh_50'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_40'] = compute_rh(df_bin_truth['gedi_rh_40'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_30'] = compute_rh(df_bin_truth['gedi_rh_30'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_25'] = compute_rh(df_bin_truth['gedi_rh_25'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_20'] = compute_rh(df_bin_truth['gedi_rh_20'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['gedi_rh_10'] = compute_rh(df_bin_truth['gedi_rh_10'],
                                      df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_100'] = compute_rh(df_bin_truth['truth_rh_100'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_98'] = compute_rh(df_bin_truth['truth_rh_98'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_90'] = compute_rh(df_bin_truth['truth_rh_90'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_80'] = compute_rh(df_bin_truth['truth_rh_80'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_75'] = compute_rh(df_bin_truth['truth_rh_75'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_70'] = compute_rh(df_bin_truth['truth_rh_70'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_60'] = compute_rh(df_bin_truth['truth_rh_60'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_50'] = compute_rh(df_bin_truth['truth_rh_50'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_40'] = compute_rh(df_bin_truth['truth_rh_40'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_30'] = compute_rh(df_bin_truth['truth_rh_30'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_25'] = compute_rh(df_bin_truth['truth_rh_25'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_20'] = compute_rh(df_bin_truth['truth_rh_20'],
                                        df_bin_truth['truth_ground_median'])
    df_bin_truth['truth_rh_10'] = compute_rh(df_bin_truth['truth_rh_10'],
                                        df_bin_truth['truth_ground_median'])
        
    df_bin_truth['truth_trad_cc'] = df_bin_truth['truth_n_canopy'] /\
        (df_bin_truth['truth_n_canopy'] + df_bin_truth['truth_n_ground'])


    print(dt - time.time())
    dt = time.time()
    print('part6')
    
    df_bin_truth = pd.merge(df_bin_08, df_bin_truth, on="bin_id",how='left')       
    
    return df_bin_truth

'''
Function to read copy ATL08 and append 30m segment information
'''    
def copy_and_append_atl08(atl03filepath, atl08filepath, out_dir, in_res=30):
    atl08file = atl08filepath.split('/')[-1]
    newfilename = 'ATL08_30m_' + atl08file.split('ATL08_')[1]
    newatl08file = os.path.join(out_dir, newfilename)
    
    copyfile(atl08filepath,newatl08file)
    h5f = h5py.File(newatl08file,'a')
    gt_list = ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']
    for gt in gt_list:
        base_key = gt + '/land_segments/30m_segment/'
        print('ATL03 Heights')
        atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath)
        
        print('ATL08 Land Segments')
        atl08 = get_atl08_struct(atl08filepath, gt, atl03)
        
        print('Match ATL08 to ATL03 by segment')
        upsampled_atl08_bin = create_atl08_bin(atl03, atl08, res_at = in_res)
        
        print('Append Data to ATL08')
        all_cols = upsampled_atl08_bin.columns
        
        for i in range(0,len(all_cols)):
            a = np.array(upsampled_atl08_bin[all_cols[i]])
            h5f[base_key + all_cols[i]] = a
    
    h5f.close()
    
def main(args):
    print(args.atl08)
    atl08file = args.atl08[0].split('/')[-1]
    newfilename = 'ATL08_30m_' + atl08file.split('ATL08_')[1]
    newatl08file = os.path.join(args.outdir[0], newfilename)
    
    copyfile(args.atl08[0],newatl08file)
    h5f = h5py.File(newatl08file,'a')
    gt_list = ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']
    for gt in gt_list:
        base_key = gt + '/land_segments/30m_segment/'
        print('ATL03 Heights')
        atl03 = get_atl03_struct(args.atl03[0], gt, args.atl08[0])
        
        print('ATL08 Land Segments')
        atl08 = get_atl08_struct(args.atl08[0], gt, atl03)
        
        print('Match ATL08 to ATL03 by segment')
        upsampled_atl08_bin = create_atl08_bin(atl03, atl08, res_at = 30)
        
        print('Append Data to ATL08')
        all_cols = upsampled_atl08_bin.columns
        
        for i in range(0,len(all_cols)):
            a = np.array(upsampled_atl08_bin[all_cols[i]])
            h5f[base_key + all_cols[i]] = a
    
    h5f.close()

    
if __name__ == "__main__":
    """ Command line entry point """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("atl03", metavar="atl03",
                         nargs="+", help="Input ATL03 File")
    
    # Required positional argument
    parser.add_argument("atl08", metavar="atl08",
                         nargs="+", help="Input ATL08 File")
    
    # Required positional argument
    parser.add_argument("outdir", metavar="outdir",
                         nargs="+", help="Output ATL08 Location")

    args = parser.parse_args()
    main(args)