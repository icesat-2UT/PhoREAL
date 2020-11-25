# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 19:21:30 2020

@author: eguenther
"""
import pandas as pd
import os
from icesatUtils import getNameParts
import numpy as np
import matplotlib.pyplot as plt
import pickle as pl


# List Files
def create_master(pickle_folder, pickle_list):
    # Load Pkl
    idx = 0
    for file_name in pickle_list:
        year, month, day = identify_date(file_name)
        file = pickle_folder + "/" + file_name
        gtstr = file_name.split('_')[-3]
        spot, strongbeam, flying, _, _, _ = identify_spot_beam_str(gtstr, 
                                                                   year, month, 
                                                                   day)
        
        try:
            print(file)
            
            if idx == 0:
                master = pd.read_pickle(file)
                master.dropna(subset=['truth_n_ground'],inplace=True)
                master['year'] = year
                master['month'] = month
                master['day'] = day
                master['gt'] = gtstr
                master['strongbeam'] = strongbeam
                master['spot'] = spot
                master['flying'] = flying
                idx = idx + 1
            else:
                df = pd.read_pickle(file)
                df.dropna(subset=['truth_n_ground'],inplace=True)
                df['year'] = year
                df['month'] = month
                df['day'] = day
                df['gt'] = gtstr
                df['strongbeam'] = strongbeam
                df['spot'] = spot
                df['flying'] = flying
                master = pd.concat([master, df])
                idx = idx + 1
        except:
            print('Does not exist')
    return master

def identify_strong_beam(atlname):
    atlname = atlname.split('.pkl')[0]
    date = atlname.split('_')[1]
    gtNum = atlname.split('_')[5]
    year = int(date[0:4])
    month = int(date[4:6])
    day = int(date[6:8])
    # ID Strong beam
    strongbeam = 1
    if (((year == 2019) and (month > 9)) or
        ((year == 2019) and (month == 9) and (day >= 9)) or
        ((year == 2018) and (month == 12) and (day <= 28)) or
        ((year == 2018 and month < 12))):
        if gtNum in ['gt1l','gt2l','gt3l']:
            strongbeam = 0
        else:
            strongbeam = 1
    else:
        if gtNum in ['gt1r','gt2r','gt3r']:
            strongbeam = 0        
    return strongbeam

def identify_strong_beam_str(gtNum, year, month, day):

    # ID Strong beam
    strongbeam = 1
    if (((year == 2019) and (month > 9)) or
        ((year == 2019) and (month == 9) and (day >= 9)) or
        ((year == 2018) and (month == 12) and (day <= 28)) or
        ((year == 2018 and month < 12))):
        if gtNum in ['gt1l','gt2l','gt3l']:
            strongbeam = 0
        else:
            strongbeam = 1
    else:
        if gtNum in ['gt1r','gt2r','gt3r']:
            strongbeam = 0        
    return strongbeam

def identify_spot_beam(atlname):
    forward = {
        'gt3r' : 1,
        'gt3l' : 2,
        'gt2r' : 3,
        'gt2l' : 4,
        'gt1r' : 5,
        'gt1l' : 6
        }
    
    backward = {
        'gt3l' : 1,
        'gt3r' : 2,
        'gt2l' : 3,
        'gt2r' : 4,
        'gt1l' : 5,
        'gt1r' : 6
        }

    atlname = atlname.split('.pkl')[0]
    gtNum = atlname.split('_')[5]
    year, month, day = identify_date(atlname)
    # ID Strong beam
    
    if (((year == 2019) and (month > 9)) or
        ((year == 2019) and (month == 9) and (day >= 9)) or
        ((year == 2018) and (month == 12) and (day <= 28)) or
        ((year == 2018 and month < 12))):
        flying = 'forward'
        spot = forward[gtNum]
        if gtNum[-1] == 'r':
            strongbeam = 1
        else:
            strongbeam = 0
        
    else:
        flying = 'backward'
        spot = backward[gtNum]
        if gtNum[-1] == 'l':
            strongbeam = 1
        else:
            strongbeam = 0       
    return spot, strongbeam, flying, year, month, day

def identify_spot_beam_str(gtNum, year, month, day):
    forward = {
        'gt3r' : 1,
        'gt3l' : 2,
        'gt2r' : 3,
        'gt2l' : 4,
        'gt1r' : 5,
        'gt1l' : 6
        }
    
    backward = {
        'gt3l' : 1,
        'gt3r' : 2,
        'gt2l' : 3,
        'gt2r' : 4,
        'gt1l' : 5,
        'gt1r' : 6
        }

    # ID Strong beam
    
    if (((year == 2019) and (month > 9)) or
        ((year == 2019) and (month == 9) and (day >= 9)) or
        ((year == 2018) and (month == 12) and (day <= 28)) or
        ((year == 2018 and month < 12))):
        flying = 'forward'
        spot = forward[gtNum]
        if gtNum[-1] == 'r':
            strongbeam = 1
        else:
            strongbeam = 0
        
    else:
        flying = 'backward'
        spot = backward[gtNum]
        if gtNum[-1] == 'l':
            strongbeam = 1
        else:
            strongbeam = 0       
    return spot, strongbeam, flying, year, month, day

def identify_date(atlname):
    atlname = atlname.split('.pkl')[0]    
    date = atlname.split('_')[1]
    year = int(date[0:4])
    month = int(date[4:6])
    day = int(date[6:8])       
    return year, month, day

def date_list(atllist):
    yearlist = []
    monthlist = []
    daylist = []
    for atlname in atllist:
        year, month, day = identify_date(atlname)
        yearlist.append(year)
        monthlist.append(month)
        daylist.append(day)
    return yearlist, monthlist, daylist
        

def filter_strong_beam(pickle_folder, target):
    pickle_list = os.listdir(pickle_folder)
    sb_list = [identify_strong_beam(item) for item in pickle_list]
    outlist = []
    for i in range(0,len(sb_list)):
        if sb_list[i] == target:
            outlist.append(pickle_list[i])
    return outlist

def filter_spot_beam(pickle_folder, spot, flying = ''):
    pickle_list = os.listdir(pickle_folder)
    item_list =\
        [identify_spot_beam(item) for item in pickle_list]
    spot_list = []
    flying_list = []
    for item in item_list:
        spot_list.append(item[0])
        flying_list.append(item[2])        
        
    outlist = []
    for i in range(0,len(spot_list)):
        if flying == '':
            if spot_list[i] == spot:
                outlist.append(pickle_list[i])            
        else:
            if spot_list[i] == spot and flying_list[i] == flying:
                outlist.append(pickle_list[i])
    return outlist

def compute_error_metrics_sigma(measured, truth, sigma):
    measured[measured > 9000] = np.nan
    nan1 = np.isnan(measured)
    measured = measured[~nan1]
    truth = truth[~nan1]
    nan2 = np.isnan(truth)    
    measured = measured[~nan2]
    truth = truth[~nan2]
    
    res = truth - measured
    
    if sigma > 0:
        mean = np.nanmean(res)
        std = np.nanstd(res)
        mask = ~np.isnan(res)
        sigmathresh = mean + (std * sigma)
        mask[mask] &= res[mask] > sigmathresh
        res[mask] = np.nan
        truth[mask] = np.nan
        mask = ~np.isnan(res)
        sigmathresh = mean - (std * sigma)
        mask[mask] &= res[mask] < sigmathresh
        res[mask] = np.nan   
        truth[mask] = np.nan
        
    me = np.nanmean(res)
    mae = np.nanmean(np.abs(res))
    rmse = np.sqrt(np.nanmean(np.square(res)))
    n = np.count_nonzero(~np.isnan(res))
    
    return me, mae, rmse, n

def compute_error_metrics(measured, truth, cutoff):
    measured[measured > 9000] = np.nan
    nan1 = np.isnan(measured)
    measured = measured[~nan1]
    truth = truth[~nan1]
    nan2 = np.isnan(truth)    
    measured = measured[~nan2]
    truth = truth[~nan2]
    
    res = truth - measured
    
    if cutoff[0] == 0 and cutoff[1] == 0:
        print('No Threshold')
    else:
        mask = ~np.isnan(res)
        sigmathresh = cutoff[1]
        mask[mask] &= res[mask] > sigmathresh
        res[mask] = np.nan
        truth[mask] = np.nan
        mask = ~np.isnan(res)
        sigmathresh = cutoff[0]
        mask[mask] &= res[mask] < sigmathresh
        res[mask] = np.nan   
        truth[mask] = np.nan
        
    me = np.nanmean(res)
    mae = np.nanmean(np.abs(res))
    rmse = np.sqrt(np.nanmean(np.square(res)))
    n = np.count_nonzero(~np.isnan(res))
    
    return me, mae, rmse, n

def compute_percent_error_sigma(measured, truth, sigma):
    measured[measured > 9000] = np.nan
    nan1 = np.isnan(measured)
    measured = measured[~nan1]
    truth = truth[~nan1]
    nan2 = np.isnan(truth)    
    measured = measured[~nan2]
    truth = truth[~nan2]
    
    res = truth - measured
    
    if sigma > 0:
        mean = np.nanmean(res)
        std = np.nanstd(res)
        mask = ~np.isnan(res)
        sigmathresh = mean + (std * sigma)
        mask[mask] &= res[mask] > sigmathresh
        res[mask] = np.nan
        truth[mask] = np.nan
        mask = ~np.isnan(res)
        sigmathresh = mean - (std * sigma)
        mask[mask] &= res[mask] < sigmathresh
        res[mask] = np.nan   
        truth[mask] = np.nan

    me = np.nanmean(res)

    rmse = np.sqrt(np.nanmean(np.square(res)))        
    mep = (me / np.nanmean(truth)) * 100
    rmsp = (rmse / np.nanmean(truth)) * 100

    return mep, rmsp

def compute_percent_error(measured, truth, cutoff):
    measured[measured > 9000] = np.nan
    nan1 = np.isnan(measured)
    measured = measured[~nan1]
    truth = truth[~nan1]
    nan2 = np.isnan(truth)    
    measured = measured[~nan2]
    truth = truth[~nan2]
    
    res = truth - measured
    
    if (cutoff[0] == 0) and (cutoff[1] == 0):
        pass
    else:
        mean = np.nanmean(res)
        std = np.nanstd(res)
        mask = ~np.isnan(res)
        sigmathresh = cutoff[1]
        mask[mask] &= res[mask] > sigmathresh
        res[mask] = np.nan
        truth[mask] = np.nan
        mask = ~np.isnan(res)
        sigmathresh = cutoff[0]
        mask[mask] &= res[mask] < sigmathresh
        res[mask] = np.nan   
        truth[mask] = np.nan

    me = np.nanmean(res)

    rmse = np.sqrt(np.nanmean(np.square(res)))        
    mep = (me / np.nanmean(truth)) * 100
    rmsp = (rmse / np.nanmean(truth)) * 100

    return mep, rmsp

    
def compute_perc_canopy_sigma(truth_ground, truth_canopy, atl03_canopy, sigma):
    
    truth_canopy = truth_canopy - truth_ground
    atl03_canopy = atl03_canopy - truth_ground

    truth_canopy[truth_canopy > 50] = np.nan
    truth_canopy[truth_canopy < 0] = np.nan
    nan1 = np.isnan(truth_canopy)
    truth_canopy = truth_canopy[~nan1]
    atl03_canopy = atl03_canopy[~nan1]   
    atl03_canopy[atl03_canopy > 50] = np.nan
    atl03_canopy[atl03_canopy < 0] = np.nan
    nan1 = np.isnan(atl03_canopy)
    truth_canopy = truth_canopy[~nan1]
    atl03_canopy = atl03_canopy[~nan1]
    
    canopy_res = truth_canopy - atl03_canopy
    
    if sigma > 0:
        mean = np.nanmean(canopy_res)
        std = np.nanstd(canopy_res)
        mask = ~np.isnan(canopy_res)
        sigmathresh = mean + (std * sigma)
        mask[mask] &= canopy_res[mask] > sigmathresh
        canopy_res[mask] = np.nan
        truth_canopy[mask] = np.nan
        atl03_canopy[mask] = np.nan
        mask = ~np.isnan(canopy_res)
        sigmathresh = mean - (std * sigma)
        mask[mask] &= canopy_res[mask] < sigmathresh
        canopy_res[mask] = np.nan   
        truth_canopy[mask] = np.nan
        atl03_canopy[mask] = np.nan 
        
    nan1 = np.isnan(truth_canopy)
    truth_canopy = truth_canopy[~nan1]
    atl03_canopy = atl03_canopy[~nan1]  

    nan1 = np.isnan(atl03_canopy)
    truth_canopy = truth_canopy[~nan1]
    atl03_canopy = atl03_canopy[~nan1]
    
    atl03_canopy_median = np.median(atl03_canopy)
    truth_canopy_median = np.median(truth_canopy)
    
    canopy_perc_dif = (1 - (atl03_canopy_median / truth_canopy_median)) * 100
    n_canopy = len(truth_canopy)

    me = np.nanmean(canopy_res)
    mae = np.nanmean(np.abs(canopy_res))
    rmse = np.sqrt(np.nanmean(np.square(canopy_res)))    
    mep = (me / np.nanmean(truth_canopy)) * 100
    maep = (mae / np.nanmean(truth_canopy)) * 100
    rmsp = (rmse / np.nanmean(truth_canopy)) * 100
    
    return canopy_perc_dif, n_canopy, me, mae, rmse, mep, maep, rmsp

def compute_perc_canopy(truth_ground, truth_canopy, atl03_canopy, cutoff):
    
    truth_canopy = truth_canopy - truth_ground
    atl03_canopy = atl03_canopy - truth_ground

    truth_canopy[truth_canopy > 50] = np.nan
    truth_canopy[truth_canopy < 0] = np.nan
    nan1 = np.isnan(truth_canopy)
    truth_canopy = truth_canopy[~nan1]
    atl03_canopy = atl03_canopy[~nan1]   
    atl03_canopy[atl03_canopy > 50] = np.nan
    atl03_canopy[atl03_canopy < 0] = np.nan
    nan1 = np.isnan(atl03_canopy)
    truth_canopy = truth_canopy[~nan1]
    atl03_canopy = atl03_canopy[~nan1]
    
    canopy_res = truth_canopy - atl03_canopy
  
    if (cutoff[0] == 0) and (cutoff[1] == 0):
        pass
    else:
        mask = ~np.isnan(canopy_res)
        sigmathresh = cutoff[1]
        mask[mask] &= canopy_res[mask] > sigmathresh
        canopy_res[mask] = np.nan
        truth_canopy[mask] = np.nan
        atl03_canopy[mask] = np.nan
        mask = ~np.isnan(canopy_res)
        sigmathresh = cutoff[0]
        mask[mask] &= canopy_res[mask] < sigmathresh
        canopy_res[mask] = np.nan   
        truth_canopy[mask] = np.nan
        atl03_canopy[mask] = np.nan 
        
    nan1 = np.isnan(truth_canopy)
    truth_canopy = truth_canopy[~nan1]
    atl03_canopy = atl03_canopy[~nan1]  

    nan1 = np.isnan(atl03_canopy)
    truth_canopy = truth_canopy[~nan1]
    atl03_canopy = atl03_canopy[~nan1]
    
    atl03_canopy_median = np.median(atl03_canopy)
    truth_canopy_median = np.median(truth_canopy)
    
    canopy_perc_dif = (1 - (atl03_canopy_median / truth_canopy_median)) * 100
    n_canopy = len(truth_canopy)

    me = np.nanmean(canopy_res)
    mae = np.nanmean(np.abs(canopy_res))
    rmse = np.sqrt(np.nanmean(np.square(canopy_res)))    
    mep = (me / np.nanmean(truth_canopy)) * 100
    maep = (mae / np.nanmean(truth_canopy)) * 100
    rmsp = (rmse / np.nanmean(truth_canopy)) * 100
    
    return canopy_perc_dif, n_canopy, me, mae, rmse, mep, maep, rmsp

def find_cutoffs(res, sigma_list = [0,1,2,3]):
    cut_off_list = []
    for sigma in sigma_list:
        if sigma == 0:
            cut_off = [0,0]
        else:
            mean = np.nanmean(res)
            std = np.nanstd(res)
            sigma_0 = mean - (std * sigma)        
            sigma_1 = mean + (std * sigma)   
            cut_off = [sigma_0, sigma_1]
        cut_off_list.append(cut_off)
    return cut_off_list

def noise_rate(df_p):
    
    tot_shots = df_p.gc_nshots_unclassed + df_p.gc_nshots_ground +\
        df_p.gc_nshots_canopy
    # gc_noise = (df.gc_nshots_unclassed + df.gc_nshots_ground +\
    #     df.gc_nshots_canopy) - df.gc_n_ground - df.gc_n_canopy
        
    # pc_noise = (df.gc_nshots_unclassed + df.gc_nshots_ground +\
    #     df.gc_nshots_canopy) - df.pc_n_ground - df.pc_n_canopy
        
    pr_ground = np.nanmean(df_p.gc_n_ground) / np.nanmean(tot_shots)
    pr_canopy = np.nanmean(df_p.gc_n_canopy) / np.nanmean(tot_shots)    
    pr_noise = np.nanmean(df_p.gc_n_unclassified) / np.nanmean(tot_shots)
    
    pr_pc_ground = np.nanmean(df_p.pc_n_ground) / np.nanmean(tot_shots)
    pr_pc_canopy = np.nanmean(df_p.pc_n_canopy) / np.nanmean(tot_shots)    
    pr_pc_noise = np.nanmean(df_p.gc_n_unclassified) / np.nanmean(tot_shots)
    
    return pr_ground, pr_canopy, pr_noise, pr_pc_ground, pr_pc_canopy, pr_pc_noise

def noise_rate_canopy_filter(df_in, cutoff):
    df_p = copy.deepcopy(df_in)
    df_p.atl03_canopy_height = df_p.atl03_canopy_max98 - df_p.truth_ground_median
    df_p.truth_canopy_height = df_p.truth_canopy_max98 - df_p.truth_ground_median
    df_p['res_canopy_height'] = df_p.truth_canopy_height - df_p.atl03_canopy_height
    print(len(df))
    if ((cutoff[0] == 0) and (cutoff[1] == 0)):
        print('ok')
    else:
        df_p.res_canopy_height[df_p.res_canopy_height < cutoff[0]] = np.nan
        df_p.res_canopy_height[df_p.res_canopy_height > cutoff[1]] = np.nan
        df_p.dropna(subset=['res_canopy_height'], inplace=True)
    df_p.dropna(subset=['res_canopy_height'], inplace=True)
    df_p.dropna(subset=['gc_nshots_ground'], inplace=True)
    df_p.dropna(subset=['gc_nshots_canopy'], inplace=True)
    df_p.dropna(subset=['gc_nshots_unclassed'], inplace=True)
    pr_ground, pr_canopy, pr_noise, pr_pc_ground, pr_pc_canopy, pr_pc_noise =\
        noise_rate(df_p)
    return pr_ground, pr_canopy, pr_noise, pr_pc_ground, pr_pc_canopy, pr_pc_noise
    
def work(df, ct_atl03_gm, ct_atl08_gm, ct_atl08_gbf, ct_atl03_cm, 
         ct_atl03_chm):
    # atl03_ground_median
    # h_te_median
    # h_te_best_fit
    # pc_ground_median
    # atl03_canopy_max98
    # h_max_canopy_abs
    # pc_canopy_max98
    atl03_ground = df.atl03_ground_median
    truth_ground = df.truth_ground_median
    truth_canopy = df.truth_canopy_max98
    atl03_canopy = df.atl03_canopy_max98
    atl08_median = df.h_te_median
    atl08_best_fit = df.h_te_best_fit
    pc_ground = df.pc_ground_median
    pc_canopy = df.pc_canopy_max98
    gme, gmae, grmse, gn = compute_error_metrics(atl03_ground, 
                                                              truth_ground,
                                                 ct_atl03_gm)
    gmem, gmaem, grmsem, gnm = compute_error_metrics(atl08_median, 
                                                     truth_ground, ct_atl08_gm)
    gmebf, gmaebf, grmsebf, gnbf = compute_error_metrics(atl08_best_fit,
                                                         truth_ground, ct_atl08_gbf)
    gmepc, gmaepc, grmsepc, gnpc = compute_error_metrics(pc_ground,
                                                         truth_ground, ct_atl03_gm)
    cme, cmae, crmse, cn = compute_error_metrics(atl03_canopy, truth_canopy,
                                                 ct_atl03_cm)
    cmepc, cmaepc, crmsepc, cnpc = compute_error_metrics(pc_canopy, 
                                                         truth_canopy, ct_atl03_cm)
    canopy_perc_dif, n_canopy, reme, remae, rermse, remep, remaep, rermsp\
        = compute_perc_canopy(truth_ground, truth_canopy, atl03_canopy, ct_atl03_chm)

    pccan_perc_dif, n_pccan, pccan_reme, pccan_remae, pccan_rermse, pccan_remep,\
    pccan_remaep, pccan_rermsp\
        = compute_perc_canopy(truth_ground, truth_canopy, pc_canopy, ct_atl03_chm)
    
    gmep, grmsp = compute_percent_error(atl03_ground, truth_ground, ct_atl03_gm)
    gmepm, grmspm = compute_percent_error(atl08_median, 
                                                     truth_ground, ct_atl08_gm)
    gmepbf, grmspbf  = compute_percent_error(atl08_best_fit,
                                                         truth_ground, ct_atl08_gbf)
    gmeppc, grmsppc = compute_percent_error(pc_ground,truth_ground, ct_atl03_gm)
    cmep, crmsp = compute_percent_error(atl03_canopy, truth_canopy,
                                                 ct_atl03_cm)
    cmeppc, crmsppc = compute_percent_error(pc_canopy, truth_canopy,
                                                 ct_atl03_cm)    
    pr_ground, pr_canopy, pr_noise, pr_pc_ground, pr_pc_canopy, pr_pc_noise =\
        noise_rate_canopy_filter(df, ct_atl03_chm)
        
    out_string = str(gme) + "," + str(gmae) + "," + str(grmse) + "," +\
        str(gn) + "," + str(gmep) + "," + str(grmsp) + ","
    out_string = out_string + str(gmem) + "," + str(gmaem) + "," +\
        str(grmsem) + "," +\
        str(gnm) + "," + str(gmepm) + "," + str(grmspm) + ","
    out_string = out_string + str(gmebf) + "," + str(gmaebf) + "," +\
        str(grmsebf) + "," +\
        str(gnbf) + "," + str(gmepbf) + "," + str(grmspbf) + ","
    out_string = out_string + str(gmepc) + "," + str(gmaepc) + "," +\
        str(grmsepc) + "," +\
        str(gnpc) + "," + str(gmeppc) + "," + str(grmsppc) + ","
    out_string = out_string + str(cme) + "," + str(cmae) + "," + str(crmse) +\
        "," + str(cn) + ","  + str(cmep) + "," + str(crmsp) + ","
    out_string = out_string + str(cmepc) + "," + str(cmaepc) + "," +\
        str(crmsepc) + "," + str(cnpc) + "," + str(cmeppc) + "," +\
            str(crmsppc) + ","
    out_string = out_string + str(canopy_perc_dif) + "," + str(n_canopy) +\
        "," + str(reme) + "," + str(remae) + "," + str(rermse) + "," +\
            str(remep) + "," +  str(remaep) + "," + str(rermsp) + "," +\
                str(pccan_remep) + "," + str(pccan_remaep) + "," + str(pccan_rermsp)
                
    out_string = out_string + "," +\
        str(pr_ground) +\
            "," + str(pr_canopy) +\
                "," + str(pr_noise) + "," + str(pr_pc_ground) +\
                        "," + str(pr_pc_canopy) +\
                            "," + str(pr_pc_noise)
    return out_string

def compute_residual(measured, truth, sigma):
    measured[measured > 9000] = np.nan
    nan1 = np.isnan(measured)
    measured = measured[~nan1]
    truth = truth[~nan1]
    nan2 = np.isnan(truth)    
    measured = measured[~nan2]
    truth = truth[~nan2]
    
    res = truth - measured
    
    if sigma > 0:
        mean = np.nanmean(res)
        std = np.nanstd(res)
        mask = ~np.isnan(res)
        sigmathresh = mean + (std * sigma)
        mask[mask] &= res[mask] > sigmathresh
        res[mask] = np.nan
        truth[mask] = np.nan
        mask = ~np.isnan(res)
        sigmathresh = mean - (std * sigma)
        mask[mask] &= res[mask] < sigmathresh
        res[mask] = np.nan   
        truth[mask] = np.nan
        
    return res

def pickle_residuals(df, sigma, out_file_ground, out_file_canopy):
    atl03_ground = df.atl03_ground_median
    truth_ground = df.truth_ground_median
    truth_canopy = df.truth_canopy_max98
    atl03_canopy = df.atl03_canopy_max98
    
    ground_res = compute_residual(atl03_ground, truth_ground, sigma)
    canopy_res = compute_residual(atl03_canopy, truth_canopy, sigma)
    
    ground_res_df = pd.DataFrame(ground_res)
    canopy_res_df = pd.DataFrame(canopy_res)
    ground_res_df.to_pickle(out_file_ground)
    canopy_res_df.to_pickle(out_file_canopy)
    
    return ground_res, canopy_res

def compute_relative_canopy(df, min_val = 0, max_val = 50):
    df_out = copy.deepcopy(df)
    df_out['atl03_canopy_height'] = df_out.atl03_canopy_max98 -\
        df.truth_ground_median
    df_out['truth_canopy_height'] = df_out.truth_canopy_max98 -\
        df.truth_ground_median
    df_out['pc_canopy_height'] = df_out.pc_canopy_max98 -\
        df.truth_ground_median

    # df_out['atl03_canopy_height'][df_out['atl03_canopy_height'] > max_val] =\
    #     np.nan 
    df_out['truth_canopy_height'][df_out['truth_canopy_height'] > max_val] =\
        np.nan 
    # df_out['pc_canopy_height'][df_out['pc_canopy_height'] > max_val] = np.nan 
    # df_out['atl03_canopy_height'][df_out['atl03_canopy_height'] < min_val] =\
    #     np.nan
    df_out['truth_canopy_height'][df_out['truth_canopy_height'] < min_val] =\
        np.nan 
    # df_out['pc_canopy_height'][df_out['pc_canopy_height'] < min_val] = np.nan 
    
    # df_out.dropna(subset=['atl03_canopy_height'],inplace=True)
    df_out.dropna(subset=['truth_canopy_height'],inplace=True)
    # df_out.dropna(subset=['pc_canopy_height'],inplace=True)

    
    return df_out

def graph_errors(atl08_folder, atl08_file, truth_folder, atl03_folder,
                 output_folder):
    in_file = atl08_folder + '/' + atl08_file
    df = pd.read_pickle(in_file)
    

    
    truth_file = 'truth_' + atl08_file.split('ATL08_')[1]
    truth_file = truth_folder + '/' + truth_file
    atl03_file = 'ATL03_' + atl08_file.split('ATL08_')[1]
    atl03_file = atl03_folder + '/' + atl03_file
    
    df_truth = pd.read_pickle(truth_file)
    df_atl03 = pd.read_pickle(atl03_file)
     
    tgy = df_truth.alongtrack[df_truth.classification == 2]
    tgz = df_truth.z[df_truth.classification == 2]
    tcy = df_truth.alongtrack[df_truth.classification == 4]
    tcz = df_truth.z[df_truth.classification == 4]

    agy = df_atl03.alongtrack[df_atl03.classification == 1]
    agz = df_atl03.h_ph[df_atl03.classification == 1]
    acy = df_atl03.alongtrack[df_atl03.classification == 2]
    acz = df_atl03.h_ph[df_atl03.classification == 2]
    ahy = df_atl03.alongtrack[df_atl03.classification == 3]
    ahz = df_atl03.h_ph[df_atl03.classification == 3]
    auy = df_atl03.alongtrack[df_atl03.classification > 1]
    auz = df_atl03.h_ph[df_atl03.classification > 1]
    
    ratio = int(np.ceil(len(tcy)/5000))
    
    f1 = plt.figure()
    title_name = (atl08_file.split('.pkl')[0] + ' Canopy Height Error')
    # f1 = plt.title(title_name)
    ax1 = plt.subplot(2,1,1)
    ax1.set_title(title_name)
    ax1.plot(tcy[::ratio], tcz[::ratio], '.', color = [0.6,0.6,0.6],
             label = 'Truth Canopy')
    ax1.plot(tgy[::ratio], tgz[::ratio], '.', color = [0.3,0.3,0.3],
             label = 'Truth Ground')
    
    ax1.plot(auy, auz, '.', color = [0.8,0.8,1],
             label = 'ATL03 Unclassified')    
    ax1.plot(ahy, ahz, '.', color = [0,0.8,0],
             label = 'ATL03 High Canopy')
    ax1.plot(acy, acz, '.', color = [0,0.5,0],
             label = 'ATL03 Canopy')    
    ax1.plot(agy, agz, '.', color = [0.8,0.5,0],
             label = 'ATL03 Ground')        

    
    ax1.plot(df.alongtrack, df.atl03_canopy_max98,
             'o', color = [0.9,0.1,0.0], label = 'ATL08 Max 98 Canopy Height')    
    ax1.plot(df.alongtrack, df.truth_canopy_max98, 'x', color = [0.0,0.3,0.9],
             label = 'Truth Max 98 Canopy Height')    
    ax1.set_ylabel('Height (m)')
    ax1.set_xlabel('Along Track (m)')    
    ax1.legend()
    
    ax2 = plt.subplot(2,1,2, sharex=ax1)
    ax2.set_title('Canopy Height Residuals')
    
    df = compute_relative_canopy(df)
    res_ch = df.truth_canopy_height - df.atl03_canopy_height
    res_pc = df.truth_canopy_height - df.pc_canopy_height
    
    mean_ch = np.nanmean(res_ch)
    mean_pc = np.nanmean(res_pc)
    
    ax2.plot(df.alongtrack, res_ch,'bo', label='Canopy Height Residuals')
    ax2.plot(df.alongtrack, res_pc,'ro', label='PC Canopy Height Residuals')   

    ax2.plot([np.min(df.alongtrack), np.max(df.alongtrack)], 
             [mean_ch, mean_ch],'b--', label = 'Canopy Mean Error')
    ax2.plot([np.min(df.alongtrack), np.max(df.alongtrack)], 
             [mean_pc, mean_pc],'r--', label = 'PC Canopy Mean Error') 
    ax2.plot([np.min(df.alongtrack), np.max(df.alongtrack)], 
             [0, 0],'k-.', label = 'Best Fit', color = [0.8, 0.8, 0.8]) 
    ax2.set_ylabel('Residuals (m)')
    ax2.set_xlabel('Along Track (m)')
    ax2.legend()
    
    out_pickle = output_folder + '/pkl/' + atl08_file
    out_png = output_folder + '/png/' + atl08_file.split('.pkl')[0] + '.png'
    
    pl.dump(f1,open(str(out_pickle),'wb'))
    f1.savefig(out_png)
    
def rmsep_graphs(measured_in, truth_in, cutoff, min_lvl = 0, max_lvl = 40, freq = 5,
                 xlab ='Truth',ylab ='Measured',
                 tlab='RMSE Plot', outpng = 'out_rmse.png', outpkl = 'out_rmse.pkl'):
    mep_list = []
    maep_list = []
    rmsep_list = []
    point_list = []

    for i in range(min_lvl,max_lvl,freq):
        truth = copy.deepcopy(truth_in)
        measured = copy.deepcopy(measured_in)
        truth[truth < i ] == np.nan
        truth[truth > i + freq] = np.nan
        measured[truth < i] == np.nan
        measured[truth > i + freq] = np.nan
        res = truth_in - measured_in
        if (cutoff[0] == 0) and (cutoff[1] == 0):
            pass
        else:
            mean = np.nanmean(res)
            std = np.nanstd(res)
            mask = ~np.isnan(res)
            sigmathresh = cutoff[1]
            mask[mask] &= res[mask] > sigmathresh
            res[mask] = np.nan
            truth[mask] = np.nan
            mask = ~np.isnan(res)
            sigmathresh = cutoff[0]
            mask[mask] &= res[mask] < sigmathresh
            res[mask] = np.nan   
            truth[mask] = np.nan
        res = truth - measured        
        me = np.nanmean(res)
        mae = np.nanmean(np.abs(res))    
        rmse = np.sqrt(np.nanmean(np.square(res)))        
        mep = (me / np.nanmean(truth)) * 100
        maep = (mae / np.nanmean(truth)) * 100
        rmsp = (rmse / np.nanmean(truth)) * 100    
        
        mep_list.append(mep)
        maep_list.append(maep)
        rmsep_list.append(rmsp)
        point_list.append(i)
    f = plt.figure(figsize=(8,5))
    f = plt.plot(point_list, mep_list,'-o',label = 'Mean Error %')
    f = plt.plot(point_list, maep_list,'-o',label = 'Mean Absolute Error %')
    f = plt.plot(point_list, rmsep_list,'-o',label = 'Root Mean Square Error %')
    f = plt.xlim(0,35)
    f = plt.ylim(-20,100)
    f = plt.legend()
    f = plt.xlabel('Height (m)')
    f = plt.ylabel('Percent Error')
    f = plt.title('Percent Error by Height (All beams)')
    f = plt.savefig('/LIDAR/server/USERS/eric/1_experiment/finland_analysis3/graphs/rmse_by_height/png/all.png')
    
    
def rmse_graphs(measured_in, truth_in, cutoff, min_lvl = 0, max_lvl = 40, freq = 5,
                 xlab ='Truth',ylab ='Measured',
                 tlab='RMSE Plot', outpng = 'out_rmse.png', outpkl = 'out_rmse.pkl'):
    me_list = []
    mae_list = []
    rmse_list = []
    point_list = []

    for i in range(min_lvl,max_lvl,freq):
        truth = copy.deepcopy(truth_in)
        measured = copy.deepcopy(measured_in)
        truth[truth < i ] == np.nan
        truth[truth > i + freq] = np.nan
        measured[truth < i] == np.nan
        measured[truth > i + freq] = np.nan
        res = truth_in - measured_in
        if (cutoff[0] == 0) and (cutoff[1] == 0):
            pass
        else:
            mean = np.nanmean(res)
            std = np.nanstd(res)
            mask = ~np.isnan(res)
            sigmathresh = cutoff[1]
            mask[mask] &= res[mask] > sigmathresh
            res[mask] = np.nan
            truth[mask] = np.nan
            mask = ~np.isnan(res)
            sigmathresh = cutoff[0]
            mask[mask] &= res[mask] < sigmathresh
            res[mask] = np.nan   
            truth[mask] = np.nan
        res = truth - measured        
        me = np.nanmean(res)
        mae = np.nanmean(np.abs(res))    
        rmse = np.sqrt(np.nanmean(np.square(res)))        
        
        me_list.append(me)
        mae_list.append(mae)
        rmse_list.append(rmse)
        point_list.append(i)
    f = plt.figure(figsize=(8,5))
    f = plt.plot(point_list, me_list,'-o',label = 'Mean Error %')
    f = plt.plot(point_list, mae_list,'-o',label = 'Mean Absolute Error %')
    f = plt.plot(point_list, rmse_list,'-o',label = 'Root Mean Square Error %')
    f = plt.xlim(0,35)
    f = plt.ylim(-10,10)
    f = plt.legend()
    f = plt.xlabel('Height (m)')
    f = plt.ylabel('Error (m)')
    f = plt.title('Error by Height (Strong beams, Target Only)')
    f = plt.savefig('/LIDAR/server/USERS/eric/1_experiment/finland_analysis3/graphs/rmse_by_height/png/all.png')
    
if __name__ == "__main__":
    
    from sklearn.metrics import mean_squared_error
    from math import sqrt
    import numpy as np
    import csv
    import copy
    
    folder08 = '/LIDAR/server/USERS/eric/1_experiment/finland_analysis5/' +\
                      'pkl/atl08'
    folder = '/LIDAR/server/USERS/eric/1_experiment/finland_analysis5/' +\
                          'pkl/master'

    csv_filepath = folder + '/' + 'csv_test55c.csv'
    beam_type = 'All'
    strong_list = os.listdir(folder08)                          
    # strong_list = filter_strong_beam(folder08, 0)

    sb = create_master(folder08, strong_list)

    # sb_file = folder + '/' + 'all_dates.pkl'
    # sb_file = '/LIDAR/server/USERS/eric/1_experiment/finland_analysis3/pkl/'+\
    #     'master/all_20200526.pkl'
    # sb = pd.read_pickle(sb_file)
    
    
    res_atl03_gm = sb.truth_ground_median - sb.atl03_ground_median
    res_atl08_gm = sb.truth_ground_median - (sb.h_te_median + sb.zoffset)
    res_atl08_gbf = sb.truth_ground_median - (sb.h_te_median + sb.zoffset)
    res_atl03_cm = sb.truth_canopy_max98 - sb.atl03_canopy_max98
    truth_canopy_height = sb.truth_canopy_max98 - sb.truth_ground_median
    atl03_canopy_height = sb.atl03_canopy_max98 - sb.truth_ground_median
    res_atl03_chm = truth_canopy_height - atl03_canopy_height
    

    res_atl03_gm[res_atl03_gm   < -1000] = np.nan
    res_atl08_gm[res_atl08_gm   < -1000] = np.nan
    res_atl08_gbf[res_atl08_gbf < -1000] = np.nan
    res_atl03_cm[res_atl03_cm   < -1000] = np.nan
    res_atl03_chm[res_atl03_chm < -1000] = np.nan
    res_atl03_gm[res_atl03_gm   > 1000] = np.nan
    res_atl08_gm[res_atl08_gm   > 1000] = np.nan
    res_atl08_gbf[res_atl08_gbf > 1000] = np.nan
    res_atl03_cm[res_atl03_cm   > 1000] = np.nan
    res_atl03_chm[res_atl03_chm > 1000] = np.nan
    ct_atl03_gm = find_cutoffs(res_atl03_gm)
    ct_atl08_gm = find_cutoffs(res_atl08_gm)
    ct_atl08_gbf = find_cutoffs(res_atl08_gbf)
    ct_atl03_cm = find_cutoffs(res_atl03_cm)
    ct_atl03_chm = find_cutoffs(res_atl03_chm)  
    
    header = 'Beam,Sigma,Solar Elevation,Snow Cover,MODIS LC,LC Label,Forest,Season,'+\
        'ATL03 Ground ME,ATL03 Ground MAE,ATL03 Ground RMSE,ATL03 Ground N,' +\
            'ATL03 Ground % Mean Error Bias,' + 'ATL03 Ground % RMSE,' +\
        'ATL08 Ground Median ME,ATL08 Ground Median MAE,' +\
            'ATL08 Ground Median RMSE,ATL08 Ground Median N,' +\
                'ATL8 Ground Median % Mean Error Bias,' +\
                    'ATL08 Ground Median% RMSE,' +\
        'ATL08 Ground Best Fit ME,ATL08 Ground Best Fit MAE,'+\
            'ATL08 Ground Best Fit RMSE,ATL08 Ground Best Fit N,' +\
            'ATL08 Ground Best Fit % Mean Error Bias,' +\
                'ATL08 Ground Best Fit  % RMSE,' +\
        'PC Ground ME,PC Ground MAE,PC Ground RMSE,PC Ground N,' +\
         'PC Ground % Mean Error Bias,' + 'PC Ground % RMSE,' +\
        'ATL03 Canopy ME,ATL03 Canopy MAE,ATL03 Canopy RMSE,ATL03 Canopy N,' +\
            'ATL03 Canopy % Mean Error Bias,' +\
                'ATL03 Canopy  % RMSE,' +\
        'PC Canopy ME,PC Canopy MAE,PC Canopy RMSE,PC Canopy N,' +\
                    'PC Canopy % Mean Error Bias,' +\
                'PC Canopy  % RMSE,' +\
                'Canopy Perc Dif,Canopy Perc Dif N,Canopy Height ME,' +\
                    'Canopy Height MAE,Canopy Height RMSE,' +\
                        'Canopy Height %ME,Canopy Height %MAE,' +\
                            'Canopy Height %RMSE, ' +\
                                'PC Canopy Height %ME,PC Canopy Height %MAE,' +\
                            'PC Canopy Height %RMSE,PR Ground, PR Canopy,PR Noise,'+\
                                'PC PR Ground,PC PR Canopy,PC PR Noise'
    # out_string = work(sb, ct_atl03_gm, ct_atl08_gm, ct_atl08_gbf, ct_atl03_cm,
    #                   ct_atl03_chm)
       # [0, 2, 3, 4, 6, 7, 9, 10, 11, 12,13, 14, 15, 16, 17, 254, 255]
    csvfile = open(csv_filepath,'w')
    csvfile.write(header)
    csvfile.write('\n')
    lc_list = list(np.unique(sb.segment_landcover))
    lc_list = [-1, -2, -3] + lc_list
    lc_dict = {
        -1 : "All",
        -2 : "Target",
        -3 : "Non-target",
        0 : "NA",
        1 : "Evergreen Needleleaf Forest",
        2 : "Evergreen Broadleaf Forest",
        3 : "Deciduous Needleleaf Forest",
        4 : "Deciduous Broadleaf Forest",
        5 : "Mixed Forests",
        6 : "Closed Shrublands",
        7 : "Open Shrublands",
        8 : "Woody Savannas",
        9 : "Savannas",
        10 : "Grasslands",
        11 : "Permanent Wetlands",
        12 : "Croplands",
        13 : "Urban and Built-Up",
        14 : "Cropland-Natural Vegetation Mosaic",
        15 : "Snow and Ice",
        16 : "Barren or Sparsley Vegetated",
        17 : "IGBP Water Bodies*",
        254 : "Unclassified",
        255 : "Fill Value"
        }
  # Corine_LC   
 
    corine_list = list(np.unique(sb.Corine_LC))
    corine_list = [-1] + lc_list
    corine_dict = {
        -1 : "All",
        111 : "Continuous urban fabric",
        112 : "Discontinuous urban fabric",
        121 : "Industrial or commercial units",
        122 : "Road and rail networks and associated land",
        123 : "Port areas",
        124 : "Airports",
        131 : "Mineral extraction sites",
        132 : "Dump sites",
        133 : "Construction sites",
        141 : "Green urban areas",
        142 : "Sport and leisure faclities",
        211 : "Non-irrigated arable land",
        212 : "Croplands",
        213 : "Rice fields",
        221 : "Vineyards",
        222 : "Fruit trees and berry plantations",
        223 : "Olive groves",
        231 : "Pastures",
        241 : "Annual crops associated with permanent crops",
        242 : "Complex cultivation patterns",
        243 : "Land principally occupied by agricultural/vegetation",
        244 : "Agro-forestry areas",
        311 : "Broad-leaved forest",
        312 : "Coniferous forest",
        313 : "Mixed forest",
        321 : "Natural grasslands",
        322 : "Moors and heathland",
        323 : "Sclerophyllous vegetation",
        324 : "Transitional woodland-shrub",
        331 : "Beaches, dunes, sands",
        332 : "Bare rocks",
        333 : "Sparsly vegetated areas",
        334 : "Burnt areas",
        335 : "Glaciers and perpetual snow",
        411 : "Inland marshes",
        412 : "Peat bogs",
        421 : "Salt marshes",
        422 : "Salines",
        423 : "Intertidal flats",
        511 : "Water courses",
        512 : "Water bodies",
        521 : "Coastal lagoons",
        522 : "Estuaries",
        523 : "Sea and oceans"
        }
    
     #For Sigma 0 to 3
    beam_list = ['All','Strong','Weak']
    summer_list = ['All','Summer','Nonsummer']
    corine_forest_list = ['All','Forest','Nonforest']
    # beam_list = ['All']
    # summer_list = ['All']
    # corine_forest_list = ['All']
    for beam in beam_list:          
        for sigma in range(0,4):
            for season in summer_list:
                for forest in corine_forest_list:
                    for solar_elevation in range(0,3):
                        for snow_cover in range(0,3):
                            for lc in lc_list:
                                df = copy.deepcopy(sb)
                                #Filter beam type
                                if beam == 'Strong':
                                    df.strongbeam[df.strongbeam == 0] = np.nan
                                    df.dropna(subset=['strongbeam'],inplace=True) 
                                elif beam == 'Weak':
                                    df.strongbeam[df.strongbeam == 1] = np.nan
                                    df.dropna(subset=['strongbeam'],inplace=True)  
                                # Filter Solar Elevation
                                if solar_elevation == 1:
                                    df.night_flag[df.night_flag == 1] = np.nan
                                    df.dropna(subset=['night_flag'], inplace=True)
                                    solar_label = 'Day'
                                elif solar_elevation == 2:
                                    df.night_flag[df.night_flag == 0] = np.nan
                                    df.dropna(subset=['night_flag'], inplace=True)
                                    solar_label = 'Night'
                                else:
                                    solar_label = 'All' 
                                # Filter Snow Cover
                                if snow_cover == 1:
                                    df.segment_snowcover[df.segment_snowcover == 2] = np.nan
                                    df.segment_snowcover[df.segment_snowcover == 3] = np.nan
                                    df.dropna(subset=['segment_snowcover'], inplace=True)
                                    snow_label = 'No snow'
                                elif snow_cover == 2:
                                    df.segment_snowcover[df.segment_snowcover == 0] = np.nan
                                    df.segment_snowcover[df.segment_snowcover == 1] = np.nan
                                    df.dropna(subset=['segment_snowcover'], inplace=True)
                                    snow_label = 'Snow'    
                                else:
                                    snow_label = 'All'
                                # Filter Summer
                                if season == 'Summer':
                                    df.month[df.month < 5] = np.nan
                                    df.month[df.month > 9] = np.nan
                                    df.dropna(subset=['month'], inplace=True)
                                elif season == 'Nonsummer':
                                    df.month[(df.month >= 5) & (df.month <= 9)] = np.nan
                                    df.dropna(subset=['month'], inplace=True)    
                                # Filter Forest (Corine)
                                if forest == 'Forest':
                                    df.Corine_LC[df.Corine_LC < 310] = np.nan
                                    df.Corine_LC[df.Corine_LC == 321] = np.nan 
                                    df.Corine_LC[df.Corine_LC == 322] = np.nan 
                                    df.Corine_LC[df.Corine_LC == 323] = np.nan 
                                    df.Corine_LC[df.Corine_LC > 325] = np.nan
                                    df.dropna(subset=['Corine_LC'], inplace=True)
                                elif forest == 'Nonforest':
                                    df.Corine_LC[df.Corine_LC == 311] = np.nan
                                    df.Corine_LC[df.Corine_LC == 312] = np.nan 
                                    df.Corine_LC[df.Corine_LC == 313] = np.nan 
                                    df.Corine_LC[df.Corine_LC == 324] = np.nan 
                                    df.dropna(subset=['Corine_LC'], inplace=True)
                                # Filter Landcover
                                if lc == -1:
                                    lc_label = lc_dict[lc]
                                elif lc == -2:
                                    target = [0, 2, 3, 4, 6, 7, 9, 10, 11, 12, 13, 14, 15,
                                              16, 17, 254, 255]
                                    lc_label = lc_dict[lc]
                                    df.segment_landcover[np.isin(df.segment_landcover,
                                                                 target)] = np.nan     
                                    df.dropna(subset=['segment_landcover'],inplace=True)
                                elif lc == -3:
                                    target = [1, 5, 8]
                                    lc_label = lc_dict[lc]
                                    df.segment_landcover[np.isin(df.segment_landcover,
                                                                 target)] = np.nan  
                                    df.dropna(subset=['segment_landcover'],inplace=True)
                                else:
                                    df.segment_landcover[df.segment_landcover != int(lc)]\
                                        = np.nan
                                    df.dropna(subset=['segment_landcover'],inplace=True)
                                    lc_label = lc_dict[lc]
            
                                print(str(sigma) + "," + str(solar_label) + "," +
                                      str(snow_cover) + "," + str(lc))
                                label = (beam + "," + str(sigma) + "," + str(solar_label) + "," +
                                      str(snow_label) + "," + str(lc) + ',' + lc_label +
                                      ',' + forest + ',' + season)
                                
                                out_file_ground = folder + "/res/" + beam_type + "_"  +\
                                    solar_label + "_"  + snow_label + "_" + str(lc_label) +\
                                        "_" + str(sigma) + "_ground.pkl"
                                out_file_canopy = folder + "/res/" + beam_type + "_"  +\
                                    solar_label + "_"  + snow_label + "_" + str(lc_label) +\
                                        "_" + str(sigma) + "_canopy.pkl"
                                ground_res, canopy_res = pickle_residuals(df, sigma, 
                                                                    out_file_ground, 
                                                 out_file_canopy)
                                # result_string = work(df, sigma)
                                result_string = work(df, ct_atl03_gm[sigma], 
                                                     ct_atl08_gm[sigma], 
                                                     ct_atl08_gbf[sigma], 
                                                     ct_atl03_cm[sigma], 
                                                     ct_atl03_chm[sigma])
                                out_string = label + "," + result_string + '\n'
                                csvfile.write(out_string)
    
    
    
    csvfile.close()
    
    # import matplotlib.pyplot as plt
    # my_cmap = matplotlib.cm.get_cmap('plasma')
    # my_cmap.set_under('w')   
    # p = plt.hist2d(truth_canopy, atl03_canopy, bins = 500, cmap = my_cmap, cmin = 5)
    # p = plt.axis('equal')
    # p = plt.axis('tight')
    # p = plt.xlabel('Airborne Lidar Relative Canopy Height (m)')
    # p = plt.ylabel('ATL08 Relative Canopy Height (m)')
    # p = plt.title('Relative Canopy Height, ATL08 vs Airborne Lidar')
    # p = plt.colorbar()
    # x = np.array([-5, 40])
    # y = np.array([-5, 40])
    # p = plt.plot(x,y,'--k')

    
    
    # me, mae, rmse, n = compute_error_metrics(measured, truth, sigma)
    
    # canopy_perc_dif, n_canopy = compute_perc_canopy(truth_ground, truth_canopy,
    #                                                 atl03_canopy, sigma)


    # class_list = np.unique(sb.segment_landcover)
    # out_arr = np.zeros((len(class_list),35))
    # import copy
    # sb_original = copy.deepcopy(sb)
    # for j in class_list:
    #     sb = copy.deepcopy(sb_original)
    #     sb.segment_landcover[j] = np.nan
    #     sb_new = sb.dropna(subset=['segment_landcover'], inplace=False)
    #     truth_canopy = sb_new.truth_canopy_max98
    #     truth_ground = sb_new.truth_ground_median
    #     atl03_canopy = sb_new.atl03_canopy_max98
        
    #     truth_canopy = truth_canopy - truth_ground
    #     atl03_canopy = atl03_canopy - truth_ground
    #     sigma = 1
        
    #     truth_canopy[truth_canopy > 50] = np.nan
    #     truth_canopy[truth_canopy < 0] = np.nan
    #     nan1 = np.isnan(truth_canopy)
    #     truth_canopy = truth_canopy[~nan1]
    #     atl03_canopy = atl03_canopy[~nan1]   
    #     atl03_canopy[atl03_canopy > 50] = np.nan
    #     atl03_canopy[atl03_canopy < 0] = np.nan
    #     nan1 = np.isnan(atl03_canopy)
    #     truth_canopy = truth_canopy[~nan1]
    #     atl03_canopy = atl03_canopy[~nan1]
        
    #     canopy_res = truth_canopy - atl03_canopy
        
    #     if sigma > 0:
    #         mean = np.nanmean(canopy_res)
    #         std = np.nanstd(canopy_res)
    #         mask = ~np.isnan(canopy_res)
    #         sigmathresh = mean + (std * sigma)
    #         mask[mask] &= canopy_res[mask] > sigmathresh
    #         canopy_res[mask] = np.nan
    #         truth_canopy[mask] = np.nan
    #         atl03_canopy[mask] = np.nan
    #         mask = ~np.isnan(canopy_res)
    #         sigmathresh = mean - (std * sigma)
    #         mask[mask] &= canopy_res[mask] < sigmathresh
    #         canopy_res[mask] = np.nan   
    #         truth_canopy[mask] = np.nan
    #         atl03_canopy[mask] = np.nan 
            
    #     nan1 = np.isnan(truth_canopy)
    #     truth_canopy = truth_canopy[~nan1]
    #     atl03_canopy = atl03_canopy[~nan1]  
    
    #     nan1 = np.isnan(atl03_canopy)
    #     truth_canopy = truth_canopy[~nan1]
    #     atl03_canopy = atl03_canopy[~nan1]
    
        

        
    #     for i in range(0,35):
    #         print(i)
    #         test = ((truth_canopy >= i ) & (truth_canopy <= i + 1))
    #         perc_dif = 1 - (np.median(atl03_canopy[test]) / np.median(truth_canopy[test]))
    #         print(perc_dif)
    #         # perc_dif_list_median.append(perc_dif)
    #         out_arr[j,i] = perc_dif

