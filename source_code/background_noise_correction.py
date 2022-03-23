#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 08:26:42 2021

@author: eguenther
"""

import argparse
import os
from shutil import copyfile    
import h5py
import numpy as np
import pandas as pd

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

def ismember(a_vec, b_vec, methodType = 'normal'):
    
    """ MATLAB equivalent ismember function """
    
    # Combine multi column arrays into a 1-D array of strings if necessary
    # This will ensure unique rows when using np.isin below
    if(methodType.lower() == 'rows'):
        
        # Turn a_vec into an array of strings
        a_str = a_vec.astype('str')
        b_str = b_vec.astype('str')
        
        # Concatenate each column of strings with commas into a 1-D array
        for i in range(0,np.shape(a_str)[1]):
            a_char = np.char.array(a_str[:,i])
            b_char = np.char.array(b_str[:,i])
            if(i==0):
                a_vec = a_char
                b_vec = b_char
            else:
                a_vec = a_vec + ',' + a_char
                b_vec = b_vec + ',' + b_char
    
    # Find which values in a_vec are present in b_vec
    matchingTF = np.isin(a_vec,b_vec)
    common = a_vec[matchingTF]
    common_unique, common_inv  = np.unique(common, return_inverse=True)    
    b_unique, b_ind = np.unique(b_vec, return_index=True)
    common_ind = b_ind[np.isin(b_unique, common_unique, assume_unique=True)]
    matchingInds = common_ind[common_inv]
    
    return matchingTF, matchingInds

def calc_radiometry(atl03filepath, atl08filepath, gt):
    # Read ATL03 .h5 file
    f = h5py.File(atl03filepath, 'r')
    segment_ph_count = np.array(f[gt + '/geolocation/segment_ph_cnt'])
    atl03_segment_id = np.array(f[gt + '/geolocation/segment_id'])
    atl03_ph_index_beg = np.array(f[gt + '/geolocation/ph_index_beg'])    
    bihr = np.asarray(f[gt + '/bckgrd_atlas/bckgrd_int_height_reduced'])
    bcr = np.asarray(f[gt + '/bckgrd_atlas/bckgrd_counts_reduced'])
    bapmc = np.asarray(f[gt + '/bckgrd_atlas/pce_mframe_cnt'])
    h_ph = np.asarray(f[gt + '/heights/h_ph'])
    hpmc = np.asarray(f[gt + '/heights/pce_mframe_cnt'])
    delta_time = np.asarray(f[gt + '/heights/delta_time'])   

    # Get segment ID to photon level
    h_seg = np.zeros(len(h_ph))
    for i in range(0,len(atl03_segment_id)):
        if atl03_ph_index_beg[i] > 0:
            h_seg[atl03_ph_index_beg[i]-1:atl03_ph_index_beg[i]-1 +\
                      segment_ph_count[i]] = atl03_segment_id[i]
    h_seg = np.int32(h_seg) 
   
    # Calculate rate and assign to the photon level
    rate = bcr / bihr

    # Assign bckgrd_atlas attributes to photon level    
    tf, inds = ismember(hpmc, bapmc)
    ph_bihr = bihr[inds]
    ph_bcr = bcr[inds]
    ph_rate = rate[inds]

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
    
    # Put everything in ATL03 into pandas df
    df = pd.DataFrame({'h_ind': h_ind, 
                       'delta_time': delta_time,
                       'bckgrd_int_height_reduced':ph_bihr,
                       'bckgrd_counts_reduces': ph_bcr,
                       'bckgrd_rate': ph_rate})
    
    # Remove ATL03 files that do not fit into ATL08 segment
    df = df[df['h_ind'] >= 0]
    
    # Average Background Rate to ATL08
    zgroup = df.groupby('h_ind')
    zout = zgroup.aggregate(np.nanmean)
    f08_rate = np.asarray(zout['bckgrd_rate'])
    # f08_bihr = np.asarray(zout['bckgrd_int_height_reduced'])
    # f08_bcr = np.asarray(zout['bckgrd_counts_reduces'])
    # f08_rh = f08[gt + '/land_segments/segment_id_beg']
    h_max_canopy = np.asarray(f08[gt + '/land_segments/canopy/h_max_canopy'])
    h_max_canopy[h_max_canopy > 100000] = np.nan
    canopy_noise_count = (f08_rate) * h_max_canopy
    
    # Count the number of photons at the ATL08 segment
    n_ca_photons = np.asarray(f08[gt + '/land_segments/canopy/n_ca_photons'])
    n_toc_photons = np.asarray(f08[gt + '/land_segments/canopy/n_toc_photons'])
    n_te_photons = np.asarray(f08[gt + '/land_segments/terrain/n_te_photons'])
    
    # Calcualte 'noise adjusted' number of canopy photons at ATL08 segment 
    n_ca_photons_nr = (n_ca_photons + n_toc_photons) - canopy_noise_count
    
    # Calculate photon rate (radiometry)
    zgroup = df.groupby('h_ind')
    zout = zgroup.aggregate(get_len_unique)
    f08_shots = np.asarray(zout['delta_time'])

    zgroup = df.groupby('h_ind')
    zout = zgroup.aggregate(get_len)
    # n_photons = np.asarray(zout['delta_time'])
    
    photon_rate_can_nr = n_ca_photons_nr / f08_shots
    photon_rate_can = (n_ca_photons + n_toc_photons) / f08_shots
    photon_rate_can[np.isnan(photon_rate_can_nr)] = np.nan
    photon_rate_signal = (n_te_photons + n_ca_photons + n_toc_photons)\
        / f08_shots
    photon_rate_ground = n_te_photons / f08_shots
    # photon_rate_total  = (n_photons) / f08_shots
    # photon_rate_noise = (n_photons -\
        # (n_te_photons + n_ca_photons +\ n_toc_photons)) / f08_shots
    
    return n_ca_photons_nr, photon_rate_can_nr, photon_rate_can,\
        photon_rate_signal, photon_rate_ground
        
def main(atl03filepath, atl08filepath, out_dir):
    atl08file = os.path.basename(atl08filepath)
    print(atl08file)
    newfilename = 'ATL08_nr_' + atl08file.split('ATL08_')[1]
    newatl08file = os.path.join(out_dir, newfilename)
    copyfile(atl08filepath,newatl08file)
    h5f = h5py.File(newatl08file,'a')
    gt_list = ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']
    for gt in gt_list:
        base_key = gt + '/land_segments/'

        n_ca_photons_nr, photon_rate_can_nr, photon_rate_can,\
            photon_rate_signal, photon_rate_ground =\
            calc_radiometry(atl03filepath, atl08filepath, gt)
        
        h5f[base_key + 'n_ca_photons_nr'] = n_ca_photons_nr
        h5f[base_key + 'photon_rate_can_nr'] = photon_rate_can_nr
        h5f[base_key + 'photon_rate_can'] = photon_rate_can
        h5f[base_key + 'photon_rate_signal'] = photon_rate_signal
        h5f[base_key + 'photon_rate_ground'] = photon_rate_ground
            
if __name__ == '__main__':
    """ Command line entry point """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("atl03_filepath", metavar="atl03",
                        help="Input ATL03 filepath (.h5)")

    parser.add_argument("atl08_filepath", metavar="atl08",
                        help="Input ATL08 filepath (.h5)")

    parser.add_argument("output_directory", metavar="out_dir",
                        help="Output directory of file")

    args = parser.parse_args()
    main(args.atl03_filepath, args.atl08_filepath, args.output_directory)
    # if os.name == 'nt':
    #     basepath03 = 'Z:/data/release/002/ATL03_r002/Finland/'
    #     basepath08 = 'Z:/data/release/002/ATL08_r002/Finland/'
    # else:
    #     basepath03 = '/laserpewpew/data/release/004/ATL03_r004/Finland_nsidc/'
    #     basepath08 = '/laserpewpew/data/release/004/ATL08_r004/Finland_nsidc/'  

    # outFilePath = '/LIDAR/server/USERS/eric/for_amy/finland_lai_100/'
    # try:
    #     os.mkdir(outFilePath)
    # except:
    #     print('Folder already exists')

    # atl03_list = os.listdir(basepath03)

    # i = 0
    # atl03file = atl03_list[i]
    # atl08file = 'ATL08' + atl03_list[i].split('ATL03')[1]
    # atl03filepath =  os.path.join(basepath03, atl03file)
    # atl08filepath =  os.path.join(basepath08, atl08file)
    # gt = 'gt1l'
    # n_ca_photons_nr, photon_rate_can_nr, photon_rate_can, photon_rate_signal,\
    #     photon_rate_ground = calc_radiometry(atl03filepath, atl08filepath, gt)