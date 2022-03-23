#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 13:26:11 2021

@author: eguenther
"""

import os
import time
import numpy as np
import pandas as pd
import scipy
# from getAtlTruthSwath_auto import getAtlTruthSwath
# from getMeasurementError_auto import getMeasurementError, offsetsStruct


# readGeoidFile , getGeoidHeight 
from icesatReader import get_atl03_struct
from icesatReader import get_atl08_struct
from icesatReader import convert_atl03_to_legacy
from icesatReader import get_atl_alongtrack
from icesatReference import legacy_get_truth_swath, legacy_get_meas_error
from icesatIO import getTruthHeaders
from icesatIO import readGeoidFile
from icesatUtils import getGeoidHeight
from icesatUtils import superFilter
from icesatBin import create_atl08_bin
from icesatBin import create_truth_bin
from icesatIO import read_geotiff
from icesatIO import find_intersecting_values
from icesatUtils import getAtl08Mapping
from icesatUtils import getAtl08Mapping_seg
from icesatIO import readAtl03DataMapping
from icesatIO import readAtl08DataMapping
from icesatUtils import ismember
from icesatBin import get_len_unique
import h5py

def calc_radiometry(atl03filepath, atl08filepath, gt):

    
    atl03_ph_index_beg, atl03_segment_id = \
    readAtl03DataMapping(atl03filepath,gt)
    
    # Read ATL08 for class mapping
    atl08_classed_pc_indx, atl08_classed_pc_flag, atl08_segment_id = \
    readAtl08DataMapping(atl08filepath, gt)
    
    allph_classed = getAtl08Mapping(atl03_ph_index_beg, atl03_segment_id, 
                                    atl08_classed_pc_indx, atl08_classed_pc_flag, 
                                    atl08_segment_id)
    
    allph_seg = getAtl08Mapping_seg(atl03_ph_index_beg, atl03_segment_id, 
                                    atl08_classed_pc_indx, atl08_classed_pc_flag, 
                                    atl08_segment_id)
    
    
    f = h5py.File(atl03filepath,'r')
    bihr = np.asarray(f[gt + '/bckgrd_atlas/bckgrd_int_height_reduced'])
    bcr = np.asarray(f[gt + '/bckgrd_atlas/bckgrd_counts_reduced'])
    
    badt = np.asarray(f[gt + '/bckgrd_atlas/delta_time'])
    bapmc = np.asarray(f[gt + '/bckgrd_atlas/pce_mframe_cnt'])
    h_ph = np.asarray(f[gt + '/heights/h_ph'])
    lat = np.asarray(f[gt + '/heights/lat_ph'])
    lon = np.asarray(f[gt + '/heights/lon_ph'])
    across = np.asarray(f[gt + '/heights/dist_ph_across'])
    along = np.asarray(f[gt + '/heights/dist_ph_along'])
    hpmc = np.asarray(f[gt + '/heights/pce_mframe_cnt'])
    delta_time = np.asarray(f[gt + '/heights/delta_time'])
    rate = bcr / bihr
    
    tf, inds = ismember(hpmc, bapmc)
    ph_bihr = bihr[inds]
    ph_badt = badt[inds]
    ph_bcr = bcr[inds]
    ph_rate = rate[inds]
    
    fill = np.zeros(len(h_ph) - len(allph_classed))
    fill = fill - 1
    allph_classed = np.concatenate([allph_classed,fill])
    
    atl03_seg = np.zeros(len(h_ph))
    atl03_seg[atl03_ph_index_beg] = atl03_segment_id
    df = pd.DataFrame(data=atl03_seg, columns = ['A'])
    h_seg = np.asarray(df.replace(to_replace=0,method='ffill')).astype(int)
    h_seg[0] = h_seg[1]
    h_seg = h_seg.T[0]
    
    
    f08 = h5py.File(atl08filepath,'r')
    seg_beg = np.asarray(f08[gt + '/land_segments/segment_id_beg']).astype(int)
    seg_end = np.asarray(f08[gt + '/land_segments/segment_id_end']).astype(int)
    longitude = np.asarray(f08[gt + '/land_segments/longitude'])
    latitude = np.asarray(f08[gt + '/land_segments/latitude'])
    # Read ATL03 metrics for class mapping
    
    start = np.min(seg_beg)
    
    
    
    h_ind = np.zeros(len(h_seg))
    h_ind[:] = -1
    
    for i in range(0,len(seg_beg)):
        h_ind[np.logical_and(h_seg >= seg_beg[i], h_seg <= seg_end[i])] = i
    
    
    # f08_ind = seg_beg - start
    # f08_ind = np.floor(f08_ind / 5).astype(int)
    
    
    # Put everything in ATL03 into pandas df
    df = pd.DataFrame({'h_ph': h_ph, 'h_ind': h_ind, 
                       'c': allph_classed, 'along': along,
                       'delta_time': delta_time,
                       'bckgrd_int_height_reduced':ph_bihr,
                       'bckgrd_counts_reduces': ph_bcr,
                       'bckgrd_rate': ph_rate})
    
    df = df[df['h_ind'] >= 0]
    
    
    # Average Background Rate
    zgroup = df.groupby('h_ind')
    zout = zgroup.aggregate(np.nanmean)
    f08_rate = np.asarray(zout['bckgrd_rate'])
    f08_bihr = np.asarray(zout['bckgrd_int_height_reduced'])
    f08_bcr = np.asarray(zout['bckgrd_counts_reduces'])
    
    f08_rh = f08[gt + '/land_segments/segment_id_beg']
    h_max_canopy = np.asarray(f08[gt + '/land_segments/canopy/h_max_canopy'])
    # h_max_canopy[h_max_canopy > 100000] = np.nan
    
    # canopy_noise_count = (f08_rate/50) * h_max_canopy
    canopy_noise_count = (f08_rate) * h_max_canopy
    
    n_ca_photons = np.asarray(f08[gt + '/land_segments/canopy/n_ca_photons'])
    n_toc_photons = np.asarray(f08[gt + '/land_segments/canopy/n_toc_photons'])
    n_te_photons = np.asarray(f08[gt + '/land_segments/terrain/n_te_photons'])
    n_ca_photons_nr = (n_ca_photons + n_toc_photons) - canopy_noise_count
    
    zgroup = df.groupby('h_ind')
    zout = zgroup.aggregate(get_len_unique)
    f08_shots = np.asarray(zout['delta_time'])
    
    photon_rate_can_nr = n_ca_photons_nr / f08_shots
    photon_rate_can = (n_ca_photons + n_toc_photons) / f08_shots
    photon_rate_can[np.isnan(photon_rate_can_nr)] = np.nan
    photon_rate_all = (n_te_photons + n_ca_photons + n_toc_photons) / f08_shots
    photon_rate_ground = n_te_photons / f08_shots
    
    return canopy_noise_count, n_ca_photons_nr, f08_shots, f08_bcr, f08_bihr, f08_rate,\
        photon_rate_can_nr, photon_rate_can, photon_rate_all, photon_rate_ground

if __name__ == "__main__":    
    if os.name == 'nt':
        basepath03 = 'Z:/data/release/002/ATL03_r002/Finland/'
        basepath08 = 'Z:/data/release/002/ATL08_r002/Finland/'
    else:
        basepath03 = '/laserpewpew/data/release/004/ATL03_r004/North_Carolina_nsidc/'
        basepath08 = '/laserpewpew/data/release/004/ATL08_r004/North_Carolina_nsidc/' 
    
    # atl03file = 'ATL03_20181118120428_07770103_002_01.h5'
    # atl08file = 'ATL08_20181118120428_07770103_002_01.h5'
    
    atl03file_list = os.listdir(basepath03)
    # Inputs
    atl03file = atl03file_list[1]
    for atl03file in atl03file_list:
        try:
            atl03filepath = basepath03 + atl03file
            atl08file = 'ATL08_' + atl03file.split('ATL03_')[1]
            atl08filepath = basepath08 + atl08file
            gt_list = ['gt1r','gt1l','gt2r','gt2l','gt3r','gt3l']
            shortname = atl08file.split('.h5')[0]
            for gt in gt_list:
                
                header_file_path =\
                    '/LIDAR/server/USERS/eric/1_experiment/North_Carolina_HeaderData.mat'
                    
                kml_bounds_txt1 = '/LIDAR/server/USERS/eric/2_production/kmlBounds.txt'
                   
                print('Generate ATL03 Struct')
                atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath, 
                                         epsg = '32617', kml_bounds_txt = kml_bounds_txt1, 
                                         header_file_path = header_file_path)
                
                print('Generate ATL08 Struct')
                atl08 = get_atl08_struct(atl08filepath, gt, atl03)
                
                canopy_noise_count, n_ca_photons_nr, f08_shots, f08_bcr, f08_bihr, f08_rate, photon_rate_can_nr,\
            photon_rate_can, photon_rate_all, photon_rate_ground = calc_radiometry(atl03filepath, atl08filepath, gt)
                
                atl08.df['photon_rate_can'] = photon_rate_can
                atl08.df['photon_rate_can_nr'] = photon_rate_can_nr
                atl08.df['photon_rate_ground'] = photon_rate_ground
                # atl08.df['atl03_radiometry_ground'] = photon_rate_ground
                atl08.df['photon_rate_all'] = photon_rate_all
                atl08.df['avg_bckgrd_calc_rate'] = f08_rate
                atl08.df['avg_bckgrd_counts_reduces'] = f08_bcr
                atl08.df['avg_bckgrd_int_height_reduced'] = f08_bihr
            
                #atl03.df = atl03.df[atl03.df['time'] < 12]
                
                df, rotation_data = get_atl_alongtrack(atl03.df)
                
                atl03.df = df
                atl03.rotationData = rotation_data
                
                print('Convert Struct to Legacy')    
                atl03legacy, rotationData = convert_atl03_to_legacy(atl03)
                
                # Legacy Truth Swath Inputs
                buffer = 50                 # Distance in cross-track (meters) around ATL03 track to look for truth data 
                useExistingTruth = False     # Option to use existing truth data if it exists
                # truthSwathDir = 'N:/data/TanDEMX//MississippiValley_TDR_DSM/Indiana/TDT_N37W086_02_DEM.TIF'
                # truthSwathDir = '/laser1/validation/data/Finland/LAZ_UTM'
                truthSwathDir = '/laser1/validation/data/NorthCarolina/LAS_UTM2'
                outFilePath = '/LIDAR/server/USERS/eric/1_experiment/ecosystem_nc/'
                createTruthFile = True      # Option to create output truth .las file
                
                
                timeStart = time.time()
                
                # Call getAtlMeasuredSwath
                # print('RUNNING getAtlMeasuredSwath...\n')
                # atl03Data, atl08Data, headerData, rotationData = getAtlMeasuredSwath(atl03FilePath, atl08FilePath, outFilePath, gtNum, trimInfo, createAtl03LasFile, createAtl03KmlFile, createAtl08KmlFile, createAtl03CsvFile, createAtl08CsvFile)
                    
                # Call getAtlTruthSwath
                print('Run Legacy Truth Swath')
                truthFileType = '.las'
                atlTruthData = legacy_get_truth_swath(atl03legacy, atl03.rotationData, 
                                                      truthSwathDir, truthFileType, 
                                        outFilePath, buffer = 20, useExistingTruth = False,
                                        createTruthFile = True)
                
                
                
                # End timer
                
                
                
                geoidDataFile = '/laserpewpew/data/validation/geoid/GEOID12B/geoid12b_latlon.mat'
                geoidData = readGeoidFile(geoidDataFile)
                atlTruthData = getGeoidHeight(geoidData,atlTruthData)
                timeEnd = time.time()
                timeElapsedTotal = timeEnd - timeStart
                timeElapsedMin = np.floor(timeElapsedTotal / 60)
                timeElapsedSec = timeElapsedTotal % 60
                    
                # Print completion message
                print('   Script Completed in %d min %d sec.' % (timeElapsedMin, 
                                                                  timeElapsedSec))
                print('\n')
                
                print('   ATL Corrections.')
                
                atlCorrections = legacy_get_meas_error(atl03legacy, atlTruthData, 
                                                        rotationData, outFilePath)
                
                print('Superfilter Complete')
                atlTruthData, atl03legacy = superFilter(atl03legacy, atlTruthData, 
                                                        xBuf = 7, classCode = [], verbose=False)
                
                ##Filter Truth Data by Include
                alongtrack = atlTruthData.alongTrack.flatten()
                crosstrack = atlTruthData.crossTrack.flatten()
                z = atlTruthData.z.flatten()
                easting = atlTruthData.easting.flatten()
                northing = atlTruthData.northing.flatten()
                classification = atlTruthData.classification.flatten()
                intensity = atlTruthData.intensity.flatten()
                
                
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
                
                # Compute ATL08 Bin Heights
                steps = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
                
                for i in steps:
                    df_bin = create_atl08_bin(atl03, atl08, res_at = i)
                    
                    # Compute Truth Bin Heights
                    df_truth_bin = create_truth_bin(df_truth, atl08, res_at = i)
                    
                    df_bin_out = pd.merge(df_bin, df_truth_bin, on="bin_id",how='left')
                    
                    #Write out Pickle File
                    pkl_file = outFilePath + 'pkl/' + shortname + '_' + gt + '_' + str(i) +'m.pkl'
                    df_bin_out.to_pickle(pkl_file)
                    
                    #Write out Mat file
                    a_dict = {col_name : df_bin_out[col_name].values for col_name \
                              in df_bin_out.columns.values}
                    mat_file = outFilePath + 'mat/' + shortname + '_' + gt + '_' + str(i) +'m.mat'
                    scipy.io.savemat(mat_file,{'struct':a_dict})
        except:
                print('fail')