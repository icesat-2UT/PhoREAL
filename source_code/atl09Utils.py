#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that contains functions for the processing and analysis of ATL09

Copyright 2020 Applied Research Laboratories, University of Texas at Austin

This package is free software; the copyright holder gives unlimited
permission to copy and/or distribute, with or without modification, as
long as this notice is preserved.

Authors:
    Jonathan Markel

Date: March 24, 2020
"""

# %% Imports
from icesatReader import get_atl03_struct
from icesatReader import get_atl08_struct
from icesatReader import get_atl09_struct
from icesatReader import read_atl03_geolocation
from icesatReader import append_atl03_geolocation

from icesatIO import readAtlH5
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.stats import skew, kurtosis, shapiro
from tqdm import tqdm


def get_beam_feature(atl09):

    if atl09.orbit_info.sc_orient == 0:
        beam = np.full(atl09.df.shape[0], int(atl09.gtNum[2]))

    else:
        if atl09.gtNum[2] == '1':
            beam = np.full(atl09.df.shape[0], 3)
        elif atl09.gtNum[2] == '2':
            beam = np.full(atl09.df.shape[0], 2)
        else:
            beam = np.full(atl09.df.shape[0], 1)

    return beam


def plot_cab(atl09, x_column='northing', remove_surface=False,
             show_layers=True, show_detected_surface=True, show_dem=True,
             overlay_alpha=1):

    if atl09.df.solar_elevation.mean() > -7:
        cmap_range = (1e-10, 1.1e-5)
    else:
        cmap_range = (1e-10, 0.9e-5)

    if remove_surface:
        cab = np.copy(atl09.cab_no_surf)
    else:
        cab = np.copy(atl09.cab_prof)

    cab[cab > 1e30] = 0
    begin_val = atl09.df.loc[:, x_column][0]
    end_val = atl09.df.loc[:, x_column][-1:].iloc[0]

    if begin_val > end_val:
        x = np.flip(atl09.df.loc[:, x_column])
        y = atl09.ds_va_bin_h
        z = np.fliplr(cab.T)

    else:
        x = atl09.df.loc[:, x_column]
        y = atl09.ds_va_bin_h
        z = (cab.T)

    fig, ax = plt.subplots(figsize=(12, 6))

    cmap = plt.get_cmap('viridis')  # plasma, inferno also good
    im = plt.pcolormesh(x,
                        y,
                        z,
                        vmin=cmap_range[0],
                        vmax=cmap_range[1],
                        cmap=cmap)

    fig.colorbar(im, ax=ax)

    styles = ['w.', 'w.', 'w.',
              'w.', 'w.', 'w.',
              'w.', 'w.', 'w.',
              'r.']
    if show_layers:
        for layer_num, style in zip(np.arange(10), styles):
            top_col = 'layer_top_' + str(layer_num)
            bot_col = 'layer_bot_' + str(layer_num)
            mask = atl09.df.loc[:, top_col] < 1e30
            if layer_num == 0:
                plt.plot(atl09.df.loc[mask, x_column],
                         atl09.df.loc[mask, top_col],
                         style,
                         alpha=overlay_alpha,
                         label='DDA Layers')
            else:
                plt.plot(atl09.df.loc[mask, x_column],
                         atl09.df.loc[mask, top_col],
                         style,
                         alpha=overlay_alpha)

            plt.plot(atl09.df.loc[mask, x_column],
                     atl09.df.loc[mask, bot_col],
                     style,
                     alpha=overlay_alpha)

    if show_dem:
        plt.plot(atl09.df.loc[:, x_column],
                 atl09.df.dem_h,
                 '.',
                 color='lightsteelblue',
                 alpha=overlay_alpha,
                 label='DEM')

    if show_detected_surface:
        mask = atl09.df.surface_bin < 1e3
        heights = atl09.df.loc[mask, 'surface_bin'].apply(lambda x: atl09.ds_va_bin_h[x])
        plt.plot(atl09.df.loc[mask, x_column],
                 heights,
                 '.',
                 color='red',
                 alpha=1,
                 label='ATL09 Detected Surface')

    plt.ylim([-1000, 17500])

    if begin_val > end_val:
        ax = plt.gca()
        xlim = ax.get_xlim()
        plt.xlim([xlim[1], xlim[0]])

    plt.xlabel(x_column.upper())
    plt.ylabel('Elevation (m)')
    plt.title(atl09.atlFileName + ', Profile ' + atl09.gtNum[2])
    plt.legend()
    ax = plt.gca()

    return ax


def smooth(x, window_len=11, window='hanning'):
    """
    From scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    """

    if x.ndim != 1:
        print("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        print("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        print("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]

    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y


def get_cab_no_surf(atl09,
                    show_plots=False,
                    idxing=None,
                    column_search_radius=5,
                    track_search_radius=5,
                    window_len=7,
                    peak_threshold=2e-5):

    if idxing is None:
        idxing = atl09.df.index

    cab_no_surf = pd.DataFrame(np.copy(atl09.cab_prof), index=atl09.df.index)
    cab_no_surf[cab_no_surf > 1e30] = 0
    # Surface estimate algorithm
    cab = pd.DataFrame(atl09.cab_prof,
                       index=atl09.df.index)
    cab[cab > 1e30] = 0

    cutoff_bins = np.full(atl09.df.shape[0], 700)
    surf_removal_method = np.zeros_like(cutoff_bins)

    # Don't accidentally terrorize matplotlib
    if idxing.shape[0] > 20:
        show_plots = False

    for idx in tqdm(idxing):
        # Get surface bin and CAB profile
        surf_bin_i = atl09.df.loc[idx, 'surface_bin'] - 1
        cab_i = cab.loc[idx, :]
        lowest_layer_h = min(atl09.df.loc[idx,
                                          ['layer_bot_' + str(x)
                                            for x in np.arange(10)]])
        lowest_layer_bin = np.around((max(atl09.ds_va_bin_h) - lowest_layer_h) / 30) + 2

        if surf_bin_i < 701:
            # If surface was found by ATL09, confirm location with CAB profile

            # too overzealous with low clouds
            # surf_bin_check = cab_i.iloc[surf_bin_i - column_search_radius:\
            #                             surf_bin_i + column_search_radius].idxmax()
            # cutoff_bin = int(surf_bin_check - 1)  # 30m buffer

            cutoff_bin = int(surf_bin_i - 3)
            surf_removal_method[idx] = 1

            if show_plots:
                # Visualize
                plt.figure()
                plt.title(str(idx) + ', ATL09 Located surface')
                plt.plot(cab_i)
                # plt.axvline(surf_bin_check, color='green', label='Found surface')
                plt.axvline(surf_bin_i, color='blue', label='ATL09 surface')
                plt.axvline(cutoff_bin, color='red', label='Cutoff point')
                plt.axvline(surf_bin_i + column_search_radius, color='yellow')
                plt.axvline(surf_bin_i - column_search_radius, color='yellow')
                plt.xlim([surf_bin_i - 10, surf_bin_i + 10])
                plt.legend()

        else:
            # If ATL09 couldn't find surface
            # See if ATL09 found the surface for nearby neighbors
            local_surf_bin = atl09.df.loc[idx - track_search_radius:
                                          idx + track_search_radius, 'surface_bin'] - 1
            min_local_bin = local_surf_bin[local_surf_bin < 701].min()
            max_local_bin = local_surf_bin[local_surf_bin < 701].max()
            mean_local_bin = local_surf_bin[local_surf_bin < 701].mean()

            if np.isnan(mean_local_bin):
                # IF no peaks identified,
                # OR no singular peak identified
                # OR no nearby surface detected
                # Use DEM to estimate surface location
                # Check granule-level surface height difference
                # Use 2 sigma difference above DEM for cutoff
                # Check DEM for estimated surface location
                dem_bin = np.around((max(atl09.ds_va_bin_h) - atl09.df.loc[idx, 'dem_h']) / 30) + 1
                error_h = atl09.df.dem_h[atl09.df.surface_bin < 800] - \
                    atl09.df.surface_height[atl09.df.surface_bin < 800]
                if len(error_h) == 0:
                    # ATL09 did not find any surface points in this granule
                    cutoff_bin = int(dem_bin)
                else:
                    stdev = np.std(error_h)
                    dem_above_bins = np.around((error_h.mean() + 2 * stdev) / 30)
                    dem_below_bins = np.around((error_h.mean() - 2 * stdev) / 30)

                    bin_2sigma_above_dem = dem_bin + dem_below_bins

                    cutoff_bin = int(bin_2sigma_above_dem)
                    surf_removal_method[idx] = 4

                if show_plots:
                    plt.figure()
                    plt.title(str(idx) + ', DEM Approach')
                    plt.plot(cab_i, label='CAB')
                    plt.axvline(dem_bin, label='DEM Surface', color='magenta')
                    plt.axvline(cutoff_bin, label='Cutoff', color='red')
                    plt.xlim([min([cutoff_bin, dem_bin]) - 30, max([cutoff_bin, dem_bin]) + 30])
                    plt.legend()

            else:
                # Use neighbors to define window to search for peak in CAB
                cab_smooth = smooth(cab_i, window_len=window_len)
                bin_smooth = np.linspace(0, 700 + window_len / 3, 699 + window_len)
                cab_search = cab_smooth[min_local_bin - column_search_radius:
                                        max_local_bin + column_search_radius]

                pk_idx, pk_h = signal.find_peaks(cab_search, height=peak_threshold)

                if len(pk_idx) == 1:
                    # IF a peak is found based on CAB, define this as surface
                    peak = np.around(bin_smooth[min_local_bin - column_search_radius + pk_idx - 1]) - 1
                    cutoff_bin = int(peak)
                    surf_removal_method[idx] = 3
                    if show_plots:
                        plt.figure()
                        plt.title(str(idx) + ', Smoothed peak finding')
                        plt.plot(cab_i, label='CAB')
                        plt.plot(bin_smooth, cab_smooth, label='Smoothed CAB profile')
                        plt.axvline(np.around(mean_local_bin), color='cyan', label='Avg nearby surface location')
                        plt.axvline(cutoff_bin, color='green', label='Cutoff bin')
                        plt.axvline(min_local_bin - column_search_radius, color='red', label='Search window')
                        plt.axvline(min_local_bin + column_search_radius, color='red')
                        plt.axvline(surf_bin_i - column_search_radius, color='red')
                        plt.axvline(lowest_layer_bin, color='black')
                        plt.xlim([min([cutoff_bin, lowest_layer_bin, min_local_bin, surf_bin_i]) - 100,
                                  max([cutoff_bin, lowest_layer_bin, min_local_bin, surf_bin_i]) + 100])
                        plt.legend()

                else:
                    # IF no clear peak identified, USE DEM approach
                    dem_bin = np.around((max(atl09.ds_va_bin_h) - atl09.df.loc[idx, 'dem_h']) / 30) + 1
                    error_h = atl09.df.dem_h[atl09.df.surface_bin < 800] - \
                        atl09.df.surface_height[atl09.df.surface_bin < 800]
                    if len(error_h) == 0:
                        # ATL09 did not find any surface points in this granule
                        cutoff_bin = int(dem_bin)
                    else:
                        stdev = np.std(error_h)
                        dem_above_bins = np.around((error_h.mean() + 2 * stdev) / 30)
                        dem_below_bins = np.around((error_h.mean() - 2 * stdev) / 30)

                        bin_2sigma_above_dem = dem_bin + dem_below_bins

                        cutoff_bin = int(bin_2sigma_above_dem)
                        surf_removal_method[idx] = 4

                    if show_plots:
                        plt.figure()
                        plt.title(str(idx) + ', DEM approach, finding peaks w smooth CAB')
                        plt.plot(cab_i, label='CAB')
                        found_peaks = np.around(bin_smooth[min_local_bin - column_search_radius + pk_idx - 1])
                        if len(found_peaks) > 1:
                            for pk in found_peaks:
                                if pk == found_peaks[0]:
                                    plt.axvline(pk, color='cyan', label='Peaks found')
                                else:
                                    plt.axvline(pk, color='cyan')

                        plt.axvline(lowest_layer_bin, color='black', label='Lowest cloud layer')
                        plt.axvline(dem_bin, label='DEM Surface', color='magenta')
                        plt.axvline(cutoff_bin, label='Cutoff', color='red')
                        plt.xlim([min([min(found_peaks), cutoff_bin, lowest_layer_bin]) - 15,
                                  max([max(found_peaks), cutoff_bin, lowest_layer_bin]) + 15])
                        plt.legend()

        if np.abs(lowest_layer_bin - cutoff_bin) < 5:
            # safety check if clouds are low
            cutoff_bin = int(lowest_layer_bin + 5)
            surf_removal_method[idx] = 5

        cutoff_bins[idx] = cutoff_bin
        cab_no_surf.iloc[idx, cutoff_bin:] = 0


    return cab_no_surf.iloc[idxing], surf_removal_method
def get_ml_features(atl09):
    df_ml = pd.DataFrame(index=atl09.df.index,
                         columns=['profile', 'beam', 'file', 'cab_int',
                                  'cab_skew', 'cab_kurt', 'cab_shapiro'])

    df_ml.loc[:, 'profile'] = np.full(df_ml.shape[0], int(atl09.gtNum[2]))
    df_ml.loc[:, 'file'] = np.full(df_ml.shape[0], atl09.atlFileName)
    df_ml.loc[:, 'beam'] = get_beam_feature(atl09)
    df_ml.loc[:, 'cab_skew'] = skew(atl09.cab_no_surf, axis=1)
    df_ml.loc[:, 'cab_kurt'] = kurtosis(atl09.cab_no_surf, axis=1)
    df_ml.loc[:, 'cab_int'] = -np.trapz(atl09.cab_no_surf, x=atl09.ds_va_bin_h)
    df_ml.loc[:, 'cab_shapiro'] = [i[0] for i in atl09.cab_no_surf.apply(shapiro, axis=1)]

    return df_ml


def get_rgt_matches(input_file_name, filenames_to_search, use_full_rgtcs=False):
    # Find the set of 8 numbers that is RGT, cycle, and segment data, in the filename, store rgtcc
    if use_full_rgtcs:
        rgtcc_str = '_(\d{8})_'
    else:
        rgtcc_str = '_(\d{6})\d{2}_'

    rgtcc_match = re.search(rgtcc_str, input_file_name)
    matches = dict()
    if rgtcc_match == None:
        print('Could not identify rgt for input file')
        matches = []

    else:
        # Get the 4-digit RGT CC combo
        rgtcc = rgtcc_match.group(1)

        # Search for RGTCC set in input list
        if use_full_rgtcs:
            search_str = '_(' + rgtcc + ')_'
        else:
            search_str = '_(' + rgtcc + '\d{2})_'

        # Returns list of matching files
        found = [x for x in filenames_to_search if re.search(search_str, x)]
        matches[input_file_name] = found

    return matches


def check_for_matching_segments(atl09file, atl03file_list, atl09basepath, atl03basepath):
    atl03_overlap_with_atl09 = pd.DataFrame(columns=['starting_seg', 'ending_seg'])

    if len(atl03file_list) == 0:
        print('WARNING No matching atl03 file found')
    else:
        atl09_segment_ids = readAtlH5(atl09basepath + atl09file,
                                      'profile_1/high_rate/segment_id')
        sc_orient = readAtlH5(atl09basepath + atl09file,
                              'orbit_info/sc_orient')

        if len(sc_orient) == 1 and sc_orient[0] == 0:
            gt = 'gt1l'
        elif len(sc_orient) == 1 and sc_orient[0] == 1:
            gt = 'gt1r'
        else:
            print('WARNING the spacecraft transitions orientation during this granule')

        for atl03file in atl03file_list:
            atl03_segment_ids = readAtlH5(atl03basepath + atl03file,
                                          gt + '/geolocation/segment_id')

            segments_in_both = np.isin(atl09_segment_ids, atl03_segment_ids)
            if np.any(segments_in_both):
                starting_seg = atl09_segment_ids[segments_in_both][0]
                ending_seg = atl09_segment_ids[segments_in_both][-1]
                atl03_overlap_with_atl09.loc[atl03file, 'starting_seg'] = starting_seg
                atl03_overlap_with_atl09.loc[atl03file, 'ending_seg'] = ending_seg

    return atl03_overlap_with_atl09


def atl09_overlap_info(atl09file, atl09basepath, atl08basepath, atl03basepath):
    """Returns dataframe of where the ATL09 file overlaps with ATL03/08 files"""

    matches_03 = get_rgt_matches(input_file_name=atl09file,
                                 filenames_to_search=np.unique(os.listdir(atl03basepath)),
                                 use_full_rgtcs=False)

    if len(matches_03[atl09file]) != 0:
        atl03file_list = matches_03[atl09file]
        overlap_with_atl09_df = check_for_matching_segments(atl09file,
                                                            atl03file_list,
                                                            atl09basepath,
                                                            atl03basepath)

        # Find corresponding atl08 file for each atl03 file
        for atl03file in overlap_with_atl09_df.index:
            matches_08 = get_rgt_matches(input_file_name=atl03file,
                                         filenames_to_search=np.unique(os.listdir(atl08basepath)),
                                         use_full_rgtcs=True)

            if len(matches_08) != 0:
                print('Found ATL08 file for ' + atl03file)
                overlap_with_atl09_df.loc[atl03file, 'atl08file'] = matches_08[atl03file]

            else:
                print('UNABLE TO FIND MATCHING ATL08 FOR ' + atl03file)
    else:
        print('UNABLE TO FIND MATCHING ATL03 FILE')

    overlap_with_atl09_df.reset_index(inplace=True)
    overlap_with_atl09_df.rename({'index': 'atl03file'},
                                 axis='columns',
                                 inplace=True)

    return overlap_with_atl09_df


# %%
if __name__ == "__main__":
    print('Test')
