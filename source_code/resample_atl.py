
"""
This requires h5repack for updating h5 output files

ensure python 3.7 is loaded:
https://www.nccs.nasa.gov/nccs-users/instructional/adapt-instructional/python
This loads 3.7.2:
  $ source /att/opt/other/centos/modules/init/bash
  $ module load stretch/anaconda3

Future updates:
    try to generalize before adding in 09 or up/downsampling options
    add 03 to 08 sample rate support
    add 03 upsampling, if dx < ~1m
    add 08 downsampling, if dx > ~100m
    add variable sampling, if dx == ~dx_08
    add 09 support
    add create.txt file that allows for creating
        datasets in a standard way

    add parallel processing

    30m canopy datasets are relative to WGS84
    08 h5s are split into canpoy/height groups;
        this could be used to determine how to
        interpolate... then maybe placed back
        into groups?? for h5..

    alongtrack is defined where 03 begins, but
    this is not where 08 begins. Normally, this
    is okay, but if there's a large gap
    between the start of 03 and 08, then this
    code will create a 30m binned domain all the
    way from the 03 start that contains zero
    08 data, which just creates a ton of
    dependent-variable nans... could just
    trim df_final at the end for the first
    dependent variable non-nan value?

"""

import matplotlib.pyplot as plt
import numpy as np
import os, sys

if sys.version_info.major < 3:
    print('python version must be > 3.0')
    sys.exit()

import pandas as pd
pd.options.mode.chained_assignment = None


DIR_ICEPY = os.path.normpath('/LIDAR/server/poseidon_files/USERS/jsipps/icesat2_icepy')
sys.path.append(DIR_ICEPY)
import icesatUtils as iu
import icesatPlot as ip
import icesatReader as ir
import region_detect as rd

import h5py as h5

tqdm_found = True
try:
    from tqdm import tqdm
except ImportError:
    tqdm_found = False


INT_MAX = np.iinfo(int).max


def raster1d(x, x0, dx, domain_dataset, df, agg):

    c = np.floor((x - x0) / dx).astype(int)
    c_range = np.arange(min(c), max(c)+1, 1)
    x_dim = x0 + c_range*dx + dx/2.0

    df_group = pd.DataFrame()
    for dataset in df.columns:
        if dataset in agg:
            df_group[dataset] = df[dataset]

    c0 = c_range[0]
    df_group['c'] = c
    df_group = df_group.groupby(['c']).agg(agg)
    col_ind_new = np.array(df_group.index - c0).astype(int)

    x_dim = x_dim[col_ind_new]
    df_group[domain_dataset] = x_dim

    return df_group


def interp_datasets(df, df_ds, domain_dataset, datasets):
    A = np.array(df[domain_dataset])
    A_dim = np.array(df_ds[domain_dataset])
    for key in df.columns:
        if key in datasets:
            y = np.interp(A_dim, A, df[key])
            df_ds[key] = y

    return df_ds


def append_h5(file_h5, df, group_new, OUT_DIR, tag=None, overwrite=False, datasets='all'):

    if type(tag) != type(None):
        if tag[0] != '_':
            tag = '_%s' % tag
    else:
        tag = ''

    file_h5_base = os.path.basename(file_h5)
    i = file_h5_base.index('ATL')
    file_h5_base = file_h5_base[i+6:]
    file_h5_new = os.path.join(OUT_DIR, file_h5_base)
    file_h5_new = '.'.join(file_h5_new.split('.')[:-1]) + tag + '.h5'

    if os.path.exists(file_h5) and not os.path.exists(file_h5_new):
        # first time making the file
        # make new file regardless of overwrite, since
        # there's nothing to overwrite
        #   initialize
        os.system('cp {} {}'.format(file_h5, file_h5_new))

    elif not os.path.exists(file_h5) and not os.path.exists(file_h5_new):
        # first time making new file, and original file doesn't exist
        print('warning: cannot make h5 file %s as the original (%s) does not exist' % (file_h5_new, file_h5))
        return 0
        # can't do anything

    elif os.path.exists(file_h5) and os.path.exists(file_h5_new):
        # new file already made, and original file exists
        if overwrite:
            # overwrite new file
            os.system('cp {} {}'.format(file_h5, file_h5_new))

    elif not os.path.exists(file_h5) and os.path.exists(file_h5_new):
        # new file already made, and original file is gone
        if overwrite:
            # cannot overwrite since original file is gone,
            # can only remove _%dm groups, if anything
            print('warning: cannot overwrite %s as original file (%s) does not exist' % (file_h5_new, file_h5))


    filt_datasets = False
    if datasets != 'all':
        filt_datasets = True
        datasets = list(datasets)

    if not os.path.exists(file_h5_new):
        print('error: %s not found' % file_h5_new)
        return 0

    repack = False
    with h5.File(file_h5_new, 'a') as fp:
        if group_new in fp:
            del fp[group_new]
            repack = True

        fp_df = fp.create_group(group_new)
        dset = {}
        for key in df:
            if filt_datasets:
                if not (key in datasets):
                    continue

            data = np.array(df[key])
            dtype = data.dtype
            if dtype == 'O':
                data = np.array([np.string_(val) for val in data])
                dtype = data.dtype
            n = len(data)
            dset[key] = fp_df.create_dataset(key, (n,), dtype=dtype, compression=1) # [1,9]
            dset[key][:] = data

    if repack:
        file_h5_new_copy = ''.join(os.path.basename(file_h5_new).split('.')[:-1]) + '_copy.h5'
        file_h5_new_copy = os.path.join(os.path.dirname(file_h5_new), file_h5_new_copy)
        os.system('cp {} {}'.format(file_h5_new, file_h5_new_copy))
        os.system('h5repack {} {}'.format(file_h5_new_copy, file_h5_new))
        os.system('rm {}'.format(file_h5_new_copy))


##################################
# user options
dx = 30 # meters, resolution
write_output = False # True
plot_debug = False # True
overwrite_h5 = False # if True, will make new h5s each run

output_types = ['mat', 'pkl', 'h5']

"""
h5 output:
overwrite_h5 is if you'd like to add multiple
resolutions to a single h5 file. Set
overwrite_h5 = False if this is the case, then
run each resolution.

Note that mat/pkl files do not stack similarly;
i.e. they will be remade every time, unless
only 'h5' is specified in output_types.
"""

# DATA_DIR contains 08_datasets.txt
DATA_DIR = os.path.normpath('/LIDAR/server/poseidon_files/USERS/jsipps/scripts_lidar/icesat_analysis/Alaska')

# will have output mat/pkl files (if write_output == True)
OUT_DIR = os.path.normpath('/LIDAR/server/poseidon_files/USERS/jsipps/scripts_lidar/icesat_analysis/Alaska/output')

# 03 files
DIR_03 = os.path.normpath('/bigtex_data/data/release/003/ATL03_r003/Alaska')

# 08 files
DIR_08 = os.path.normpath('/bigtex_data/data/release/003/ATL08_r003/Alaska')


"""
08 upsampling description:
This is broken into a few steps, since there are various
datatypes:
1. floating-point: interpolate floating point to 30m bins
    a.) remove domain datasets, i.e. lon/lat/easting/northing/etc,
        since a new domain is made anyway, for 30m bins
    b.) intp datasets
    c.) if any dataset has only 1 point in a given region,
        attempt to extrapolate data with a linear fit;
        this only works for height datasets, so, for
        a single point:
            i) if 'h_' (or '_h') in dataset name (key), fit a line using
                03 data as x/y points for the slope
                j) if dataset is veg/canopy based, use zero slope;
                    this is b/c canopy is relative to ground
            ii) else, use zero slope

2. integer/str (flags): set int data to 30m bins that overlap 100m bins

"""


##################################
# advanced options/checks (do not modify)
domain_dataset = 'alongtrack'
# dx = 30
dx_08 = 100
# domain_dataset = 'delta_time'
# dx = 0.0021296686652447563
# dx_08 = 0.01422280089504886 # ~100m dt

if type(dx) != int:
    dx = int(np.round(dx))
    print('warning: dx != int, dx -> %d' % dx)

print('dx == %dm' % dx)
if not (1 <= dx <= dx_08):
    print('warning: dx ==', dx)
    if dx < 1:
        print('artifacts may exist in 03 downsampling')
    elif dx > dx_08:
        print('artifacts may exist in 08 upsampling')

output_types_f = []
output_types_allowed = ['mat', 'pkl', 'h5']
for i in range(len(output_types)):
    temp = output_types[i].lower()
    if temp in output_types_allowed:
        output_types_f.append(temp)
    else:
        print('warning: output_type %s not recognized' % temp)

if not write_output:
    print('write_output=False: not writing output')
else:
    print('write_output=True: writing', output_types_f, 'output')

if 'h5' in output_types_f:
    if overwrite_h5:
        # print('warning: h5 in output_types and overwrite_h5=True (not appending data)')
        print('overwrite_h5=True: making new h5 files')
    else:
        # print('warning: h5 in output_types and overwrite_h5=False (appending data)')
        print('overwrite_h5=False: appending to old h5 files')



discrete_datasets = []
# discrete_datasets = ['classification', 'signal_conf_ph']

# intp_datasets:
#   only floating point data that should be intp to 30m bins
#   domain_dataset is automatically removed from intp_datasets
intp_datasets = ['delta_time', 'lat_ph', 'lon_ph', 'easting', 'northing']
if domain_dataset in intp_datasets:
    intp_datasets.remove(domain_dataset)

# define names for new 03 datasets to create at dx resolution
ground_median_dataset = '%dm_ground_median' % dx
canopy_dataset0 = '%dm_canopy_rh_' % dx
canopy_datasets = []
pc = list(np.arange(10,100+5,5)) + [98]
for val in pc:
    canopy_datasets.append(canopy_dataset0 + '%d' % val)

canopy_openness = 'canopy_openness'
canopy_roughness = 'canopy_roughness'
n_te_photons = 'n_te_photons'
n_can_photons = 'n_can_photons'

# read in 08 datasets to keep
datasets_08 = []
with open(os.path.join(DATA_DIR, '08_datasets.txt')) as fp:
    for line in fp:
        line = line.split()
        datasets_08.append(line[0])

# read in 03 files and ensure they match to 08 files
files_03_cmp = []
for file_03_sub in os.listdir(DIR_03):
    if file_03_sub.startswith('ATL03') and file_03_sub.endswith('.h5'):
        files_03_cmp.append(os.path.normpath(DIR_03 + '/' + file_03_sub))

files_03_match, files_08_match = [], []
for f, file_03 in enumerate(files_03_cmp):
    # file_03_sub = file_03.split('/')[-1]
    file_03_sub = os.path.basename(file_03)
    temp = list(file_03_sub)
    temp[4] = '8'
    file_08_sub = ''.join(temp)
    search_str = file_08_sub[:29] # ATL08_20181016010218_02660103
    for file_08_sub in os.listdir(DIR_08):
        files_08_cmp = []
        if file_08_sub.startswith(search_str) and file_08_sub.endswith('.h5'):
            files_08_cmp.append(os.path.normpath(DIR_08 + '/' + file_08_sub))

        if len(files_08_cmp) > 1:
            print('warning: files_08_cmp')
        elif len(files_08_cmp) == 0:
            continue

        file_08 = files_08_cmp[0]
        files_03_match.append(file_03)
        files_08_match.append(file_08)


# f_debug = {191: 'gt3r', 139: 'gt2l', 8: 'gt1r'}
# f_debug = {191: 'gt3r'} #, 139: 'gt2l'}
# f_debug = {232: 'gt2l'}
# f_debug = {244: 'gt1r'}
# Main calc loop
ovrh5 = False

gt_all = ['gt1r', 'gt2r', 'gt3r', 'gt1l', 'gt2l', 'gt3l']
num_files = len(files_03_match)
init = False # used only for debugging
for f in range(num_files):

    # 192 gt3r, 140 gt2l, 9, gt1r
    # if f > 6:
    #   # continue
    #   sys.exit()

    # if init:
    #   break
    # if not (f in f_debug):
    #   continue

    init = True

    # file_03/8 are full-path files
    file_03 = files_03_match[f]
    file_08 = files_08_match[f]

    # output progress to prompt
    print('%d/%d' % (f+1, num_files))
    print(file_03)
    print(file_08)

    if overwrite_h5:
        ovrh5 = True
        # per ground track

    # if 1:
    #     gt = f_debug[f]
    for gt in gt_all:
        # if gt != 'gt1r':
        #   continue

        print(gt)

        try:
            # reads and classifies photons, if valid
            atl03 = ir.get_atl03_struct(file_03, gt, file_08)
        except (KeyError, UnboundLocalError) as e:
            print(file_03)
            print('warning:', e)
            # usually when "heights" dataset isn't found, or atl file is corrupted
            continue

        if not atl03.dataIsMapped:
            print('warning: data not mapped')
            # cannot do analysis without classified photons
            continue

        df_03 = atl03.df

        df1 = df_03[df_03.classification == 1] # ground
        df2 = df_03[df_03.classification == 2] # low veg
        df3 = df_03[df_03.classification == 3] # high veg

        df_ground = df_03[df_03.classification == 1] # ground
        df_canopy = df_03[(df_03.classification == 2) | (df_03.classification == 3)] # low+high veg == canopy
        if df_ground.empty and df_canopy.empty:
            print('warning: df_ground and df_canopy empty')
            continue

        x0 = min(df_03[domain_dataset]) # common origin for all data; for rastering in 1d

        # this interpolates all "domain" datasets to 30m bins
        #   alongtrack, crosstrack, lon/lat, E/N, delta_time/time, etc
        df_03_ds = raster1d(df_03[domain_dataset], x0, dx, domain_dataset, df_03, {domain_dataset: np.mean})
        df_03_ds = interp_datasets(df_03, df_03_ds, domain_dataset, intp_datasets)

        # additionally, add delta_time start/end to the bin
        A = df_03[domain_dataset]
        A_dim = df_03_ds[domain_dataset]
        delta_time_start = np.interp(A_dim - dx/2.0, A, df_03['delta_time'])
        delta_time_end = np.interp(A_dim + dx/2.0, A, df_03['delta_time'])

        df_03_ds['delta_time_start'] = delta_time_start
        df_03_ds['delta_time_end'] = delta_time_end


        # Must have a different raster for every classification subset, since
        # the domain is changing. However, rasters could be trimmed down within
        # a single classification subset, such as with df_canopy.

        if len(discrete_datasets) > 0:
            agg_discrete = {}
            for dataset in discrete_datasets:
                agg_discrete[dataset] = iu.mode
            df_03_discrete = raster1d(df_03[domain_dataset], x0, dx, domain_dataset, df_03, agg_discrete)
            df_03_ds = df_03_ds.join(df_03_discrete.drop(columns=[domain_dataset]))

        if not df2.empty:
            # canopy_openness
            df2_co = raster1d(df2[domain_dataset], x0, dx, domain_dataset, df2, {'h_ph': np.std})
            df2_co = df2_co.rename(columns={'h_ph': canopy_openness})
            df_03_ds = df_03_ds.join(df2_co.drop(columns=[domain_dataset]))

        if not df3.empty:
            # canopy_roughness
            df3_co = raster1d(df3[domain_dataset], x0, dx, domain_dataset, df3, {'h_ph': np.std})
            df3_co = df3_co.rename(columns={'h_ph': canopy_roughness})
            df_03_ds = df_03_ds.join(df3_co.drop(columns=[domain_dataset]))

        if not df_ground.empty:
            # 30m_ground_median
            df_ground_ds = raster1d(df_ground[domain_dataset], x0, dx, domain_dataset, df_ground, {'h_ph': np.median})
            df_ground_ds = df_ground_ds.rename(columns={'h_ph': ground_median_dataset})
            df_03_ds = df_03_ds.join(df_ground_ds.drop(columns=[domain_dataset]))

            # n_te_photons
            df_ground_ds = raster1d(df_ground[domain_dataset], x0, dx, domain_dataset, df_ground, {'h_ph': np.size})
            df_ground_ds = df_ground_ds.rename(columns={'h_ph': n_te_photons})
            df_03_ds = df_03_ds.join(df_ground_ds.drop(columns=[domain_dataset]))

        if not df_canopy.empty:
            # 30m_canopy_rh_##
            #   could be made into a single raster?.. maybe trimmed
            for k, c_dataset in enumerate(canopy_datasets):
                agg_canopy = {'h_ph': lambda x: np.nanpercentile(x, pc[k])}
                df_canopy_ds = raster1d(df_canopy[domain_dataset], x0, dx, domain_dataset, df_canopy, agg_canopy)
                df_canopy_ds = df_canopy_ds.rename(columns={'h_ph': c_dataset})
                df_03_ds = df_03_ds.join(df_canopy_ds.drop(columns=[domain_dataset]))

            # n_can_photons
            df_canopy_ds = raster1d(df_canopy[domain_dataset], x0, dx, domain_dataset, df_canopy, {'h_ph': np.size})
            df_canopy_ds = df_canopy_ds.rename(columns={'h_ph': n_can_photons})
            df_03_ds = df_03_ds.join(df_canopy_ds.drop(columns=[domain_dataset]))


        df_03_ds = df_03_ds.reset_index(drop=True)
        entries = [ground_median_dataset, canopy_openness, canopy_roughness, n_te_photons, n_can_photons] + canopy_datasets
        for dataset in entries:
            if not dataset in df_03_ds:
                df_03_ds[dataset] = np.nan


        if plot_debug:
            # check df_03_ds
            fig, ax = ip.make_fig()
            ax.plot(df_ground[domain_dataset], df_ground['h_ph'], '.', color=ip.BROWN, label='df_ground h_ph')
            ax.plot(df_canopy[domain_dataset], df_canopy['h_ph'], '.', color='C2', label='df_canopy h_ph')
            ax.plot(df_03_ds[domain_dataset], df_03_ds[ground_median_dataset], '.', color='r', label=ground_median_dataset)
            ax.plot(df_03_ds[domain_dataset], df_03_ds['%dm_canopy_rh_98' % dx], '.', color='C1', label='%dm_canopy_rh_98' % dx)
            ax.legend(fontsize=9)
            ax.set_title('check df_03_ds vs. df_ground/canopy')
            fig.show()


        # load in 08, match to 03 ttg
        try:
            atl08 = ir.get_atl08_struct(file_08, gt, atl03)
        except (KeyError, UnboundLocalError) as e:
            print(file_08)
            print('warning:', e)
            continue

        df_08 = atl08.df
        df_08 = df_08.sort_values(by=[domain_dataset]) # rarely domain can reverse

        # assign 1e30 error to np.nan for
        # consistency (helps with np.interp)
        for key in df_08.columns:
            dtype = str(df_08[key].dtype)
            # if 'int' in dtype:
            #   continue
            if not ('float' in dtype):
                continue
            b = abs(df_08[key]) > 1e20
            df_08[key][b] = np.nan

        # double-check that all floating values
        # have been assigned to nans properly
        for key in df_08.columns:
            dtype = str(df_08[key].dtype)
            # if 'int' in dtype:
            #   continue
            if not ('float' in dtype):
                continue
            vals = abs(np.array(df_08[key]))
            b = ~np.isnan(vals)
            if (vals[b] > 1e20).any():
                print(key)
                print('warning: vals > 1e20')



        x_03_ds = np.array(df_03_ds[domain_dataset])
        y_03_ds = np.array(df_03_ds[ground_median_dataset])
        # ground_03_ds = np.array(df_03_ds['atl03_ground_median'])
        # canopy_03_ds = np.array(df_03_ds['atl03_canopy_max98'])
        j_03_ds = np.array(df_03_ds.index).astype(int)

        x_08 = np.array(df_08[domain_dataset])

        # interpolate all 08 datasets to dx resolution
        #   done by splitting up the data into regions; 
        #   where dx res overlaps with 08 is where 
        #   interpolation is possible
        key_int = [] # anything that isn't a float, actually inc str
        df_08_us = pd.DataFrame()
        for key in df_08.columns:

            # save flag/integer data for later treatment
            # by binning
            if not (key in datasets_08):
                # user choice
                continue

            dtype = str(df_08[key].dtype)
            # if 'int' in dtype:
            #   key_int.append(key)
            #   continue
            if not ('float' in dtype):
                key_int.append(key)
                continue

            # remove domain datasets as 30m bins have
            # intp domain anyway
            if key in df_03_ds.columns:
                continue
            if key == 'longitude' or key == 'latitude':
                continue
            if 'delta_time' in key:
                continue

            y_08 = np.array(df_08[key])
            x_08_y = x_08[~np.isnan(y_08)]
            bins_08 = []
            if len(x_08_y) > 0:
                bins_08 = rd.bin_reg(x_08_y, dx_08)
                regions_08 = rd.combine_region(bins_08)

            intp_total = np.full(x_03_ds.shape, np.nan)
            for reg in regions_08:
                b_08_reg = (x_08 >= reg[0]) & (x_08 <= reg[1])
                x_08_reg = x_08[b_08_reg]
                y_08_reg = y_08[b_08_reg]

                b_03_reg = (x_03_ds >= reg[0]) & (x_03_ds <= reg[1])
                x_03_reg = x_03_ds[b_03_reg]
                j_03_reg = j_03_ds[b_03_reg]

                if len(x_08_reg) > 1:
                    # multiple points to interpolate thru
                    intp = np.interp(x_03_reg, x_08_reg, y_08_reg)

                elif len(x_08_reg) == 1:
                    # only 1 pt, so intp not possible; must extrapolate
                    y_03_reg = y_03_ds[b_03_reg]

                    b_nan = np.isnan(x_03_reg) | np.isnan(y_03_reg)
                    x_03_reg, y_03_reg, j_03_reg = x_03_reg[~b_nan], y_03_reg[~b_nan], j_03_reg[~b_nan]

                    if len(x_03_reg) != len(y_03_reg):
                        print('error: length')
                        iu.pause()

                    if len(x_03_reg) <= 1:
                        # only 1 point of 08 and 1 or 0 points of 03,
                        # no way to interpolate
                        continue

                    x_regr = np.copy(np.reshape(x_03_reg, (len(x_03_reg),1)))
                    y_regr = np.copy(np.reshape(y_03_reg, (len(y_03_reg),1)))
                    x_regr -= x_regr[0] # for cond(H), will mess up intercept

                    H = np.hstack((x_regr, np.ones(x_regr.shape)))
                    beta = np.matmul(np.linalg.pinv(H), y_regr)
                    slope, _ = float(beta[0]), float(beta[1])

                    x0 = x_08_reg[0] # only one point
                    y0 = y_08_reg[0]

                    if (('h_' in key) or ('_h' in key)) and not ('canopy' in key):
                        """
                        slope is in units of height / domain
                        The only datasets that can be adjusted by a non-zero
                        slope are those that correspond to spatial data, i.e.
                        height metrics;
                        Further, 08 canopy datasets are relative to ground,
                        so their slope is zero b/c the ground creates the
                        sloping effect
                        """
                        intp = slope*(x_03_reg - x0) + y0

                    else:
                        # non-spatial floating-point datasets, 
                        # i.e. asr, snr, toc_roughness, etc
                        #   give these constant slope
                        intp = np.full(x_03_reg.shape, y0)

                else:
                    continue

                intp_total[j_03_reg] = intp

            # key += '_08'
            df_08_us[key] = intp_total
            
        df_08_us.reset_index(drop=True)

        # join 03 downsample with 08 upsample
        df_final = pd.concat([df_03_ds, df_08_us], axis=1)


        """
        Indexing for flags:
            Some discrete data cannot be interpolated;
            it has to be set to the value of the lower
            sampling rate. That is, for a given 08
            100m bin, the 30m bins that overlap must be
            given the flag number that the 08 bin has.
        """     
        bins_08 = rd.bin_reg(x_08, dx_08)

        try:
            # fast, robust
            #   finds region overlaps using a hash

            bins_03 = rd.bin_reg(x_03_ds, dx)
            n_03_ds = len(x_03_ds)

            # from scipy.interpolate import CubicSpline, interp1d
            from scipy.interpolate import InterpolatedUnivariateSpline as scipy_intp

            def calc_err_cs(cs, N, N_s):

                err_cs = []
                n = len(N)
                for j in range(n):
                    j_test = int(cs(N_s[j]))
                    if j_test < 0:
                        j_test = 0
                    elif j_test >= n:
                        j_test = n-1
                    val = N[j_test]
                    err = val - N_s[j]
                    err_cs.append(err)

                return np.array(err_cs)

            def calc_cs(N, dx_mu, calc_err=False):
                index = np.array([k for k in range(len(N))])
                N_s, index_s = iu.sort_first([N, index])
                index_s = index_s.astype(int)
                # Fs_mu = 1.0 / dx_mu
                # N_ds, _, (index_ds) = cm.downsample(Fs_mu, N_s, index_s)
                N_ds, index_ds = np.copy(N), np.copy(index).astype(int)

                # cs = CubicSpline(N_ds, index_ds)
                cs = scipy_intp(N_ds, index_ds, k=1) # linear intp
                err_cs = np.array([])
                if calc_err:
                    err_cs = calc_err_cs(cs, N, N_s)
                return cs, err_cs

            cs, err_cs = calc_cs(x_08, dx_08, calc_err=True)
            # j1_est = int(cs1(N2_08[j2]))

            if (abs(err_cs) > 0.0).any():
                print('warning: err_cs')

            if plot_debug:
                x = np.arange(min(x_08), max(x_08))
                index_08 = np.array([k for k in range(len(x_08))])
                fig, ax = ip.make_fig()
                ax.plot(x_08, index_08, '.', label='index_08')
                ax.plot(x_03_ds, cs(x_03_ds), '.', label='cs(x_03_ds)')
                ax.plot(x, cs(x), label='cs(x)')
                ax.legend(fontsize=9)
                ax.set_title('hash')
                fig.show()

            i_index = np.full(n_03_ds, -1).astype(int)

            n_08 = len(bins_08)
            di = 10 # large safety factor
            for j in range(n_03_ds):
                bin_03 = bins_03[j]

                mu = np.mean(bin_03)
                i_est = int(np.round(cs(mu)))

                lower = i_est-di
                if lower < 0:
                    lower = 0
                elif lower >= n_08:
                    continue

                upper = i_est+di+1
                if upper <= 0:
                    continue
                elif upper > n_08:
                    upper = n_08

                p_08 = []
                for i in range(lower, upper):
                    bin_08 = bins_08[i]
                    reg_inner = rd.inner_region(bin_08, bin_03)
                    p = 0.0
                    if reg_inner != []:
                        # percentage of overlap
                        p = (reg_inner[1] - reg_inner[0]) / (bin_03[1] - bin_03[0])
                    # else:
                    #   # no overlap
                    #   # p = 0.0
                    #   pass
                    p_08.append([p, i])

                if len(p_08) == 0:
                    continue

                p_08 = sorted(p_08, reverse=True)
                p_max, i_max = p_08[0]
                if p_max < 0.5:
                    continue

                i_index[j] = i_max

        except KeyboardInterrupt:
            raise

        except:

            print('warning: hashing failed, trying alternate method')
            # slow but robust
            #   non-hash version of above

            bins_03 = rd.bin_reg(x_03_ds, dx)
            n_03_ds = len(x_03_ds)

            i_index = np.full(n_03_ds, -1).astype(int)
            iterator = range(n_03_ds)
            if tqdm_found:
                iterator = tqdm(iterator)
            # for j in tqdm(range(n_03_ds)):
            for j in iterator:
                bin_03 = bins_03[j]
                p_08 = []
                for i, bin_08 in enumerate(bins_08):
                    reg_inner = rd.inner_region(bin_08, bin_03)
                    p = 0.0
                    if reg_inner != []:
                        # percentage of overlap
                        p = (reg_inner[1] - reg_inner[0]) / (bin_03[1] - bin_03[0])
                    else:
                        # no overlap
                        # p = 0.0
                        pass

                    p_08.append([p, i])

                if len(p_08) == 0:
                    continue

                p_08 = sorted(p_08, reverse=True)
                p_max, i_max = p_08[0]
                if p_max < 0.5:
                    continue

                i_index[j] = i_max


        if plot_debug:
            b = i_index == -1

            fig, ax = ip.make_fig()
            ax.plot(df_03_ds[domain_dataset], df_03_ds[ground_median_dataset], '.', color='C0', label=ground_median_dataset)
            ax.plot(df_03_ds[domain_dataset], df_03_ds['%dm_canopy_rh_98' % dx], '.', color='C1', label='%dm_canopy_rh_98' % dx)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            dy = 10.0
            for i, reg in enumerate(bins_08):
                ax.plot([reg[0], reg[0]], ylim, 'g--')
                ax.plot([reg[1], reg[1]], ylim, 'r')

            ax.plot(df_03_ds[domain_dataset][b], df_03_ds[ground_median_dataset][b], 'r.', label='%s (invalid)' % ground_median_dataset)
            ax.plot(df_03_ds[domain_dataset][b], df_03_ds['%dm_canopy_rh_98' % dx][b], 'r.', label='%dm_canopy_rh_98 (invalid)' % dx)

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.legend(fontsize=9)
            ax.set_title('check i_index')
            fig.show()
            

        df_final['flag_index'] = i_index

        # upsample flags via binning/flag index
        for key in key_int:
            vec = []
            for j, i in enumerate(df_final['flag_index']):
                if i == -1:
                    if key != 'beam_type':
                        vec.append(INT_MAX)
                    else:
                        vec.append('unknown')
                else:
                    vec.append(df_08[key][i])

            if key != 'beam_type':
                df_final[key] = np.array(vec).astype(int)

            else:
                df_final[key] = np.array(vec)

        df_final = df_final.drop(columns=['flag_index'])
        df_final = df_final.reset_index(drop=True)


        year, doy = iu.get_h5_meta(file_08, meta='date', rtn_doy=True)
        try:
            sc_orient = iu.get_sc_orient(file_08, np.array(df_final['delta_time']))
            beam_number, beam_type = iu.get_beam_info(sc_orient, gt)
        except KeyError:
            print('error: could not find [%s,%s] in 08, using 08 attrs..' % (datasets[0], datasets[1]))
            sc_orient = -1
            beam_number = -1
            beam_type = 'unknown'
            with h5.File(file_08, 'r') as fp:
                fp_a = fp[gt].attrs
                if 'sc_orientation' in fp_a and 'atlas_spot_number' in fp_a and 'atlas_beam_type' in fp_a:
                    # print(fp_a['sc_orientation'])
                    sc_orient_str = (fp_a['sc_orientation']).decode().lower()
                    if sc_orient_str == 'forward':
                        sc_orient = 1
                    elif sc_orient_str == 'backward':
                        sc_orient = 0

                    beam_number = (fp_a['atlas_spot_number']).decode()
                    try:
                        beam_number = int(beam_number)
                    except ValueError:
                        print(file_08, gt)
                        print('error: beam_number (1)')
                        print(beam_number)
                        # iu.pause()

                    beam_type = (fp_a['atlas_beam_type']).decode()


            if not (sc_orient == 0 or sc_orient == 1):
                print('error: sc_orient')
                print(file_08, gt)
                print(sc_orient_str, sc_orient)
                # iu.pause()

            if not (beam_number in [1,2,3,4,5,6]):
                print(file_08, gt)
                print('error: beam_number (2)')
                # iu.pause()

            if not (beam_type in ['strong', 'weak']):
                print(file_08, gt)
                print('error: beam_type')
                # iu.pause()


        df_final['year'] = int(year)
        df_final['doy'] = int(doy)
        df_final['sc_orient'] = sc_orient
        df_final['beam_number'] = beam_number
        df_final['beam_type'] = beam_type


        # output dx res dataframe via pkl and mat files
        OUT_DIR_PKL = os.path.join(OUT_DIR, '%dm'%dx, 'pkl')
        OUT_DIR_MAT = os.path.join(OUT_DIR, '%dm'%dx, 'mat')
        OUT_DIR_H5 = os.path.join(OUT_DIR, 'h5')

        file_08_sub = os.path.basename(file_08)
        file_pkl = file_08_sub[6:-3] + '_%s.pkl' % gt
        file_mat = file_08_sub[6:-3] + '_%s.mat' % gt
        file_h5 = file_08_sub[6:-3] + '.h5'

        if write_output:
            def makedir(d):
                if not os.path.exists(d):
                    os.makedirs(d)

            if 'pkl' in output_types_f:
                makedir(OUT_DIR_PKL)
                ir.write_pickle(df_final, os.path.join(OUT_DIR_PKL, file_pkl))

            if 'mat' in output_types_f:
                makedir(OUT_DIR_MAT)
                ir.convert_df_to_mat(df_final, os.path.join(OUT_DIR_MAT, file_mat))

            if 'h5' in output_types_f:
                makedir(OUT_DIR_H5)
                append_h5(file_08, df_final, '%s_%dm/land_segments' % (gt, dx), OUT_DIR_H5, overwrite=ovrh5, datasets='all')
                ovrh5 = False

        # break
    # break


# DIR = '/LIDAR/server/poseidon_files/USERS/jsipps/scripts_lidar/icesat_analysis/Alaska/output/h5'
# fn = DIR + '/20181016010218_02660103_003_01.h5'

# fp = h5.File(fn, 'r')
# x = np.array(fp['gt1r_30m/alongtrack'])
# y = np.array(fp['gt1r_30m/n_ca_photons'])
# z = np.array(fp['gt1r_30m/h_te_best_fit'])


# import file_search as fs
# DIR = '/LIDAR/server/poseidon_files/USERS/jsipps/scripts_lidar/icesat_analysis/Alaska/output/h5'
# for fn0 in sorted(os.listdir(DIR)):
#     print(fn0)
#     fn = DIR + '/' + fn0
#     fs.show_h5(fn)
#     iu.pause()
