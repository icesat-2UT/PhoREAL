
# ensure python 3.7 is loaded:
# https://www.nccs.nasa.gov/nccs-users/instructional/adapt-instructional/python
# This loads 3.7.2:
#   $ source /att/opt/other/centos/modules/init/bash
#   $ module load stretch/anaconda3

import matplotlib.pyplot as plt
import numpy as np
import os, sys

if sys.version_info.major < 3:
    print('python version must be > 3.0')
    sys.exit()

import pandas as pd
pd.options.mode.chained_assignment = None


DIR_ICEPY = '.'
sys.path.append(DIR_ICEPY)
import icesatUtils as iu
import icesatPlot as ip
import icesatReader as ir
import region_detect as rd

import h5py as h5


# def nancount(vec):
#   return np.isnan(vec).sum()
    
def convert_df_to_mat(df,outfilename):
    from scipy import io
    comps =  outfilename.split('.')
    if comps[-1] != 'mat':
        outfilename = outfilename + ".mat"
    io.savemat(outfilename, {'struct':df.to_dict("list")})


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


def bin_reg(A, dx):

    if (np.diff(A) < 0.0).any():
        print('warning: x decreasing')

    if dx <= 0:
        print('warning: dx < 0')

    A = np.append(min(A) - dx, A)
    A = np.append(A, max(A) + dx)

    dx_lim = 1.5*dx
    bin_dim = []
    for j in range(1,len(A)-1):
        da1 = A[j] - A[j-1]
        da2 = A[j+1] - A[j]
        if da1 > dx_lim:
            da1 = dx_lim
        if da2 > dx_lim:
            da2 = dx_lim
        bin_dim.append([A[j] - da1/2.0, A[j] + da2/2.0])

    return bin_dim


def check_bin_overlap(bins, tag, loading_bar=False):

    print('checking bin overlap for %s..' % tag)
    iterator = range(len(bins))
    if loading_bar:
        from tqdm import tqdm
        iterator = tqdm(range(len(bins)))

    # for i, bin1 in enumerate(bins):
    for i in iterator:
        bin1 = bins[i]
        for j, bin2 in enumerate(bins):
            if i == j:
                continue
            # bin1, bin2 = list(np.round(bin1,8)), list(np.round(bin2,8))
            reg = rd.inner_region(bin1, bin2, equal=False)
            if reg != []:
                if reg[1]-reg[0] < 1e-6:
                    continue
                print('warning: %s overlap' % tag)
                print(i, bin1)
                print(j, bin2)
                print(reg)
                iu.pause()

##################################
# user options
dx = 30 # meters, resolution
write_output = False # True
plot_debug = False # True
debug = False # True

# DATA_DIR contains 08_datasets.txt
DATA_DIR = ''

# will have output mat/pkl files (if write_output == True)
OUT_DIR = ''

# 03 files
DIR_03 = ''

# 08 files
DIR_08 = ''


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

2. integer (flags): set int data to 30m bins that overlap 100m bins

"""


##################################
# advanced options (do not modify)
domain_dataset = 'alongtrack'
# dx = 30
dx_08 = 100
# domain_dataset = 'delta_time'
# dx = 0.0021296686652447563
# dx_08 = 0.01422280089504886 # ~100m dt

if not (1 <= dx <= dx_08):
    print('warning: dx ==', dx)
    if dx < 1:
        print('artifacts may exist in 03 downsampling')
    elif dx > dx_08:
        print('artifacts may exist in 08 upsampling')

if type(dx) != int:
    dx = int(np.round(dx))
    print('warning: dx != int, dx -> %d' % dx)

if not write_output:
    print('warning: not writing output')

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


# for f in range(len(files_03_match)):
#   print(files_03_match[f])
#   print(files_08_match[f])
#   cm.pause()


# files_03_cmp = cm0.search(DIR_03, ['ATL03'], ext='h5')
# # files_08_cmp = cm0.search(DIR_08, ['ATL08'], ext='h5')

# files_03_match = []
# files_08_match = []
# for f, file_03 in enumerate(files_03_cmp):
#   files_08_cmp = cm0.find_match(file_03, 'ATL08', DIR_08, debug=0)
#   if len(files_08_cmp) > 1:
#       print('warning: files_08_cmp')
#   elif len(files_08_cmp) == 0:
#       continue

#   file_08 = files_08_cmp[0]
#   files_03_match.append(file_03)
#   files_08_match.append(file_08)




# f_debug = {191: 'gt3r', 139: 'gt2l', 8: 'gt1r'}
# f_debug = {191: 'gt3r'} #, 139: 'gt2l'}
# f_debug = {15: 'gt2l'}
# Main calc loop
gt_all = ['gt1r', 'gt2r', 'gt3r', 'gt1l', 'gt2l', 'gt3l']
num_files = len(files_03_match)
init = False # used only for debugging
for f in range(num_files):

    # 192 gt3r, 140 gt2l, 9, gt1r
    # if f <= 8:
    #   continue

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

            # fig, ax = ip.make_fig()
            # ax.plot(df_ground[domain_dataset], df_ground['h_ph'], '.', color=ip.BROWN, label='df_ground h_ph')
            # ax.plot(df_canopy[domain_dataset], df_canopy['h_ph'], '.', color='C2', label='df_canopy h_ph')
            # ax.plot(df_03_ds[domain_dataset], df_03_ds[ground_median_dataset], '.', color='r', label=ground_median_dataset)
            # # for c_dataset in canopy_datasets:
            # #     ax.plot(df_03_ds[domain_dataset], df_03_ds[c_dataset], '.')
            # ax.legend(fontsize=9)
            # ax.set_title('check df_03_ds vs. df_ground/canopy')
            # fig.show()


        # load in 08, match to 03 ttg
        try:
            atl08 = ir.get_atl08_struct(file_08, gt, atl03)
        except (KeyError, UnboundLocalError) as e:
            print(file_08)
            print('warning:', e)
            continue

        df_08 = atl08.df
        df_08 = df_08.sort_values(by=[domain_dataset]) # rarely domain can reverse

        # add year/doy, sc_orient, beam_number/type to 08 dataframe
        year, doy = iu.get_h5_meta(file_08, meta='date', rtn_doy=True)
        sc_orient = -1
        beam_number = -1
        beam_type = 'u'
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
                    print('warning: beam_number (1)')
                    print(beam_number)
                    iu.pause()

                beam_type = (fp_a['atlas_beam_type']).decode()


        if not (sc_orient == 0 or sc_orient == 1):
            print('warning: sc_orient')
            print(file_08)
            print(sc_orient_str, sc_orient)
            iu.pause()

        if not (beam_number in [1,2,3,4,5,6]):
            print('warning: beam_number (2)')
            iu.pause()

        if not (beam_type in ['strong', 'weak']):
            print('warning: beam_type')
            iu.pause()

        df_08['year'] = int(year)
        df_08['doy'] = int(doy)
        df_08['sc_orient'] = sc_orient
        df_08['beam_number'] = beam_number
        df_08['beam_type'] = beam_type


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
                bins_08 = bin_reg(x_08_y, dx_08)
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
                        print('warning: length')
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
                        slope is in units of height / alongtrack
                        The only datasets that can be adjusted by a non-zero
                        slope are those that correspond to spatial data, i.e.
                        height metrics;
                        Further, 08 canopy datasets are relative to ground,
                        so their slope is zero b/c the ground creates the
                        sloping effect
                        """

                        intp = slope*(x_03_reg - x0) + y0

                        # if key == 'h_canopy':
                        #   print('')
                        #   print(key)
                        #   print(x0, y0)
                        #   print('')

                        #   fig, ax = cm.make_fig()
                        #   ax.plot(x_03_reg, y_03_reg, '.')
                        #   ax.plot(x_03_reg, np.mean(y_03_reg) + intp)
                        #   fig.show()

                        #   cm.pause()
                        #   plt.close(fig)

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
        bins_08 = bin_reg(x_08, dx_08)
        if debug:
            check_bin_overlap(bins_08, 'bins_08', loading_bar=True)

        # error check on bins_08 having no overlaps
        #   if bins overlap, this can get the flag
        #   index off by 1 or more 08 bins, causing
        #   all flags past that bin to be removed

        option = 2
        if option == 1:

            n_03_ds = len(x_03_ds)
            # np.random.shuffle(x_03_ds)

            # ensure x_03_ds is always increasing for binning
            index_s = np.argsort(x_03_ds)
            index_s_rev = np.argsort(index_s)
            x_03_s = x_03_ds[index_s]

            if ((x_03_s[index_s_rev] - x_03_ds) != 0.0).any():
                # check that x_03_s is increasing, and sorting
                # is reversible
                print('warning: index_s')

            """
            loop over dx res bins, determine which 08 bin
            overlays, then assign i index, which will be
            used later to determine which flag value goes
            to which dx res bin

            There's something wrong with this..
            I know for sure, if an 08 bin exists with no
            03 point in it, or an 08 bin is zero in length,
            then it'll mess up.
            """
            i_index_s = np.full(n_03_ds, -1).astype(int)
            in_reg = 0 # in region
            i = 0 # 08 bin index
            for j in range(n_03_ds):
                if bins_08[i][0] <= x_03_s[j] <= bins_08[i][1]:
                    i_index_s[index_s_rev[j]] = i
                    in_reg = 1
                else:
                    if in_reg:
                        in_reg = 0
                        if i < len(bins_08)-1:
                            i += 1
                            if bins_08[i][0] <= x_03_s[j] <= bins_08[i][1]:
                                i_index_s[index_s_rev[j]] = i
                        else:
                            break

            i_index = i_index_s[index_s_rev]


        elif option == 2:
            try:
                # fast, robust
                #   finds region overlaps using a hash
                #   if this fails for any reason, try option 3
                #   option 3 is a non-hash version

                bins_03 = bin_reg(x_03_ds, dx)
                n_03_ds = len(x_03_ds)

                if debug:
                    check_bin_overlap(bins_03, 'bins_03 (1) (option == %d)' % option, loading_bar=True)

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

                # fig, ax = ip.make_fig()
                # ax.plot(index_08, '.')
                # ax.plot(cs(x_08))
                # fig.show()
                if (abs(err_cs) > 0.0).any():
                    print('warning: err_cs')

                if plot_debug:
                    x = np.arange(min(x_08), max(x_08))
                    index_08 = np.array([k for k in range(len(x_08))])
                    fig, ax = ip.make_fig()
                    ax.plot(x_08, index_08, '.', label='index_08')
                    # ax.plot(x_03_ds, cs(x_03_ds), '.', label='cs(x_03_ds)')
                    # ax.plot(x, cs(x), label='cs(x)')
                    ax.legend(fontsize=9)
                    ax.set_title('hash')
                    fig.show()

                i_index = np.full(n_03_ds, -1).astype(int)

                n_08 = len(bins_08)
                di = 10 # large safety factor
                # for j, bin_03 in enumerate(bins_03):
                # for j in tqdm(range(n_03_ds)):
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
                    # for i in i_range:
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
                    # if p_max == 0.0:
                    if p_max < 0.5:
                        # equivalent criteria to option 1
                        continue

                    # if bin_03[0] > 225544:
                    #   print(bin_03)
                    #   print(np.mean(bin_03))

                    #   print(bin_08)
                    #   print(np.mean(bin_08))

                    #   print(p_max, j, i_max)
                    #   n0 = min([3,len(p_08)])
                    #   print(p_08[:n0])
                    #   pause()

                    # j matches with i_max,
                    # use 08 flag at i_max

                    i_index[j] = i_max

            except KeyboardInterrupt:
                raise

            except:

                # slow but robust
                #   non-hash version of option 2

                bins_03 = bin_reg(x_03_ds, dx)
                n_03_ds = len(x_03_ds)

                if debug:
                    check_bin_overlap(bins_03, 'bins_03 (2) (option == %d)' % option, loading_bar=True)

                i_index = np.full(n_03_ds, -1).astype(int)
                # for j, bin_03 in enumerate(bins_03):
                from tqdm import tqdm
                for j in tqdm(range(n_03_ds)):
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
                    # if p_max == 0.0:
                    if p_max < 0.5:
                        # equivalent criteria to option 1
                        continue

                    # if bin_03[0] > 225544:
                    #   print(bin_03)
                    #   print(np.mean(bin_03))

                    #   print(bin_08)
                    #   print(np.mean(bin_08))

                    #   print(p_max, j, i_max)
                    #   n0 = min([3,len(p_08)])
                    #   print(p_08[:n0])
                    #   pause()

                    # j matches with i_max,
                    # use 08 flag at i_max

                    i_index[j] = i_max

                    # pause()

        # elif option == 4:
        #   bins_03 = bin_reg(x_03_ds, dx)

        #   regions_03 = rd.combine_region(bins_03)
        #   regions_08 = rd.combine_region(bins_08)

        #   regions_overlap = []
        #   for reg2 in regions_03:
        #       for reg1 in regions_08:
        #           reg = rd.inner_region(reg2, reg1)
        #           if reg != []:
        #               regions_overlap.append(reg)

        #   x_03_ds_filt, jr = rd.filt_reg_x(x_03_ds, regions_overlap)
        #   x_08_filt, ir = rd.filt_reg_x(x_08, regions_overlap)


        if plot_debug:
            b = i_index == -1

            fig, ax = ip.make_fig()
            ax.plot(df_03_ds[domain_dataset], df_03_ds[ground_median_dataset], '.', color='C0', label=ground_median_dataset)
            ax.plot(df_03_ds[domain_dataset], df_03_ds['%dm_canopy_rh_98' % dx], '.', color='C1', label='%dm_canopy_rh_98' % dx)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            dy = 10.0
            # h_08_est = np.array(df_08['h_te_median'])
            for i, reg in enumerate(bins_08):
                # mu = h_08_est[i]
                # ax.plot([reg[0], reg[0]], [mu-dy, mu+dy], 'g--')
                # ax.plot([reg[1], reg[1]], [mu-dy, mu+dy], 'r--')
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
        INT_MAX = np.iinfo(int).max
        for key in key_int:
            vec = []
            for j, i in enumerate(df_final['flag_index']):
                if i == -1:
                    # vec.append(np.nan)
                    if key != 'beam_type':
                        vec.append(INT_MAX)
                    else:
                        vec.append('unknown')
                else:
                    vec.append(df_08[key][i])

            if key != 'beam_type':
                # key += '_08'
                df_final[key] = np.array(vec).astype(int)

            else:
                # key += '_08'
                df_final[key] = np.array(vec)

        df_final = df_final.drop(columns=['flag_index'])
        df_final = df_final.reset_index(drop=True)


        # output dx res dataframe via pkl and mat files
        OUT_DIR_PKL = os.path.normpath(OUT_DIR + '/%dm/pkl' % dx)
        OUT_DIR_MAT = os.path.normpath(OUT_DIR + '/%dm/mat' % dx)
        for d in [OUT_DIR_PKL, OUT_DIR_MAT]:
            if not os.path.exists(d):
                os.makedirs(d)

        # file_08_sub = file_08.split('/')[-1]
        file_08_sub = os.path.basename(file_08)
        file_pkl = file_08_sub[6:-3] + '_%s.pkl' % gt
        file_mat = file_08_sub[6:-3] + '_%s.mat' % gt

        if write_output:
            ir.write_pickle(df_final, os.path.normpath(OUT_DIR_PKL + '/' + file_pkl))
            convert_df_to_mat(df_final, os.path.normpath(OUT_DIR_MAT + '/' + file_mat))

        # break
    # break
