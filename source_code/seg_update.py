
import numpy as np
import os, sys

import pandas as pd

import h5py as h5
from tqdm import tqdm
import icesatReader

import time


DIR = '/LIDAR/server/poseidon_files/USERS/eric/1_experiment/finland_analysis3'
# 	   /LIDAR/server/USERS/eric/1_experiment/finland_analysis5
DIR_PKL = DIR + '/pkl/atl08'
DIR_MAT = DIR + '/mat/atl08'

DIR_03 = '/bigtex_data/data/release/002/ATL03_r002/Finland'
DIR_08 = '/bigtex_data/data/release/002/ATL08_r002/Finland'

DIR_TEST = '/LIDAR/server/poseidon_files/USERS/jsipps/scripts_lidar/icesat_analysis'


overwrite = True
# if True, does something
# otherwise, code just pretends to do something

args = sys.argv[1:]
if len(args) > 0:
	DIR = args[0]
	DIR_PKL = DIR + '/pkl/atl08'
	DIR_MAT = DIR + '/mat/atl08'

	DIR_03 = args[1]
	DIR_08 = args[2]

# print(DIR)
# print(DIR_08)
# sys.exit()

# add
#	doy
#	beam number
#	weak/strong flag

##### Function to identify the segment id for each atl03 photon return
def getAtl03SegID(atl03_ph_index_beg, atl03_segment_id, atl03_heights_len):

	# Filter all data where atl03_ph_index starts at 0 (0 indicates errors)
	indsNotZero = atl03_ph_index_beg != 0
	atl03_ph_index_beg = atl03_ph_index_beg[indsNotZero]
	atl03_segment_id = atl03_segment_id[indsNotZero]

	# Subtract 1 from ph_index_beg to start at python 0th pos
	atl03_ph_index_beg = atl03_ph_index_beg - 1

	# Sometimes the ph_index_beg is not at the 0th position, it is is not,
	# add it in and then add the associated segment id
	# Warning, this is assuming that the segment id for the points are from 
	# the segment id directly before it, this assumption might fail but I have
	# not come across a case yet where it does.  If you want to play it safe
	# you could comment this section out and then if the first position is not
	# 0 then all photons before the first position will not be assigned a
	# segment id.
	# if atl03_ph_index_beg[0] != 0:
	#     atl03_ph_index_beg = np.append(0,atl03_ph_index_beg)
	#     first_seg_id = atl03_segment_id[0] -1
	#     atl03_segment_id = np.append(first_seg_id,atl03_segment_id)


	# Append atl03_height_len to end of array for final position
	atl03_ph_index_beg = np.append(atl03_ph_index_beg,atl03_heights_len)

	# Make array equal to the length of the atl03_heights photon level data
	ph_segment_id = np.zeros(atl03_heights_len).astype(int)

	# Iterate through ph_index_beg, from the first to second to last number
	# and set the photons between ph_index_beg i to ph_index_beg i + 1 to
	# segment id i
	for i in range(0,len(atl03_ph_index_beg) - 1):
		ph_segment_id[atl03_ph_index_beg[i]:atl03_ph_index_beg[i+1]] =\
			atl03_segment_id[i]

	# Return list of segment_id at the photon level
	return ph_segment_id



def write_pickle(data, filename):
	import pickle
	fp = open(filename, 'wb')
	pickle.dump(data, fp)
	fp.close()

def read_pickle(filename):
	import pickle
	fp = open(filename, 'rb')
	data = pickle.load(fp)
	fp.close()
	return data

def get_pkl_meta(file):
	file_info = file[:-4]
	gt_index = file_info.index('gt')
	file_sub = file_info[:gt_index-1] + '.h5'
	other = file_info[gt_index:]

	try:
		other_split = other.split('_')
		gt, t_start, t_end = other_split[0], float(other_split[1]), float(other_split[2])
	except ValueError:
		gt = other_split[0]
		t_start, t_end = -1e12, 1e12
		print('warning: ValueError')

	return file_sub, gt, t_start, t_end


def get_date(*args, debug=0):
	
	"""
	Example:
	import cmglam as cm
	year, doy = cm.get_date(year, month, day)
	or
	year, month, day = cm.get_date(year, doy)

	"""
	# def help():
	# 	print("""
	# Example:
	# import cmglam as cm
	# year, doy = cm.get_date(year, month, day)
	# or
	# year, month, day = cm.get_date(year, doy))
	# 		""")


	import datetime

	if len(args) == 2:
		y = int(args[0]) #int(sys.argv[1])
		d = int(args[1]) #int(sys.argv[2])
		yp1 = datetime.datetime(y+1, 1, 1)
		y_last_day = yp1 - datetime.timedelta(days=1)
		if d <= y_last_day.timetuple().tm_yday:
			date = datetime.datetime(y, 1, 1) + datetime.timedelta(days=d-1)
			return date.month, date.day
		else:
			print("error")
			print("Michael says, \"There aren't {} days in {}. Git good.\"".format(d, y))
			# help()

	elif len(args) == 3:
		y, m, d = args
		date = datetime.datetime(y, m, d)
		doy = int(date.timetuple().tm_yday)
		# print("doy = {}".format(date.timetuple().tm_yday))
		return str(y), str(doy).zfill(3)

	else:
		print("error: incorrect number of args")
		# help()

def get_h5_meta(h5_file, meta='date', rtn_doy=False, rtn_hms=True, file_start='', debug=0):
	# ATL03_20181016000635_02650109_200_01.h5

	h5_file = h5_file.split('/')[-1]

	meta = meta.lower()

	i0 = 0
	try:
		i0 = h5_file.index(file_start)
	except ValueError:
		if debug:
			print('warning: substring %s not found in %s' % (file_start, h5_file))

	if meta == 'date':
		year = int(h5_file[i0+6:i0+10])
		month = int(h5_file[i0+10:i0+12])
		day = int(h5_file[i0+12:i0+14])

		if rtn_doy:
			year0, doy0 = get_date(year, month, day)
			return str(year0), str(doy0).zfill(3)

		return str(year), str(month).zfill(2), str(day).zfill(2)

	else:
		print('error: unknown meta=%s' % meta)
		return 0


def get_new_attrs(file_08, gt):

	spot_number = 'unknown'
	beam_type = 'unknown'
	err = 0
	with h5.File(file_08, 'r') as fp_08:
		attrs = fp_08[gt].attrs
		if 'atlas_spot_number' in attrs and 'atlas_beam_type' in attrs:
			spot_number = attrs['atlas_spot_number'].decode()
			beam_type = attrs['atlas_beam_type'].decode()
		else:
			print('warning: atlas_spot_number or atlas_beam_type not found')
			# continue
			err = 1

	year, doy = get_h5_meta(file_08, meta='date', rtn_doy=True, file_start='ATL', debug=0)

	return doy, spot_number, beam_type, err

def convert_df_to_mat(df,outfilename):
	from scipy import io
	comps =  outfilename.split('.')
	if comps[-1] != 'mat':
		outfilename = outfilename + ".mat"
	# scipy.io.savemat(outfilename, {'struct':df.to_dict("list")})
	io.savemat(outfilename, {'struct':df.to_dict("list")})


def get_index(t, t0, tol=2.0, err_msg=True):
	"""
	Gets closest index of t0 in array t, i.e. the index of time t0
	in time array t

	tol is used in case t0 is not found; default is tol = 2.0 sec
	"""
	if type(t) != np.ndarray:
		t = np.array(t)
	index = np.argmin(abs(t - t0))
	if abs(t[index] - t0) > tol:
		if err_msg:
			print('warning: index not found: looking for %4.2f vs. found %4.2f' % (t0, t[index]))
		# return -1
	return index



print('%s: adding meta data to ATL08 files..' % sys.argv[0])
if not overwrite:
	print('%s: overwrite == False, so code is not actually doing anything btw' % sys.argv[0])


files_08 = sorted(os.listdir(DIR_08))
files_08_pkl = sorted(os.listdir(DIR_PKL))

num_files = len(files_08_pkl)
# for f, file_pkl in enumerate(files_08_pkl):
init = False
dt_iter = []
for f in tqdm(range(len(files_08_pkl))):
# for f in range(len(files_08_pkl)):

	if init:
		dt_iter.append(time.time() - t1)
		# print(dt_iter[-1])
	t1 = time.time()

	# if f < 4950:
	# 	continue

	# if init:
	# 	sys.exit()

	file_pkl = files_08_pkl[f]
	init = True
	# print(file_pkl)
	# if file_pkl.startswith('ATL08_20181016133637_02740103_002_01'):
	# 	continue
	# else:
	# 	init = True
	# print('%d/%d' % (f+1, num_files))

	if not ('gt' in file_pkl):
		print('gt not found')
		print(file_pkl)
		# input('')
		continue

	# find original 08 file_pkl
	file_08_sub, gt_08, t_start, t_end = get_pkl_meta(file_pkl) #, ftype='pkl')
	file_08 = DIR_08 + '/' + file_08_sub
	gt = gt_08
	if not os.path.exists(file_08):
		print('atl08 file not found')
		print(file_08)
		# input('')
		continue

	# find original 03 file
	file_as_03 = list(file_pkl)
	file_as_03[4] = '3'
	file_as_03 = ''.join(file_as_03)
	file_03_sub, gt_03, t_start_03, t_end_03 = get_pkl_meta(file_as_03) #, ftype='pkl')
	file_03 = DIR_03 + '/' + file_03_sub
	if not os.path.exists(file_03):
		print('atl03 file not found')
		print(file_03)
		# continue
		# if not found, try matching up through date, but not release/extra stuff
		b = 0
		for file in os.listdir(DIR_03):
			if file.startswith(file_03_sub[:29]):
				b = 1
				break
		if not b:
			continue

		file_03 = DIR_03 + '/' + file
		# file_03_test = cm.search(DIR_03, [file_03_sub[:29]], ext='h5')
		# if len(file_03_test) == 0:
		# 	continue
		# file_03 = file_03_test[0]
		print('warning: using alt file_03')
		print(file_03)
		print(file_08)

	# if file_03 != '/bigtex_data/data/release/002/ATL03_r002/Finland/ATL03_20181019013011_03120105_002_01.h5':
	# 	continue

	if gt_03 != gt_08:
		print('gt_03 != gt_08')
		print(file_03)
		print(file_08)
		continue

	atl03 = icesatReader.get_atl03_struct(file_03, gt, file_08)
	atl03_g = icesatReader.read_atl03_geolocation(file_03, gt)
	atl08 = icesatReader.get_atl08_struct(file_08, gt, atl03)
	# atl08 = icesatReader.get_atl08_struct(file_08, gt)


	# fp = h5.File(file_08, 'r')
	# delta_time_08_h5 = np.array(fp[gt + '/land_segments/delta_time'])
	# t_08_h5 = delta_time_08_h5 - delta_time_08_h5[0]
	# fp.close()

	seg_id_beg_08 = np.array(atl08.df.segment_id_beg).astype(int)
	seg_id_end_08 = np.array(atl08.df.segment_id_end).astype(int)

	t = np.array(atl03.df.time).astype(float)
	# t_08 = np.array(atl08.df.time).astype(float)

	# t_08 -= t_08[0]
	# t_08 = np.array(atl08.df.delta_time).astype(float)
	# t_08 = t_08 - t_08[0]

	# if len(t_08) != len(seg_id_beg_08):
	# 	print('error: t_08')
	# 	print(file_03)
	# 	print(file_08)
	# 	continue

	df_08 = read_pickle(DIR_PKL + '/' + file_pkl)
	# df_time = np.array(df_08.time) + t_start

	c = np.array(atl03.df.classification)
	nan_index = np.where(np.isnan(c))
	c[nan_index] = -1
	c = c.astype(int)

	index_beg = np.array(atl03_g.ph_index_beg)
	seg_id_beg = np.array(atl03_g.segment_id).astype(int)
	ph_count = np.array(atl03_g.segment_ph_cnt).astype(int)

	valid = index_beg != 0
	index_beg = index_beg[valid]
	seg_id_beg = seg_id_beg[valid]
	ph_count = ph_count[valid]

	index_beg -= 1

	df_time_08 = np.array(df_08.delta_time)
	delta_time_08 = np.array(atl08.df.delta_time)
	t_mu = []
	w = len(df_08.delta_time)
	n_08 = len(delta_time_08)

	for j in range(n_08-w+1):
		val = np.mean(abs(delta_time_08[j:j+w] - df_time_08))
		t_mu.append(val)

	j0 = np.argmin(t_mu)
	j_index = np.arange(j0,j0+w,1)

	dt_err = delta_time_08[j_index] - df_time_08

	seg_id_beg_08 = seg_id_beg_08[j_index]
	seg_id_end_08 = seg_id_end_08[j_index]

	if len(np.unique(dt_err)) != 1:
		print('error: dt_err')
		print(file_03)
		print(file_08)
		input('enter to continue')
		continue

	rdm_ground = np.full(seg_id_beg_08.shape, np.nan)
	rdm_veg = np.full(seg_id_beg_08.shape, np.nan)
	rdm_canopy = np.full(seg_id_beg_08.shape, np.nan)
	n_shots_unique = np.full(seg_id_beg_08.shape, np.nan)

	n_id = len(seg_id_beg)
	for s in range(len(seg_id_beg_08)):
		k0 = get_index(seg_id_beg, seg_id_beg_08[s], tol=0.1, err_msg=False)
		k1 = get_index(seg_id_beg, seg_id_end_08[s], tol=0.1, err_msg=False)

		# f numbers for edge cases
		# 74, 2253, 2261
		
		warn = False
		b_edge = False
		if seg_id_beg[k0] < seg_id_beg_08[s]:
			# left side incomplete
			# cm.pause('beg')
			k = k0
			while seg_id_beg[k] < seg_id_beg_08[s]:
				k += 1
				if k >= n_id:
					b_edge = True
					break

		elif seg_id_beg[k0] > seg_id_beg_08[s]:
			# print('warning: 03 seg id beg %d > 08 seg id beg %d' % (seg_id_beg[k0], seg_id_beg_08[s]))
			warn = True

		# else:
			# equal, totally fine

		# if seg_id_beg[k1] != seg_id_end_08[s]:
		if seg_id_beg[k1] > seg_id_end_08[s]:
			# right side incomplete
			# cm.pause('end')
			k = k1
			while seg_id_beg[k] > seg_id_end_08[s]:
				k -= 1
				if k < 0:
					b_edge = True
					break

		elif seg_id_beg[k1] < seg_id_end_08[s]:
			# print('warning: 03 seg id beg %d < 08 seg id beg %d' % (seg_id_beg[k0], seg_id_beg_08[s]))
			warn = True

		# else:
			# equal, totally fine

		if b_edge:
			# 08 segment is entirely outside of 03 segment data
			continue

		if warn:
			print('')
			print(f)
			print('03: [%d, %d]' % (seg_id_beg[k0], seg_id_beg[k1]))
			print('08: [%d, %d]' % (seg_id_beg_08[s], seg_id_end_08[s]))
			# cm.pause()

		i0, i1 = index_beg[k0], index_beg[k1] + ph_count[k1] - 1

		t_seg = t[i0:i1+1] # inclusive index
		c_seg = c[i0:i1+1]

		n_shots_total_uq = len(np.unique(t_seg))
		n_shots_ground = (c_seg == 1).sum()
		n_shots_veg = (c_seg == 2).sum()
		n_shots_canopy = (c_seg == 3).sum()

		# n_shots_total.append(n_shots_total_uq)
		# rdm_ground.append( float(n_shots_ground / n_shots_total_uq) )
		# rdm_veg.append( float(n_shots_veg / n_shots_total_uq) )
		# rdm_canopy.append( float(n_shots_canopy / n_shots_total_uq) )

		n_shots_unique[s] = n_shots_total_uq
		rdm_ground[s] = float(n_shots_ground / n_shots_total_uq)
		rdm_veg[s] = float(n_shots_veg / n_shots_total_uq)
		rdm_canopy[s] = float(n_shots_canopy / n_shots_total_uq)


	######################################
	# add doy, beam_type, beam_number to 08 pkl/mat files

	doy, spot_number, beam_type, err = get_new_attrs(file_08, gt)
	if err:
		# input('')
		continue

	# file_copy = file_pkl[:-4] + '-copy.pkl'
	# file_mat = file_pkl[:-4] + '-copy.mat'
	# file_copy = file_pkl[:-4] + '-copy.pkl'
	file_mat = file_pkl[:-4] + '.mat'
	if not os.path.exists(DIR_MAT + '/' + file_mat):
		print('mat file not found')
		print(file_mat)
		# input('')
		continue

	# df_08 = read_pickle(DIR_PKL + '/' + file_pkl)

	df_shape = df_08['time'].shape
	df_08['doy'] = np.full(df_shape, int(doy))
	df_08['beam_type'] = np.full(df_shape, beam_type)
	df_08['beam_number'] = np.full(df_shape, int(spot_number))

	df_08['n_shots_unique'] = n_shots_unique
	df_08['rdm_ground'] = rdm_ground
	df_08['rdm_veg'] = rdm_veg
	df_08['rdm_canopy'] = rdm_canopy

	if overwrite:
		# overwrites files with new info
		write_pickle(df_08, DIR_PKL + '/' + file_pkl)
		convert_df_to_mat(df_08, DIR_MAT + '/' + file_mat)
		# write_pickle(df_08, DIR_TEST + '/' + file_pkl)
		# convert_df_to_mat(df_08, DIR_TEST + '/' + file_mat)
		# break

	# break

	# input('enter to continue')

# df_test = read_pickle(DIR_TEST + '/' + file_pkl)
# mat_test = io.loadmat(DIR_TEST + '/' + file_mat)
