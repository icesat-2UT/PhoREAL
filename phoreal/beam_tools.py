import numpy as np
import utils





def build_beam_path(fileID=None, gtx=None, phorealDS_settings=None):
        
        """
            For an open H5 file reader object fileID, read the lat_ph,
            lon_ph and delta_time arrays for a given ground track string
            gtx ('gt1l', 'gt3r', etc). Extract points from the array every
            n array elements defined in the PhorealDataSet settings attribute
            under the key 'granule_sample_size'.
            
            Sometimes one ground track returns arrays that are one element
            shorter than others, so there is a helper function
            _get_final_sample_index() that ensures all ground tracks have the
            same number of elements.

            Arguments
            ---------

                fileID: Open H5 reader object.
                gtx: ground track key.

            Returns
            -------

                beam_path_dict

        """

        lats = fileID[gtx]['heights']['lat_ph'][:]
        lons = fileID[gtx]['heights']['lon_ph'][:]
        times = fileID[gtx]['heights']['delta_time'][:]
        
        every_n = phorealDS_settings['granule_sample_size']
        
        lats_every_n = int(lats.shape[0]/every_n)
        lons_every_n = int(lats.shape[0]/every_n)
        times_every_n = int(lats.shape[0]/every_n)

        lat_strt, lat_end = lats[0], lats[-1]
        mid_lats = lats[::lats_every_n]
        
        lon_strt, lon_end = lons[0], lons[-1]
        mid_lons = lons[::lons_every_n]
        
        delta_t_strt, delta_t_end = times[0], times[1]
        mid_delta_t = times[::times_every_n]
        
        if mid_lats.shape[0]-1 == every_n-1:
            new_ind = utils._get_final_sample_index(orig_array=lats, sampled_array=mid_lats)
            mid_lats = np.hstack((mid_lats, lats[new_ind]))

        if mid_lons.shape[0]-1 == every_n-1:
            mid_lons = np.hstack((mid_lons, lons[new_ind]))
        
        if mid_delta_t.shape[0]-1 == every_n-1:
            mid_delta_t = np.hstack((mid_delta_t, times[new_ind]))

        lons = np.hstack((mid_lons, np.asarray([lon_end])))
        lats = np.hstack((mid_lats, np.asarray([lat_end])))
        times = np.hstack((mid_delta_t, np.asarray([delta_t_end])))
        
        beam_path_dict = {'lons': lons,
                          'lats': lats,
                          'times': times,
                          'beam_name': gtx}

        return beam_path_dict



def get_granule_path(granule_info_dict=None, granule_name=None):
        
    """
        For beam paths stored in PhoREAL Dataset granule_info[granule_name]['beam_paths']
        extract all lons, lats and times to a 2D array and
        take the mean of each to have a mean path of a granule.

        Arguments
        ---------

            granule_name: name of the granule to use as key in
                            self.granule_info dictionary.

        Returns
        -------

            mean granule_path_dict

    """
    
    lons = []
    lats = []
    times = []
    
    for beam in granule_info_dict[granule_name]['beam_paths']:

        lons.append(beam['lons'])
        lats.append(beam['lats'])
        times.append(beam['times'])
        
    lons = np.asarray(lons)
    lats = np.asarray(lats)
    times = np.asarray(times)
    
    lons = np.mean(lons, axis=0)
    lats = np.mean(lats, axis=0)
    times = np.mean(times, axis=0)
        
    return {'lons': lons, 'lats': lats, 'times': times}
