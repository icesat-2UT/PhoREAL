import glob
import re
import numpy as np
import h5py

import utils
import icesat_io
import beam_tools



class PhorealDataSet():
    
    
    def __init__(self):
        
        self.granule_list = []
        self.granule_info = {}
        self.granule_data = {}
        self.settings = {'granule_sample_size':10}
    
    def build_granules_metadata(self):
        
        """
            Populates self.granule_info['granule_name'] for
            every granule in self.granule_list
            
            Returns
            -------

                populated self.granule_info dictionary

        """

        for granule_path in self.granule_list:

            granule_name = granule_path.split('/')[-1].split('.')[0]

            granule_beam_paths, granule_meta = icesat_io.get_ATL03_metadata(granule_path, self.settings)

            self.granule_info[granule_name] = {
                                               'granule_meta': granule_meta,
                                               'beam_paths': granule_beam_paths
                                              }
            
            self.granule_data[granule_name] = {}
            
            granule_path = beam_tools.get_granule_path(granule_info_dict=self.granule_info,
                                                       granule_name=granule_name)

            self.granule_info[granule_name]['granule_paths'] = granule_path

        return
    
    
    def get_all_granule_paths(self, convert_units=True):
        
        """
            For granule in self.granule_list create a list of 
            dictionaries containing the sampled path lons, lats,
            times. Optionally create wgs_lons, wgs_lats by setting
            convert_units=True (default).
    
            Arguments
            ---------

                convert_units: add wgs_lons and wgs_lats keys to the
                               granule path dictionary.

            Returns
            -------

                list of all granule paths in self.granule_list.

        """

        all_granule_paths = []

        for granule_fname in self.granule_list:

            granule_name = granule_fname.split('/')[-1].split('.')[0]

            path_dict = self.granule_info[granule_name]['granule_paths']
            
            gps_epoch = self.granule_info[granule_name]['granule_meta']['atlas_sdp_gps_epoch'][0]
            
            if not isinstance(path_dict['times'][0], datetime.datetime):
            
                utc_times = utils.time2UTC(gps_seconds_array=path_dict['times'] + gps_epoch)
            else:
                utc_times = path_dict['times']
                
            if convert_units:

                wgs_lons, wgs_lats = utils.wgs84_to_web_mercator(in_lons=path_dict['lons'], in_lats=path_dict['lats'])

                path_dict = {'wgs_lons': wgs_lons, 'wgs_lats': wgs_lats}
                
            path_dict['times'] = utc_times
            path_dict['granule_name'] = granule_name
            
            all_granule_paths.append(path_dict)

        return all_granule_paths
    
    
    def load_directory(self, directory=None):
        
        """
            Add all H5 file filenames in directory to the 
            self.granule_list and sort them.
    
            Arguments
            ---------

                directory: working directory with input granules.

            Returns
            -------

                populates self.granule_list.

        """
        
        self.granule_list = glob.glob(directory + '*.h5')
        self.granule_list.sort()
        
        return
