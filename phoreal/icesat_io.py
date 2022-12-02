import h5py
import re

import beam_tools


def get_ATL03_metadata(granule_fname=None, phorealDS_settings=None):
    """
        Reads ICESat-2 ATL03 Global Geolocated Photons data files.

        Arguments
        ---------
        granule_fname: full path to ATL03 file
        phorealDS_settings: PhoREAL Dataset settings dictionary

        Returns
        -------
        IS2_atl03_mds: dictionary with ATL03 variables
        IS2_atl03_attrs: dictionary with ATL03 attributes
        IS2_atl03_beams: list with valid ICESat-2 beams within ATL03 file

    """
    fileID = h5py.File(granule_fname, 'r')
    #-- Open the HDF5 file for reading

    #-- allocate python dictionaries for ICESat-2 ATL03 variables and attributes
    IS2_atl03_mds = {}

    #-- read each input beam within the file
    IS2_atl03_beams = []
    # for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
    for key in fileID.keys():
        if re.match(r'gt\d[lr]',key):
            IS2_atl03_beams.append(key)

    all_beam_paths = [beam_tools.build_beam_path(fileID=fileID, gtx=gtx,
                      phorealDS_settings=phorealDS_settings) for gtx in IS2_atl03_beams]
    
    IS2_atl03_attrs = {}
    
    IS2_atl03_attrs['ancillary_data'] = {}

    ancillary_keys = ['atlas_sdp_gps_epoch','data_end_utc','data_start_utc',
        'end_cycle','end_geoseg','end_gpssow','end_gpsweek','end_orbit',
        'end_region','end_rgt','granule_end_utc','granule_start_utc','release',
        'start_cycle','start_geoseg','start_gpssow','start_gpsweek',
        'start_orbit','start_region','start_rgt','version']
    for key in ancillary_keys:
        #-- get each HDF5 variable
        IS2_atl03_attrs['ancillary_data'][key] = fileID['ancillary_data'][key][:]

    #-- Closing the HDF5 file
    fileID.close()

    return all_beam_paths, IS2_atl03_attrs['ancillary_data']