import numpy as np
import datetime



def _get_final_sample_index(orig_array=None, sampled_array=None):

    final_sampled_index = np.argwhere(orig_array == sampled_array[-1])[0][0]
    
    additional_sample = int(((orig_array.shape[0] - final_sampled_index)/2) + final_sampled_index)

    return additional_sample


def time2UTC(gps_seconds_array=None):
    
    """
        Converts array of gps seconds to UTS.
        
        Arguments
        ---------
        gps_seconds_array
        
        Returns
        -------
        utc_time_array

    """

    # Number of Leap seconds
    # See: 'https://www.ietf.org/timezones/data/leap-seconds.list'
    # 15 - 1 Jan 2009
    # 16 - 1 Jul 2012
    # 17 - 1 Jul 2015
    # 18 - 1 Jan 2017
    leap_seconds = 18

    gps_start = datetime.datetime(year=1980, month=1, day=6)
    time_ph = [datetime.timedelta(seconds=time) for time in gps_seconds_array]
    # last_photon = datetime.timedelta(seconds=gps_seconds[-1])
    error = datetime.timedelta(seconds=leap_seconds)
    utc_time_array = [(gps_start + time - error) for time in time_ph]

    return utc_time_array

def wgs84_to_web_mercator(in_lons=None, in_lats=None):
    
    # https://stackoverflow.com/questions/57178783/how-to-plot-latitude-and-longitude-in-bokeh

    k = 6378137
    lons = in_lons * (k * np.pi/180.0)
    
    lats = np.log(np.tan((90 + in_lats) * np.pi/360.0)) * k

    return lons, lats