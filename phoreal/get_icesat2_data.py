import icepyx as ipx
from icepyx.core.visualization import Visualize
import os
import shutil
import matplotlib
import argparse


"""
Select and download specific IceSAT-2 data granule from Earthdata using icepyx.

Command line input example below:

>>> python3 get_icesat2_data.py -p ATL08 -lat 61.92 63 -lon 25 28.7 -d 2021-04-02 2021-04-04 
    -uep my_uid me@gmail.com /home/me/data_dir

You MUST include a temporal (date range) and/or orbital (cycles, trakcs) boundary for the dataset.
"""

def parse_arguments() -> dict:
    '''Command line entry point'''
    parser = argparse.ArgumentParser()
    '''Required arguments:'''
    parser.add_argument('-p', '--product'
                        , required=True
                        , type=str
                        , help='Data product of interest (i.e., ATL03)')
    
    parser.add_argument('-lat', '--latitude'
                        , required=True
                        , nargs='+'
                        , type=float
                        , help='Min/max latitude for data selection \
                                in decimal degree format')
    
    parser.add_argument('-lon', '--longitude'
                        , required=True
                        , nargs='+'
                        , type=float
                        , help='Min/max longitude for data selection \
                                (i.e., 45.930351153400515 46.02060957038562)')
        # @TODO may refactor/separate the login credentials and download path later
        # but will keep Earthdata input together for now
    parser.add_argument('-uep', '--uid_email_opath'
                        , required=True
                        , nargs='+'
                        , type=str
                        , help='Earthdata login user ID, email for \
                                status updates, and output path for \
                                granule download \
                                (dev_team your_email@dom.com /dir/).')
    '''Optional arguments:'''
    parser.add_argument('-s', '--spatial_extent'
                        , nargs ='+'
                        , type=str
                        , default='bb'
                        , help='Regional type for data selection \
                                (bb==bounding box, pv==polygon \
                                vertices, pf==polygon file). Provide \
                                full path and name of kml, shp, \
                                or gpkg polygon file (-s pf /dir/fn).')
    
    parser.add_argument('-d', '--date_range'
                        , nargs='+'
                        , type=str
                        , default=None
                        , help='Date range for data selection \
                                (YYYY-MM-DD YYYY-MM-DD).')
    
    parser.add_argument('-c', '--cycles'
                        , nargs='+'
                        , type=str
                        , default=None
                        , help='Orbital cycle(s) for data selection \
                                (i.e., 03 or 03 04 05...).')
    
    parser.add_argument('-tr', '--tracks'
                        , nargs='+'
                        , type=str
                        , default=None
                        , help='Reference Ground Track (RGT) for data \
                                selection (i.e., 0849 or 0849 0902...).')
    
    parser.add_argument('-st', '--start_time'
                        , type=str
                        , default=None
                        , help='HH:MM:SS. Only for temporally limited data selection.')
    
    parser.add_argument('-et', '--end_time'
                        , type=str
                        , default=None
                        , help='HH:MM:SS. Only for temporally limited data selection.')
    
    parser.add_argument('-v', '--version'
                        , type=str
                        , default=None
                        , help='Version of the data product to use (i.e., 002).')
    
    args = parser.parse_args()
        
        # Spatial extent needs to be redefined to fit either bounding box,
        # polygon vertices, or file.
        # Not doing a whole lot of checking here because query.py handles 
        # incompatible input.
    if args.spatial_extent == 'bb':
        spatial_extent = [  min(args.longitude)
                          , min(args.latitude)
                          , max(args.longitude)
                          , max(args.latitude)  ]
    elif args.spatial_extent == 'pv':
        spatial_extent = [   (min(args.longitude), min(args.latitude))
                           , (min(args.longitude), max(args.latitude))
                           , (max(args.longitude), max(args.latitude))
                           , (max(args.longitude), min(args.latitude))
                           , (min(args.longitude), min(args.latitude))  ]
    else:
        spatial_extent = args.spatial_extent[1]
        
        # Prefer passing args this way here. Number of key, value pairs
        # need to be consistent so that icepyx can handle errors correctly.
        # (icepyx is written to handle cases where arguments are None)
    args_dict = {  'product'        : args.product
                 , 'spatial_extent'  : spatial_extent
                 , 'latitude'        : args.latitude
                 , 'longitude'       : args.longitude
                 , 'date_range'      : args.date_range
                 , 'cycles'          : args.cycles
                 , 'tracks'          : args.tracks
                 , 'start_time'      : args.start_time
                 , 'end_time'        : args.end_time
                 , 'version'         : args.version
                 , 'uid_email_opath' : args.uid_email_opath }
    return args_dict

def main(args_dict):
    print('args_dict: ', args_dict)
    region = ipx.Query(   args_dict['product']
                        , args_dict['spatial_extent']
                        , args_dict['date_range']
                        , args_dict['start_time']
                        , args_dict['end_time']
                        , args_dict['version']
                        , args_dict['cycles']
                        , args_dict['tracks']   )
    
    region.earthdata_login(  args_dict['uid_email_opath'][0]
                           , args_dict['uid_email_opath'][1]  )
    
    region.download_granules( args_dict['uid_email_opath'][2] )
    
    #cyclemap, rgtmap = region.visualize_elevation()
    #cyclemap
    #rgtmap
    
if __name__ == '__main__':
    args_dict = parse_arguments()
    main(args_dict)  