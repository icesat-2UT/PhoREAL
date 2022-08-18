# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 09:22:31 2022

@author: ERICG
"""

import argparse
import os
import laspy
from laspy.file import File                                                                                   
import numpy as np                                                                                       
from scipy import interpolate                                                                            
                               
def interpolate_dem(dem, nan_value=-999):                                                                
    dem[dem == nan_value] = np.nan                                                                       
    # If the dem is only 1 pixel wide in one direction it will cause a problem                           
    # with the interpolation, so check for this and pad it                                               
    invalid_x_size = False                                                                               
    invalid_y_size = False                                                                               
    if dem.shape[0] == 1:  # pragma: no cover
        dem = np.concatenate((dem, dem, dem), axis=0)                                                    
        invalid_x_size = True                                                                            
    if dem.shape[1] == 1:  # pragma: no cover                                                            
        dem = np.concatenate((dem, dem, dem), axis=1)                                                    
        invalid_y_size = True                                                                            
    if invalid_x_size and invalid_y_size:  # pragma: no cover                                            
        raise RuntimeError('Found single pixel dem, cannot process')  
    # Derive x and y indexes from mesh-grid                                                              
    x_indx, y_indx = np.meshgrid(np.arange(0, dem.shape[1]),                                             
                             np.arange(0, dem.shape[0]))                                             

    # Mask all invalid values                                                                            
    zs_masked = np.ma.masked_invalid(dem)                                                                

    # Retrieve the valid, non-Nan, defined values                                                        
    valid_xs = x_indx[~zs_masked.mask]                                                                   
    valid_ys = y_indx[~zs_masked.mask]                                                                   
    valid_zs = zs_masked[~zs_masked.mask]                                                                

    # Interpolate missing values on the DEM                                                              
    dem = interpolate.griddata((valid_xs, valid_ys), valid_zs.ravel(),                                   
                               (x_indx, y_indx), method='linear')                                        

                    
    # Remove padding if it was added                                                                     
    if invalid_x_size:  # pragma: no cover                                                               
        dem = dem[1, :]                                                                                  
        dem = np.expand_dims(dem, axis=0)                                                                
    if invalid_y_size:  # pragme: no cover                                                               
        dem = dem[:, 1]                                                                                  
        dem = np.expand_dims(dem, axis=1)                                                                

    return dem                                                                                           

# Run ACE
def ace(x, y, z, c, interpolate = False):
    # Define static variables
    ground = 2 # ASPRS ground classification code
    unclassed = 1 # ASPRS unclassed classification code
    x_res = 1 # Res of 'grid' in 'x'                                                                                      
    y_res = 1 # Res of 'grid' in 'y'
    ground_buffer = 0.5 # Default buffer of ground
    
    # Grid las files                          
    x = np.floor((x - np.min(x))/x_res).astype(int)                                                     
    y = np.floor((y - np.min(y))/y_res).astype(int)      
    # dem = np.zeros([np.max(x[c == ground]) + 1, np.max(y[c == ground]) + 1])                                                       
    dem = np.zeros([np.max(x) + 1, np.max(y) + 1])                                                       
    dem[x[c == ground],y[c == ground]] = z[c == ground]                                                                                 

    # Interpolate DEM
    if interpolate:
        dem = interpolate_dem(dem, nan_value=0)                                                              
        dem[np.isnan(dem)] = 0                                                                               

    #Normalize las.z to DEM                                                                    
    z = z - dem[x,y]                                                                                   

    #Where C is 0 and Height is between -1 and 1.                                                        
    c[np.where((c == unclassed) & (z > -ground_buffer) &\
               (z < ground_buffer))] = ground
    
    return c

# Read las/laz, run ACE, write out new las/laz
def convert_las_ace(in_file, out_file, interpolate = False):                                                                                     
    # Read the las file   
    las = File(in_file, mode = "r")
    
    # Run ACE
    c = ace(las.x, las.y, las.z, las.classification, interpolate = False)
    
    # Write and close out in/out las files
    las_out = laspy.file.File(out_file, mode="w", header=las.header)
    las_out.set_points(las.get_points())
    las.close()
    las_out.classification = c
    las_out.close()
    
def main(in_file, out_file, interpolate):
    # Check if in_file is file
    if os.path.isfile(in_file):
        # Check if out_file is a .las or .laz
        if out_file[-4:] == '.las' or out_file[-3:] == '.laz':
            convert_las_ace(in_file, out_file, interpolate)
        # Check if out_file is actually a directory
        elif os.path.isdir(out_file):
            out_file = os.path.join(out_file,\
                                os.path.basename(in_file)[:-4] + '_ace.las')
            convert_las_ace(in_file, out_file, interpolate)
    # Check if in_file and out_file are actually directories
    elif os.path.isdir(in_file) and os.path.isdir(out_file):
        for f in os.listdir(in_file):
            _, ext = os.path.splitext(f)
            if ext == '.las':
                out_file = os.path.join(out_file,\
                                os.path.basename(f)[:-4] + '_ace.las')
                convert_las_ace(f, out_file, interpolate)
            elif ext == '.laz':
                out_file = os.path.join(out_file,\
                                os.path.basename(f)[:-4] + '_ace.laz')
                convert_las_ace(f, out_file, interpolate)
                                                                      
if __name__ == '__main__':
    """ Command line entry point """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("in_file",
                        help="Input .las/.laz or directory of .las/.laz files")

    parser.add_argument("out_file",
                        help="Output .las/.laz or directory of .las/.laz files")

    parser.add_argument('-interpolate', '--interpolate_dem', action='store_true',
                        help='Interpolate missing pixels in DEM')


    args = parser.parse_args()

    args.in_file = os.path.normpath(args.in_file)
    args.out_file = os.path.normpath(args.out_file)

    main(args.in_file,
        args.out_file,
        args.interpolate)                    