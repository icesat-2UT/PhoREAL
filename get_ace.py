#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 11:51:52 2020

@author: eguenther
"""
import numpy as np
import vox
import acer
from scipy import interpolate
from operator import itemgetter
import time
from icesatUtils import getRaster

def get_ace(x,y,z,c,i):
    np.warnings.filterwarnings('ignore')
#    print('    Start')    
    try:
        start_time = time.time()
        xminground = np.min(x[c == 2])
        xmaxground = np.max(x[c == 2])
        yminground = np.min(y[c == 2])
        ymaxground = np.max(y[c == 2])
        
        checkpoint = time.time()
#        print('    Split in-range data')    
        infilt = (((x >= xminground) & (x <= xmaxground)) & ((y >= yminground) & (y <= ymaxground)))
        xin = x[infilt]
        yin = y[infilt]
        zin = z[infilt]
        cin = c[infilt]
        iin = i[infilt]
    #    ein = e[infilt]
    #    nin = n[infilt]
        
        outfilt = np.logical_not(infilt)
        xout = x[outfilt]
        yout = y[outfilt]
        zout = z[outfilt]
        cout = c[outfilt]
        iout = i[outfilt]
        
        
        groundclass = 2
        res = 1
        
#        print("    %s seconds" % (time.time() - checkpoint))
        
        checkpoint = time.time()
        # Grid the ground points, 1m grids
#        print('    Create DEM')
        gx = x[c == groundclass]
        gy = y[c == groundclass]
        gz = z[c == groundclass]
        
        grid = getRaster(gx, gy, gz, res, 'median', fillValue = -999, time = [])
        
        # Set -999 to nan
        dem = grid.grid
        dem[dem == -999] = np.nan
    
    #    print("    %s seconds" % (time.time() - checkpoint))
        
    #    checkpoint = time.time
#        print('    Prepare for Interpolation')    
        # Derive x and y indexes from mesh-grid
        x_indx, y_indx = np.meshgrid(np.arange(0, dem.shape[1]),
                                 np.arange(0, dem.shape[0]))
        
        # Mask all invalid values
        zs_masked = np.ma.masked_invalid(dem)
        
        # Retrieve the valid, non-Nan, defined values
        valid_xs = x_indx[~zs_masked.mask]
        valid_ys = y_indx[~zs_masked.mask]
        valid_zs = zs_masked[~zs_masked.mask]
        
#        print("    %s seconds" % (time.time() - checkpoint))
        
        checkpoint = time.time()
        # Generate interpolated array of z-values
#        print('    Interpolate Missing Values')
        dem = interpolate.griddata((valid_xs, valid_ys), valid_zs.ravel(),
                                         (x_indx, y_indx), method='linear')
        del valid_xs
        del valid_ys
        del valid_zs
        del x_indx
        del y_indx
        del zs_masked
        
        
        xg = grid.x
        yg = grid.y
        zg = dem
        
        xg = xg.flatten()
        yg = yg.flatten()
        zg = zg.flatten()
        
        xg = np.reshape(xg,[len(xg),1])
        yg = np.reshape(yg,[len(xg),1])
        zg = np.reshape(zg,[len(xg),1])
        
        xc = np.reshape(xin,[len(xin),1])
        yc = np.reshape(yin,[len(xin),1])
        zc = np.reshape(zin,[len(xin),1])
        cc = np.reshape(cin,[len(xin),1])
        ic = np.reshape(iin,[len(xin),1])
        
#        print("    %s seconds" % (time.time() - checkpoint))
        
        checkpoint = time.time()    
#        print('    Voxelize, Part 1')
        _, dict1, _, _ = vox.voxelize_dict(xc,yc,zc,1.0)
        gvind, gdict1, _, gdigits  = vox.voxelize_dict(xg,yg,zg,1.0)
#        print("    %s seconds" % (time.time() - checkpoint))
        checkpoint = time.time()    
#        print('    Voxelize, Part 2')
        refl = [vox.ijk_to_voxelind(list(gvind[i]),gdigits) for i in range(0,len(xg))]
        
        indexc = itemgetter(*refl)(dict1)
        indexg = itemgetter(*refl)(gdict1)
#        print("    %s seconds" % (time.time() - checkpoint))
        
        checkpoint = time.time()    
#        print('    Reassign unclassed values')   
        cc[cc == 13] = 0
        cc[cc == 3] = 4
        cc[cc == 5] = 4
        cc = acer.aceReassign(indexc, indexg, zc, cc, zg)
        xc = xc.flatten()
        yc = yc.flatten()
        zc = zc.flatten()
        cc = cc.flatten()
        ic = ic.flatten()
    #    ein = ein.flatten()
    #    ein = ein.flatten()
        xc = np.concatenate((xc, xout))
        yc = np.concatenate((yc, yout))
        zc = np.concatenate((zc, zout))
        cc = np.concatenate((cc, cout))
        ic = np.concatenate((ic, iout)) 
    #    ein = np.concatenate((ein, eout)) 
    #    nin = np.concatenate((nin, nout)) 
        return xc, yc, zc,cc,ic
    except:
        print('ACE Failed, returning original values')
        return x,y,z,c,i