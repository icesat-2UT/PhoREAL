#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:46:19 2019

@author: eguenther
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 08:47:00 2019

@author: eguenther
"""

import os
import numpy as np
import time as runTime
import sys
import warnings
#from icesatIO import *
#from icesatUtils import *

from phoreal.getMeasurementError import *   
from scipy import spatial
import matplotlib.pyplot as plt
import copy
import matplotlib.pyplot as plt


#Include/exclude filter
def includeFilter(sortedMeasured, superTruth):
#    xtest2 = superTruth.alongTrack
    difarr = superTruth.alongTrack[1:len(superTruth.alongTrack)] - superTruth.alongTrack[0:(len(superTruth.alongTrack) - 1)]
    difarr  = np.append(difarr,0)
    indices = np.arange(len(difarr))
    selected = [(i,j) for (i,j) in zip(difarr,indices) if i >= 0.5]
    selectlist = list(zip(*selected))
    if len(selectlist) > 0:
        selectlist = selectlist[1]
        startindex = np.array([0])
        selectlist = np.array(selectlist)
        selectlist = selectlist + 1
        startindex2 = np.append(startindex,selectlist)
        endindex = np.array(selectlist) 
        endindex = endindex - 1
        endindex = np.append(endindex,(len(superTruth.alongTrack) - 1))
        setindex = np.column_stack([startindex2,endindex])
        delindex = np.array([-1])
        for i in range(0,setindex.shape[0]):
            indexdif= setindex[i][1] - setindex[i][0]
            if (indexdif > 1000):
                if np.sum(delindex) == -1:
                    delindex = i
                else:
                    delindex = np.append(delindex,i)
                    
        
        includeindex = setindex[delindex]
        
        
        
        #Create an 0 length array of the same size as the ATL03 data
        #Iterate through the start and end index for the truth data
        #At each iteration, whereever ATL03 is between the start and stop...
        #...make array of 1's where it is between the values
        #...add that array to a global array filter
        #...use array to create a filter
        yfilter = np.zeros([len(sortedMeasured.alongTrack),1])
        if includeindex.ndim == 1:
                start =  superTruth.alongTrack[includeindex[0]]
                stop = superTruth.alongTrack[includeindex[1]]
                yfilter[(sortedMeasured.alongTrack > start)  \
                        & (sortedMeasured.alongTrack < stop)] = 1        
        else:
            for i in range(0,len(includeindex)):
                start =  superTruth.alongTrack[includeindex[i][0]]
                stop = superTruth.alongTrack[includeindex[i][1]]
                yfilter[(sortedMeasured.alongTrack > start)  \
                        & (sortedMeasured.alongTrack < stop)] = 1
    else:
        yfilter = np.zeros([len(sortedMeasured.alongTrack),1])
        start =  np.min(superTruth.alongTrack)
        stop = np.max(superTruth.alongTrack)
        yfilter[(sortedMeasured.alongTrack > start)  \
                        & (sortedMeasured.alongTrack < stop)] = 1        
        
    
    return yfilter

# Perfect Classifer 
def perfect_classifier(atl03, truth_swath,ground = [2],canopy = [4], 
                      unclassed = [1, 6, 7, 18], keepsize = True):
    # Find max/min along track
    print('Run Perfect Classifier')
    maxy = np.max(truth_swath.alongtrack)
    miny = np.min(truth_swath.alongtrack)
    
    measy = np.array(atl03.df.alongtrack)
    measz = np.array(atl03.df.h_ph)
    measc = np.array(atl03.df.classification)
    # measz = measz - 1
    
    # Reclass everything to generic classes
    measc[measc == -1] = 0
    measc[measc == 3] = 2
    
    # Calculate filter offsets for later
    if keepsize == True:
        if np.min(measy) < miny:
            minmeasindex = next(x for x, val in enumerate(measy) if val > miny)
        else:
            minmeasindex = 0
        if np.max(measy) > maxy:
            maxmeasindex = next(x for x, val in enumerate(measy) if val > maxy)
        else:
            maxmeasindex = 0

    # Apply filter based on min/max of truth data
    measfilt = np.where((measy <= maxy) & (measy >= miny))

    measyfilt = measy[measfilt]
    measzfilt = measz[measfilt]
    meascfilt = measc[measfilt]
    classfilt = [np.isin(np.array(truth_swath.classification),ground) | np.isin(np.array(truth_swath.classification),canopy)]
    ty = np.array(truth_swath.alongtrack)[classfilt[0]]
    tz = np.array(truth_swath.z)[classfilt[0]]
    tc = np.array(truth_swath.classification)[classfilt[0]]
    #Create KDtree and then query the tree
    pts = np.vstack((measyfilt,measzfilt)).T
    print('    Building KDtree')
    tree = spatial.KDTree(list(zip(ty, tz)))
    print('    Querying KDtree')
    dist, index = tree.query(pts)
    
    # Pre-populate all photon classed array with zeroes
    measpc = np.zeros(len(index))
    measpc = measpc - 1
    
    # Populate all photon classed array from ATL08 classifications
    measpc = tc[index]

    measpc[dist > 1.5] = 0

    # Reclass everything to generic classes
    measpc[np.isin(measpc,unclassed)] = 0
    measpc[np.isin(measpc,ground)] = 1
    measpc[np.isin(measpc,canopy)] = 2
    measpc[measpc > 2] = 0    

    keepsize = True

    # If keepsize is True, append 0's to measpc
    if keepsize == True:
        if ((minmeasindex == 0) & (maxmeasindex == 0)):
            measpc = measpc
        elif ((minmeasindex > 0) & (maxmeasindex == 0)):
            frontzeros = np.zeros(minmeasindex)
            measpc = np.append(frontzeros,measpc)
        elif ((minmeasindex == 0) & (maxmeasindex > 0)):
            backzeros = np.zeros(len(measc) - maxmeasindex)
            # measpc = np.append(frontzeros,measpc)
            measpc = np.append(measpc,backzeros)
        elif minmeasindex > maxmeasindex:
            frontzeros = np.zeros(maxmeasindex)
            backzeros = np.zeros(len(measc) - minmeasindex)
            measpc = np.append(frontzeros,measpc)
            measpc = np.append(measpc,backzeros)
        elif minmeasindex <= maxmeasindex:
            frontzeros = np.zeros(minmeasindex)
            backzeros = np.zeros(len(measc) - maxmeasindex)
            measpc = np.append(frontzeros,measpc)
            measpc = np.append(measpc,backzeros)
        measpc = measpc.astype('int')
        measoc = measc
    else:
        measpc = measpc.astype('int')
        measoc = meascfilt
    
    return measpc, measoc


if __name__ == "__main__":
    pass

