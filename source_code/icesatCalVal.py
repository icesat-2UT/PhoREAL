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
from getAtlMeasuredSwath_auto import getAtlMeasuredSwath
from getAtlTruthSwath_ACE import getAtlTruthSwath
from icesatIO import *
from icesatUtils import *

from getMeasurementError_auto import *   
from scipy import spatial
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
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



#Perfect Classifier
def perfectClassifier(sortedMeasured, superTruth,ground = [2],canopy = [4], 
                      unclassed = [1, 6, 7, 18], keepsize = True):
    # Find max/min along track
    print('Run Perfect Classifier')
    maxy = np.max(superTruth.alongTrack.ravel())
    miny = np.min(superTruth.alongTrack.ravel())
    
    measy = sortedMeasured.alongTrack.ravel()
    measz = sortedMeasured.z.ravel()
    measc = sortedMeasured.classification.ravel()
    measz = measz - 1
    
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
#        if (minmeasindex == 0) and (maxmeasindex == 0):
#            minmeasindex = next(x for x, val in enumerate(measy) if val > miny) 
#            maxmeasindex = next(x for x, val in enumerate(measy) if val > maxy)
#            if (minmeasindex == 0) and (maxmeasindex == 0):
#                print("Error, no indices in range")

    # Apply filter based on min/max of truth data
    measfilt = np.where((measy <= maxy) & (measy >= miny))

    measyfilt = measy[measfilt]
    measzfilt = measz[measfilt]
    meascfilt = measc[measfilt]
    classfilt = [np.isin(superTruth.classification,ground) | np.isin(superTruth.classification,canopy)]
    ty = superTruth.alongTrack[classfilt]
    tz = superTruth.z[classfilt]
    tc = superTruth.classification[classfilt]
    #Create KDtree and then query the tree
    pts = np.vstack((measyfilt,measzfilt)).T
    print('    Building KDtree')
    tree = spatial.KDTree(list(zip(ty.ravel(), tz.ravel())))
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



def classificationReport(sortedMeasured, superTruth):
    # Run Perfect Classifier
    measpc, measoc = perfectClassifier(sortedMeasured, superTruth,ground = [2],
                                       canopy = [4], unclassed = [1, 6, 7, 18], 
                                       keepsize = True)
    # Create includeFilter
    measFilter = includeFilter(sortedMeasured, superTruth)
    measFilter = measFilter.ravel()
    measFilter = np.where(measFilter == 1)
    # Apply includeFilter
#    atl03_filter = np.where(includeFilter == 1)
    measpcfilter = measpc[measFilter] 
    measocfilter = measoc[measFilter]
    
    # Run classification Report
    
    report = classification_report(measpcfilter, measocfilter)
    print(report)
    return report



def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
#    classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[0]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')
    ax.axis('scaled')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax


if __name__ == "__main__":

    print('icesatError')
    import numpy as np
    #import time as runTime
    from getAtlMeasuredSwath import getAtlMeasuredSwath, atlRotationStruct
    from getAtlTruthSwath_ACE import getAtlTruthSwath
    from getMeasurementError_auto import getMeasurementError, offsetsStruct
    from icesatUtils import superFilter, createScalogram, createScalogramCanopy
    from icesatError import includeFilter, perfectClassifier
    import matplotlib.pyplot as plt
    import scipy.io as sio
    #import os
    import copy
    import matplotlib.pyplot as plt
    import os

    gtNum =  'gt1l'
#    atl03FilePath = '/laserpewpew/data/release/002/ATL03_r002/Sweden/ATL03_20181015140217_02590103_002_01.h5'
#    atl08FilePath = '/laserpewpew/data/release/002/ATL08_r002/Sweden/ATL08_20181015140217_02590103_002_01.h5'
#    atl03FilePath = '/laserpewpew/data/release/002/ATL03_r002/Maine/processed_ATL03_20190220020027_08190206_002_01.h5'
##    atl08FilePath = '/laserpewpew/data/release/002/ATL08_r002/Maine/processed_ATL08_20190220020027_08190206_002_01.h5'    
#    atl03FilePath = '/laserpewpew/data/release/002/ATL03_r002/Maine/processed_ATL03_20190224015209_08800206_002_01.h5'
#    atl08FilePath = '/laserpewpew/data/release/002/ATL08_r002/Maine/processed_ATL08_20190224015209_08800206_002_01.h5'    
##    atl03FilePath = '/laserpewpew/data/release/R002/ATL03_rR002/ATL03_20191024174938_04250506_R002_01_sreq_3090.h5'
##    atl08FilePath = '/laserpewpew/data/release/R002/ATL08_rR002/ATL08_20191024174938_04250506_R002_01_sreq_3090.h5'    
#    atl03FilePath = '/laserpewpew/data/release/002/ATL03_r002/Sweden/ATL03_20181015140217_02590103_002_01.h5'
#    atl08FilePath = '/laserpewpew/data/release/002/ATL08_r002/Sweden/ATL08_20181015140217_02590103_002_01.h5'
#    atl03FilePath = '/laserpewpew/data/release/002/ATL03_r002/Maine/processed_ATL03_20190220020027_08190206_002_01.h5'
##    atl08FilePath = '/laserpewpew/data/release/002/ATL08_r002/Maine/processed_ATL08_20190220020027_08190206_002_01.h5'    
#    atl03FilePath = '/laserpewpew/data/release/002/ATL03_r002/Maine/processed_ATL03_20190224015209_08800206_002_01.h5'
#    atl08FilePath = '/laserpewpew/data/release/002/ATL08_r002/Maine/processed_ATL08_20190224015209_08800206_002_01.h5'    
##    atl03FilePath = '/laserpewpew/data/release/R002/ATL03_rR002/ATL03_20191024174938_04250506_R002_01_sreq_3090.h5'
##    atl08FilePath = '/laserpewpew/data/release/R002/ATL08_rR002/ATL08_20191024174938_04250506_R002_01_sreq_3090.h5'    


#    atl03FilePath = '/laserpewpew/data/release/002/ATL03_r002/Sonoma/ATL03_20190328035403_13700206_002_01.h5'
#    atl08FilePath = '/laserpewpew/data/release/002/ATL08_r002/Sonoma/ATL08_20190328035403_13700206_002_01.h5'    

#
    atl03FilePath = '/laserpewpew/data/release/002/ATL03_r002/Sonoma/ATL03_20190228170214_09510202_002_01.h5'
    atl08FilePath = '/laserpewpew/data/release/002/ATL08_r002/Sonoma/ATL08_20190228170214_09510202_002_01.h5'    
    truthSwathDir = '/LIDAR/server/USERS/eric/ecosystem/Sonoma/'  # Path to existing truth data
    outFilePath = '/LIDAR/server/USERS/eric/ecosystem/'
    
    comps = atl03FilePath.split('/')
    filename = comps[-1].split('.')
    title = filename[0] + "_" + gtNum
    # User options
    trimInfo = 'auto'   # OPTIONS: ('none', 'auto', or 'manual')
                        # none - does not trim data at all
                        # auto - trims ATL03 track to extent of truth bounding region
                        # manual - trims ATL03 track by latitude or time
                            # Example: 'manual,lat,38,39'
                            # Only uses data between latitude 38 and 39 deg
                            # Example: 'manual,time,3,4'
                            # Only uses data between time 3 and 4 seconds
    
    createAtl03LasFile = True    # Option to create output measured ATL03 .las file
    createAtl03KmlFile = False   # Option to create output measured ATL03 .kml file
    createAtl08KmlFile = False    # Option to create output measured ATL08 .kml file
    createAtl03CsvFile = False    # Option to create output measured ATL03 .csv file
    createAtl08CsvFile = False    # Option to create output measured ATL08 .csv file
        
    ##### End Inputs for getAtlMeasuredSwath
    
    
    ##### Start Inputs for getAtlTruthSwath
    
    buffer = 50                 # Distance in cross-track (meters) around ATL03 track to look for truth data 
    useExistingTruth = False     # Option to use existing truth data if it exists
    truthSwathDir = '/LIDAR/server/USERS/eric/ecosystem/Sonoma'  # Path to existing truth data
    createTruthFile = True      # Option to create output truth .las file
    
    ##### End Inputs for getAtlTruthSwath
    
    
    ##### CODE BELOW -- DO NOT EDIT ###############################################
    
    timeStart = runTime.time()
    
    # Call getAtlMeasuredSwath
    print('RUNNING getAtlMeasuredSwath...\n')
    atl03Data, atl08Data, headerData, rotationData = getAtlMeasuredSwath(atl03FilePath, 
                                                                     atl08FilePath, 
                                                                     outFilePath, 
                                                                     gtNum, 
                                                                     trimInfo, 
                                                                     createAtl03LasFile, 
                                                                     createAtl03KmlFile, 
                                                                     createAtl08KmlFile,
                                                                     createAtl03CsvFile, 
                                                                     createAtl08CsvFile, 
                                                                     kmlBoundsFilePath = 'kmlBounds_Sonoma_las.txt')

#    bb = np.where((atl03Data.alongTrack > 221100) & (atl03Data.alongTrack < 228480))
#    ytrack = np.ravel(atl03Data.alongTrack[bb])
#    h_ph = np.ravel(atl03Data.z[bb])
#    classification = np.ravel(atl03Data.classification[bb])
#    lat = np.ravel(atl03Data.lat[bb])
#    lon = np.ravel(atl03Data.lon[bb])
    

#    from icesatIO import createHTMLChart
#    createHTMLChart(ytrack, h_ph, classification, lat, lon,
#                        classification_list = [1,2,3],output_folder = "",
#                        in_file03_name = "ATL03", blank = "Viewer_blank.html")

#    # Call getAtlTruthSwath
#    print('RUNNING getAtlTruthSwath...\n')
#    atlTruthData = getAtlTruthSwath(atl03Data, headerData, rotationData, useExistingTruth, truthSwathDir, buffer, outFilePath, createTruthFile)
#    print('done').h5'
    truthSwathDir = '/LIDAR/server/USERS/eric/ecosystem/Sonoma/'  # Path to existing truth data
    outFilePath = '/LIDAR/server/USERS/eric/ecosystem/'
    
    comps = atl03FilePath.split('/')
    filename = comps[-1].split('.')
    title = filename[0] + "_" + gtNum
    # User options
    trimInfo = 'auto'   # OPTIONS: ('none', 'auto', or 'manual')
                        # none - does not trim data at all
                        # auto - trims ATL03 track to extent of truth bounding region
                        # manual - trims ATL03 track by latitude or time
                            # Example: 'manual,lat,38,39'
                            # Only uses data between latitude 38 and 39 deg
                            # Example: 'manual,time,3,4'
                            # Only uses data between time 3 and 4 seconds
    
    createAtl03LasFile = True    # Option to create output measured ATL03 .las file
    createAtl03KmlFile = False   # Option to create output measured ATL03 .kml file
    createAtl08KmlFile = False    # Option to create output measured ATL08 .kml file
    createAtl03CsvFile = False    # Option to create output measured ATL03 .csv file
    createAtl08CsvFile = False    # Option to create output measured ATL08 .csv file
        
    ##### End Inputs for getAtlMeasuredSwath
    
    
    ##### Start Inputs for getAtlTruthSwath
    
    buffer = 50                 # Distance in cross-track (meters) around ATL03 track to look for truth data 
    useExistingTruth = False     # Option to use existing truth data if it exists
    truthSwathDir = '/LIDAR/server/USERS/eric/ecosystem/Sonoma'  # Path to existing truth data
    createTruthFile = True      # Option to create output truth .las file
    
    ##### End Inputs for getAtlTruthSwath
    
    
    ##### CODE BELOW -- DO NOT EDIT ###############################################
    

    # Call getAtlTruthSwath
    print('RUNNING getAtlTruthSwath...\n')
    atlTruthData = getAtlTruthSwath(atl03Data, headerData, rotationData, useExistingTruth, truthSwathDir, buffer, outFilePath, createTruthFile)
    print('done')
    
        
    
    atlTruthData.classification[atlTruthData.classification == 3] = 4
    atlTruthData.classification[atlTruthData.classification == 5] = 4
    # Clip Measured Data to Truth Data Alongtrack
    maxAT = np.nanmax(atlTruthData.alongTrack)
    minAT = np.nanmin(atlTruthData.alongTrack)
    measATfilt = np.where((atl03Data.alongTrack < maxAT) & (atl03Data.alongTrack > minAT))
    atl03Data.alongTrack = np.c_[atl03Data.alongTrack[measATfilt]]
    atl03Data.crossTrack = np.c_[atl03Data.crossTrack[measATfilt]]
    atl03Data.z = np.c_[atl03Data.z[measATfilt]]
    atl03Data.classification = np.c_[atl03Data.classification[measATfilt]]
    atl03Data.northing = np.c_[atl03Data.northing[measATfilt]]
    atl03Data.easting = np.c_[atl03Data.easting[measATfilt]]
    atl03Data.lat = np.c_[atl03Data.lat[measATfilt]]
    atl03Data.lon = np.c_[atl03Data.lon[measATfilt]]
    atl03Data.intensity = np.c_[atl03Data.intensity[measATfilt]]
    atl03Data.signalConf = np.c_[atl03Data.signalConf[measATfilt]]
    atl03Data.time= np.c_[atl03Data.time[measATfilt]]
    
    measAT08filt = np.where((atl08Data.alongTrack < maxAT) & (atl08Data.alongTrack > minAT))
    
    atl08Data.lat = np.c_[atl08Data.lat[measAT08filt]]
    atl08Data.lon = np.c_[atl08Data.lon[measAT08filt]]
    atl08Data.easting = atl08Data.easting[measAT08filt]
    atl08Data.northing = atl08Data.northing[measAT08filt]
    atl08Data.crossTrack = np.c_[atl08Data.crossTrack[measAT08filt]]
    atl08Data.alongTrack = np.c_[atl08Data.alongTrack[measAT08filt]]
    atl08Data.maxCanopy = np.c_[atl08Data.maxCanopy[measAT08filt]]
    atl08Data.teBestFit = np.c_[atl08Data.teBestFit[measAT08filt]]
    atl08Data.teMedian = np.c_[atl08Data.teMedian[measAT08filt]]
    atl08Data.time = np.c_[atl08Data.time[measAT08filt]]
    atl08Data.signalConf = np.c_[atl08Data.signalConf[measAT08filt]]
    atl08Data.classification = np.c_[atl08Data.classification[measAT08filt]]
    atl08Data.intensity = np.c_[atl08Data.intensity[measAT08filt]]
    
    atl08Data.landcover = np.c_[atl08Data.landcover[measAT08filt]]
    atl08Data.night_flag = np.c_[atl08Data.night_flag[measAT08filt]]
    atl08Data.solar_elevation = np.c_[atl08Data.solar_elevation[measAT08filt]]
    atl08Data.snowcover = np.c_[atl08Data.snowcover[measAT08filt]]
    atl08Data.urban_flag = np.c_[atl08Data.urban_flag[measAT08filt]]
    atl08Data.terrain_flg = np.c_[atl08Data.terrain_flg[measAT08filt]]
    atl08Data.brightness_flag = np.c_[atl08Data.brightness_flag[measAT08filt]]
    atl08Data.psf_flag = np.c_[atl08Data.psf_flag[measAT08filt]]
    atl08Data.layer_flag = np.c_[atl08Data.layer_flag[measAT08filt]]
    try:
        atl08Data.cloud_flag_asr = np.c_[atl08Data.cloud_flag_asr[measAT08filt]]
    except:
        atl08Data.cloud_flag_asr = atl08Data.cloud_flag_asr
    atl08Data.msw_flag = np.c_[atl08Data.msw_flag[measAT08filt]]
    atl08Data.asr = np.c_[atl08Data.asr[measAT08filt]]
    atl08Data.snr = np.c_[atl08Data.snr[measAT08filt]]
    atl08Data.n_seg_ph = np.c_[atl08Data.n_seg_ph[measAT08filt]]
    try:
        atl08Data.landsat_perc = np.c_[atl08Data.landsat_perc[measAT08filt]]
    except:
        atl08Data.landsat_perc = atl08Data.landsat_perc
    
    desiredAngle = 90
    atl03Data.crossTrack, atl03Data.alongTrack, R_mat, xRotPt, yRotPt, phi = getCoordRotFwd(atl03Data.easting, atl03Data.northing, [], [], [], desiredAngle)
    atl03Data.crossTrack = np.c_[atl03Data.crossTrack]
    atl03Data.alongTrack = np.c_[atl03Data.alongTrack]
    atl03Data.z = np.c_[atl03Data.z]
    atl03Data.classification = np.c_[atl03Data.classification]

    rotationData = atlRotationStruct(R_mat, xRotPt, yRotPt, desiredAngle, phi)
    
    atl08Data.crossTrack, atl08Data.alongTrack, _, _, _, _ = getCoordRotFwd(atl08Data.easting, atl08Data.northing, R_mat, xRotPt, yRotPt, desiredAngle)
    atl08Data.crossTrack = np.c_[atl08Data.crossTrack]
    atl08Data.alongTrack = np.c_[atl08Data.alongTrack]
    
    atlTruthData.crossTrack, atlTruthData.alongTrack, _, _, _, _ = getCoordRotFwd(atlTruthData.easting, atlTruthData.northing, R_mat, xRotPt, yRotPt, desiredAngle)
    atlTruthData.crossTrack = np.c_[atlTruthData.crossTrack]
    atlTruthData.alongTrack = np.c_[atlTruthData.alongTrack]
    # endIf
    
    atl03Data.alongTrack = atl03Data.alongTrack + 0
    atl03Data.crossTrack = atl03Data.crossTrack + 4
    atl03Data.z = atl03Data.z + (-6.59)
    
    
    superTruth, sortedMeasured = superFilter(atl03Data, atlTruthData, xBuf = 7, classCode = [])
    truthgroundclass  = 2
    truthcanopyclass = 4
    unclassedlist = [6, 9, 18]

    measpc, measoc = perfectClassifier(sortedMeasured, superTruth,
                                       ground = [truthgroundclass],
                                       canopy = [truthcanopyclass], 
                                       unclassed = unclassedlist, 
                                       keepsize = True)
    
    
    
#    uy = sortedMeasured.lat[atl03Data.classification == -1]
#    uz = sortedMeasured.z[atl03Data.classification == -1]
#    dy = sortedMeasured.lat[atl03Data.classification == 0]
#    dz = sortedMeasured.z[atl03Data.classification == 0]
#    gy = sortedMeasured.lat[atl03Data.classification == 1]
#    gz = sortedMeasured.z[atl03Data.classification == 1]
#    cy = sortedMeasured.lat[atl03Data.classification == 2]
#    cz = sortedMeasured.z[atl03Data.classification == 2]
#    hy = sortedMeasured.lat[atl03Data.classification == 3]
#    hz = sortedMeasured.z[atl03Data.classification == 3]
    
    measoc_c = np.c_[measoc]
    measpc_c = np.c_[measpc]
    pdy = sortedMeasured.alongTrack[measpc_c == 0]
    pdz = sortedMeasured.z[measpc_c == 0]
    pgy = sortedMeasured.alongTrack[measpc_c == 1]
    pgz = sortedMeasured.z[measpc_c == 1]
    pcy = sortedMeasured.alongTrack[measpc_c == 2]
    pcz = sortedMeasured.z[measpc_c == 2]
    
    tgy = superTruth.alongTrack[superTruth.classification == 2]
    tgz = superTruth.z[superTruth.classification == 2]
#    tcy = superTruth.z[superTruth.classification < 4]
#    tcz = 
    
    vcy = atl08Data.lat[atl08Data.maxCanopy < 10000]
    vcz = atl08Data.maxCanopy[atl08Data.maxCanopy < 10000]
    vgy = atl08Data.lat
    vgz = atl08Data.teBestFit
#    s0y = atl03Data.lat[atl03Data.signalConf == 0]
#    s0z = atl03Data.z[atl03Data.signalConf == 0]    
#    s1y = atl03Data.lat[atl03Data.signalConf == 1]
#    s1z = atl03Data.z[atl03Data.signalConf == 1]
#    s2y = atl03Data.lat[atl03Data.signalConf == 2]
#    s2z = atl03Data.z[atl03Data.signalConf == 2]
#    s3y = atl03Data.lat[atl03Data.signalConf == 3]
#    s3z = atl03Data.z[atl03Data.signalConf == 3]
#    s4y = atl03Data.lat[atl03Data.signalConf == 4]
#    s4z = atl03Data.z[atl03Data.signalConf == 4]
##
##
    import matplotlib.pyplot as plt
##    
    f = plt.figure()
    plt.grid()
    plt.plot(superTruth.alongTrack[::10],superTruth.z[::10],'.',color=[0.6,0.6,0.6])
    plt.plot(tgy[::5],tgz[::5],'.',color=[0.3,0.3,0.3])
    plt.plot(pdy[::2],pdz[::2],'.',color=[0.4,0.4,0.9])
    plt.plot(pgy[::2],pgz[::2],'.',color=[0.6,0.6,0])
    plt.plot(pcy[::2],pcz[::2],'.',color=[0,0.6,0])

##    ptc = plt.plot(tcy[::300],tcz[::300],'.',color=[0.6,0.6,0.6])
##    ptg =plt.plot(tgy[::300],tgz[::300],'.',color=[0.3,0.3,0.3])
##37.74 to 37.82
#    plt.plot(s0y,s0z,'.',color =[0.95, 0.95, 0.95],alpha=0.01,markeredgecolor="None")
#    plt.plot(s1y,s1z,'.',color =[0.3, 0.3, 0.3],alpha=0.7,markeredgecolor="None")
#    plt.plot(s2y,s2z,'.',color =[0.2, 0.2, 0.2],alpha=0.9,markeredgecolor="None",)
#    plt.plot(s3y,s3z,'.',color =[0.1, 0.1, 0.1],markeredgecolor="None",)
##    plt.plot(s4y,s4z,'.',color =[0, 0, 0],markeredgecolor="None",)
#    skip = 10
#    pmu = plt.plot(uy[::skip],uz[::skip],'.',color=[0.7, 0.82, 1])
#    pmd = plt.plot(dy[::skip],dz[::skip],'.',color=[0, 0, 0.78])
#    pmh = plt.plot(hy[::skip],hz[::skip],'.',color=[0.29, 0.6, 0])
#    pmc = plt.plot(cy[::skip],cz[::skip],'.',color=[105/255, 215/255, 0/255])
#    pmg = plt.plot(gy[::skip],gz[::skip],'.',color=[200/255, 110/255, 0/255])
#    pvc = plt.plot(vcy[::],vcz[::],'o',color=[0/255, 90/255, 0/255])
#    pvg = plt.plot(vgy[::],vgz[::],'o',color=[200/255, 0/255, 0/255])
#    pmu = plt.plot(uy[::skip],uz[::skip],'.',color =[0.95, 0.95, 0.95],alpha=0.01,fillstyle = 'full')
#    pmd = plt.plot(dy[::skip],dz[::skip],'.',color =[0.3, 0.3, 0.3],alpha=0.9,fillstyle = 'full')
#    pmh = plt.plot(hy[::skip],hz[::skip],'.',color =[0.2, 0.2, 0.2],fillstyle = 'full')
#    pmc = plt.plot(cy[::skip],cz[::skip],'.',color =[0.1, 0.1, 0.1],fillstyle = 'full')
#    pmg = plt.plot(gy[::skip],gz[::skip],'.',color =[0, 0, 0],fillstyle = 'full')
    plt.xlabel('Latitude (Decimal Degrees)')
    plt.ylabel('Elevation (m)')
    plt.title(title)
    ldg = plt.legend(['Truth Canopy','Truth Ground','ATL03 Unclassed','ATL03 Ground','ATL03 Canopy','ATL08 Canopy','ATL08 Ground'],fontsize=8)
###    plt.legend([ptg,ptc,pmg,pmc,pmh,pmd,pmu],
###               ['Truth Ground','Truth Canopy','ATL08 Ground','ATL08 Canopy','ATL08 High Canopy','ATL08 Unclassified','ATL08 Never Classified'])
###    ax.legend((ptg,ptc,pmg,pmc,pmh,pmd,pmu),
###               ('Truth Ground','Truth Canopy','ATL08 Ground','ATL08 Canopy','ATL08 High Canopy','ATL08 Unclassified','ATL08 Never Classified'))
#    plt.show()    

