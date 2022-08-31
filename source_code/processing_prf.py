# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 13:41:02 2022

@author: ERICG
"""

if __name__ == "__main__":
    import os
    # from getAtlTruthSwath_auto import getAtlTruthSwath
    # from getMeasurementError_auto import getMeasurementError, offsetsStruct
    from icesatReader import get_atl03_struct
    from icesatReader import get_atl08_struct
    # from icesatReader import convert_atl03_to_legacy
    from icesatReader import get_atl_alongtrack
    # from icesatReader import get_atl_alongtrack
    from icesatReference import getTruthFilePaths
    from icesatReference import getTruthHeaders
    from icesatReference import reprojectHeaderData
    from icesatReference import findMatchingTruthFiles
    from icesatReference import loadLasFile
    from icesatReference import perfect_classifier
    from icesatReference import make_buffer
    from icesatReference import ace

    from icesatBinner import rebin_atl08
    from icesatBinner import rebin_truth
    from icesatBinner import match_truth_fields

    from getMeasurementError import getMeasurementError

    
    import numpy as np
    import scipy.io as sio
    import pandas as pd
    import random


    # import pandas as pd


    out_folder = 'E:/0_data/is2/prf/' 
    if os.name == 'nt':
        # basepath03 = 'Z:/data/release/002/ATL03_r002/Finland/'
        # basepath08 = 'Z:/data/release/002/ATL08_r002/Finland/'
        in_atl03 = 'E:/0_data/is2/prf/ATL03/'
        in_atl08 = 'E:/0_data/is2/prf/ATL08/'
    else:
        basepath03 = '/laserpewpew/data/release/002/ATL03_r002/Finland/'
        basepath08 = '/laserpewpew/data/release/002/ATL08_r002/Finland/'
        


    # atl03file = 'ATL03_20181118120428_07770103_002_01.h5'
    # atl08file = 'ATL08_20181118120428_07770103_002_01.h5'

    atl03_list = os.listdir(in_atl03)
    atl08_list = os.listdir(in_atl08)
    # atl03file = atl03_list[2]
    # atl08file = atl08_list[2]
    random.shuffle(atl03_list)

    for i in range(0,len(atl03_list)):
            print(atl03_list[i])
            atl03file = atl03_list[i]
            atl08file = 'ATL08' + atl03_list[i].split('ATL03')[1]
            atl03filepath =  os.path.join(in_atl03, atl03file)
            atl08filepath =  os.path.join(in_atl08, atl08file)
            gt_list = ['gt1l','gt2l','gt3l','gt2r','gt3r','gt1r']
            random.shuffle(gt_list)
            for gt in gt_list:

                # Inputs
                # atl03filepath = basepath03 + atl03file
                # atl08filepath = basepath08 + atl08file
                # gt = 'gt1r'
                
                
                header_file_path =\
                    '/LIDAR/server/USERS/eric/1_experiment/Finland_HeaderData.mat'
                    
                kml_bounds_txt1 = '/LIDAR/server/USERS/eric/2_production/kmlBounds.txt'
               
                print('Generate ATL03 Struct')
                atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath, 
                                         epsg = '32618', kml_bounds_txt = kml_bounds_txt1, 
                                         header_file_path = header_file_path)    
                
                # atl03.df = atl03.df[atl03.df['time'] < 12]
            
                ###### Truth Header Section #######
                # If existing truth header, read it, if there is none, create it
                truthSwathDir = 'E:/data/2018spl_2959_6647'
                truthFileType = 'laz'
                truthFilePaths = getTruthFilePaths(truthSwathDir, truthFileType, logFileID=False)
                
                truthHeaderDF = getTruthHeaders(truthFilePaths, truthFileType, logFileID=False)
            
                ##### Reproject Truth Header ######
                epsg_atl = 'epsg:32618'
                truthHeaderNewDF = reprojectHeaderData(truthHeaderDF, epsg_atl)
                
                # Filter atl03.df by size of reproject
                atl03.df = atl03.df[atl03.df.northing > np.min(truthHeaderNewDF.ymin)]
                atl03.df = atl03.df[atl03.df.northing < np.max(truthHeaderNewDF.ymax)]
                atl03.df = atl03.df.reset_index()
                atl03.df, atl03.rotationData = get_atl_alongtrack(atl03.df)
                
                out_folder = 'E:/PRF/atl08_30m'
                res = 30
                out_name = atl03.atlFileName + '_' + gt + '_' + str(res) + '.mat'
                out_file = os.path.normpath(os.path.join(out_folder, out_name))
                if os.path.exists(out_file) == False:
                    if len(np.unique(atl03.df.classification)) > 2:
                    
                        # Get ATL08
                        atl08 = get_atl08_struct(atl08filepath, gt, atl03) 
                        atl08.df = atl08.df[atl08.df.latitude > np.min(atl03.df.lat_ph)]
                        atl08.df = atl08.df[atl08.df.latitude < np.max(atl03.df.lat_ph)]
                        if len(atl08.df) > 1:
                            # Find truth files that intersect ICESat-2 track
                            #### atlMeasuredData = atlMeasuredData.df
                            buffer = 25
                            _, matchingTruthFileInds = findMatchingTruthFiles(truthHeaderNewDF, 
                                                                               atl03.df,  atl03.rotationData, buffer)
                            matchingTruthFiles = np.array(truthFilePaths)[matchingTruthFileInds]
                            
                            if len(matchingTruthFiles) > 0:
                                print(matchingTruthFiles)
                                
                                # truthFilePath = matchingTruthFiles[1]
                                
                                truth_swath = pd.DataFrame()
                                
                                for i in range(0,len(matchingTruthFiles)):
                                    truth_df = loadLasFile(matchingTruthFiles[i], epsg_atl, atl03.rotationData, decimate_n = 5)
                                    truth_df['classification'] = ace(np.array(truth_df.easting), np.array(truth_df.northing), 
                                                 np.array(truth_df.z), np.array(truth_df.classification))
                                    truth_df = make_buffer(atl03, truth_df, buffer)
                                    truth_swath = truth_swath.append(truth_df)
                            
                                # Calculate off-set corrections for ATL03        
                                atlCorrections = getMeasurementError(atl03, truth_swath)
                                # TODO: make atlCorrections not bad
                                # if type(atlCorrections.alongTrack ) != list:
                                
                                
                                    # Apply along-track/cross-track, and height corrections to ATL03    
                                    # atl03.df.alongtrack = atl03.df.alongtrack + atlCorrections.alongTrack 
                                    # atl03.df.crosstrack = atl03.df.crosstrack + atlCorrections.crossTrack 
                                    # atl03.df.h_ph = atl03.df.h_ph + atlCorrections.z 
                                    
                                truth_swath.alongtrack = truth_swath.alongtrack - atlCorrections.alongTrack 
                                truth_swath.crosstrack = truth_swath.crosstrack - atlCorrections.crossTrack 
                                truth_swath.z = truth_swath.z - atlCorrections.z 
                                
                                # # Apply beam width measurement
                                truth_swath = make_buffer(atl03, truth_swath, 5)
                                # truth_swath.z  = truth_swath.z  +
                                # Do Perfect Classifier
                                measpc, measoc = perfect_classifier(atl03, truth_swath,ground = [2],canopy = [3,4,5], 
                                                  unclassed = [1, 6, 7, 18], keepsize = True)
                                
                                
                                res = 30
                                res_field = 'alongtrack'
                                
                                
                                # atl03.df['classification'] = measoc
                                # ATL08 30m
                                atl08_bin = rebin_atl08(atl03, atl08, gt, res, res_field)
                                
                                out_folder = 'E:/PRF/atl08_30m'
                                out_name = atl03.atlFileName + '_' + gt + '_' + str(res) + '.mat'
                                out_file = os.path.normpath(os.path.join(out_folder, out_name))
                            
                                dictionary = {}
                                dictionary['struct'] = atl08_bin.to_dict('list')
                                sio.savemat(out_file,dictionary)
                            
                                out_folder = 'E:/PRF/atl08_30m'
                                out_name = atl03.atlFileName + '_' + gt + '_' + str(res) + '.csv'
                                out_file = os.path.normpath(os.path.join(out_folder, out_name))
                                atl08_bin.to_csv(out_file)
                            
                                # Truth 30m
                                out_folder = 'E:/PRF/truth_30m'
                                out_name = 'truth_ ' + atl03.atlFileName.split('ATL03_')[1] + '_' + gt + '_' + str(res) + '.mat'
                                out_file = os.path.join(out_folder, out_name)
                                
                                truth_bin = rebin_truth(atl03, truth_swath, res, res_field)
                                truth_bin = match_truth_fields(truth_bin, atl08_bin)
                                dictionary = {}
                                dictionary['struct'] =truth_bin.to_dict('list')
                                sio.savemat(out_file,dictionary)
                                
                                out_folder = 'E:/PRF/truth_30m'
                                out_name = 'truth_ ' + atl03.atlFileName.split('ATL03_')[1] + '_' + gt + '_' + str(res) + '.csv'
                                out_file = os.path.normpath(os.path.join(out_folder, out_name))
                                truth_bin.to_csv(out_file)
                                
                                # ATL08 PC 30m
                                out_folder = 'E:/PRF/atl08_30m_pc'
                                out_name = atl03.atlFileName + '_' + gt + '_' + str(res) + '_pc.mat'
                                out_file = os.path.join(out_folder, out_name)
                                
                                if len(atl03.df['classification']) == len(measpc[0:len(atl03.df)]):
                                    atl03.df['classification'] = measpc[0:len(atl03.df)]
                                    atl08_bin_pc = rebin_atl08(atl03, atl08, gt, res, res_field)
                                    dictionary = {}
                                    dictionary['struct'] =atl08_bin_pc.to_dict('list')
                                    sio.savemat(out_file,dictionary)
                                    
                                    out_folder = 'E:/PRF/atl08_30m_pc'
                                    out_name = atl03.atlFileName + '_' + gt + '_' + str(res) + '_pc.csv'
                                    out_file = os.path.join(out_folder, out_name)
                                    atl08_bin.to_csv(out_file)
