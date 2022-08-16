# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 09:02:09 2022

@author: ERICG
"""

import pandas as pd
import scipy.io as sio

def main():
    in_file = ''
    out_file = ''
    
    # Read CSV by 
    df = pd.read(in_file)
        
    # convert DF to dictionary
    dictionary = df.to_dict('list')
    
    sio.savemat(out_file,dictionary)
    
if __name__ == "__main__":
    main()