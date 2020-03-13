#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 12:56:37 2019

@author: eguenther
"""

import numpy as np

def aceReassign(indexc, indexg, zc, cc, zg):
    for r in range(0,len(cc)):
        cindex = np.array(indexc[r])
        gindex = np.array(indexg[r])
        zcl = zc[cindex]
        ccl = cc[cindex]
        if 9 in ccl:
            continue
        zdem = zg[gindex]
        zdelta = zcl - zdem
        check1 = ((zdelta <= 0.5) & (ccl == 1))
        check2 = ((zdelta > 0.5) & (ccl == 1))
        index1 = cindex[check1[:,0]]
        index2 = cindex[check2[:,0]]
        cc[index1] = 2
        cc[index2] = 4
    return cc