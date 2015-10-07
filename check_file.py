#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

def quick_fix(filename):
    
    with h5py.File(filename, 'r+') as myf:
            
        #ascc = myf['data'].keys()
            
        for k in myf['data'].iterkeys():
            try:
                sIDsize = myf['data'][k].size
            except Exception:
                pass
            if sIDsize != myf['header'][k]['nobs']:
                print k
                
        
filelist = glob.glob('/data2/talens/3mEast/fLC_20150605LPE.hdf5')
filelist = np.sort(filelist)

for filename in filelist:
    print 'Processing', filename
    quick_fix(filename)
