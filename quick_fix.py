#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

def quick_fix(filename):
    
    with h5py.File(filename, 'r+') as myf:
            
        ### Get the sID of the stars to drop
        todrop = []
        for k in myf['data'].iterkeys():
            try:
                sIDsize = myf['data'][k].size
            except Exception:
                pass
            if sIDsize < 51:
                del myf['data'][k]
                del myf['header'][k]
                todrop.append(k)
                
        ### Delete those entries from the header table? assuming it exists!
        the_sids = np.copy(myf['header_table']['ascc'])
        to_extract = [(the_sids == t).nonzero()[0] for t in todrop]

        for k in myf['header_table'].keys():
            tmp = np.copy(myf['header_table'][k])
            tmp_ = np.delete(tmp, to_extract)
            ### and now??? how to replace it???
            del myf['header_table'][k]
            myf['header_table'].create_dataset(k, data=tmp_)
        
filelist = glob.glob('/data2/talens/3mEast/fLC_201508*LPE.hdf5')
filelist = np.sort(filelist)

for filename in filelist:
    print 'Processing', filename
    quick_fix(filename)
