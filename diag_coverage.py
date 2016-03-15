#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob

import h5py
import numpy as np

filelist = glob.glob('/data2/mascara/LaPalma/2015????LP?/fLC/fLC_*LP?.hdf5')
filelist = np.sort(filelist)

ascc = np.array([])
for filename in filelist:
    
    with h5py.File(filename, 'r') as f:
        
        try:
            ascc_ = f['header_table/ascc'].value
        except:
            pass
        else:
            ascc = np.append(ascc, ascc_)
        
print len(ascc), len(np.unique(ascc))
