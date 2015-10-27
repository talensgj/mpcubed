#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

filelist = glob.glob('/data2/talens/3mEast/fLC_201506??LPE.hdf5')
filelist = np.sort(filelist)

array = np.full((31, 13500), fill_value=np.nan)

for i in range(len(filelist)):
    with h5py.File(filelist[i], 'r') as f:
        
        try:
            lc = f['/data/807144'].value
        except:
            pass
        else:
            array[i, lc['lstidx'].astype('int')] = lc['flux0']

plt.imshow(array, aspect='auto')
plt.show()
