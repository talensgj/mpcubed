#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from package.statistics.statistics import idxstats

import matplotlib.pyplot as plt

path = '/data2/talens/inj_signals/reference/'

with h5py.File(path + 'fLC_201506ALPE.hdf5', 'r') as f:
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value

arg, = np.where(ascc == '344250')
args, = np.where((np.abs(dec - dec[arg]) < 2) & (np.abs(ra - ra[arg]) < 2))
close = ascc[args] 

filelist = glob.glob(path + 'sys0_*.hdf5')
filelist = np.sort(filelist)

ASCC = np.array([])
MAG = np.array([])

for filename in filelist:
    with h5py.File(filename, 'r') as f:
        ascc = f['data/magnitudes/ascc'].value
        mag = f['data/magnitudes/mag'].value
        
    ASCC = np.append(ASCC, ascc)
    MAG = np.append(MAG, mag)
    
ascc, idx = np.unique(ASCC, return_inverse=True)
delta = idxstats(idx, MAG, statistic=np.ptp)

for i in range(len(ascc)):
    if ascc[i] in close:
        print delta[i]

count = np.bincount(idx)
select = (count == 6)

plt.hist(delta[select], bins = np.linspace(0, .25, 51))
plt.xlabel('peak-to-peak of magnitude calibration')
plt.show()

    
