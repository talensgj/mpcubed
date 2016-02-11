#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

from package.statistics.statistics import idxstats

import matplotlib.pyplot as plt

def quarter_variation():
    
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
    
    return 
    
def colorcolor():
    
    with h5py.File('/data2/talens/2015Q2/LPE/sys0_201506BLPE.hdf5', 'r') as f:
        grp = f['data/magnitudes']
        vmag = grp['vmag'].value
        m = grp['mag'].value
    
    with h5py.File('/data2/talens/2015Q2/LPE/fLC_201506BLPE.hdf5', 'r') as f:
        grp = f['header_table']
        bmag = grp['bmag'].value
        
    print np.where(np.isnan(m))
        
    plt.plot(bmag - vmag, m - vmag, '.', alpha=.2)
    plt.xlabel('B - V')
    plt.ylabel('M - V')
    plt.show()
    
    return
    
    
def main():
    
    colorcolor()
    
    return
    
if __name__ == '__main__':
    main()
    
    
    
    
    
    
    



    
