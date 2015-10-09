#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
import TransitModels

filelist = glob.glob('/data2/talens/3mEast/fLC_2015060?LPE.hdf5')
filelist = np.sort(filelist)

plt.figure(figsize=(16,8))

for filename in filelist:
    with h5py.File(filename, 'r') as f:
    
        try:
            lc = f['data/807144'].value
        except:
            continue
    
    ax = plt.subplot(211)
    plt.title('ASCC 807144')
    plt.plot(lc['jdmid'], lc['flux0'], '.')
    plt.plot(lc['jdmid'], lc['flux1'], '.')
    plt.ylabel('Flux')
    plt.subplot(212, sharex=ax)
    #plt.plot(lc['jdmid'], lc['sky'], '.')
    plt.plot(lc['jdmid'], TransitModels.MA_Model(lc['jdmid'], 2454037.612, 2.21857312, .0995*(1.138/0.805), 215.1*0.03142/0.805, 0.0041, 85.51*np.pi/180., 90.*np.pi/180., 0.,0.))
    phase = (lc['jdmid']-2454037.612)/2.21857312
    n = np.unique(np.floor(phase))
    for i in n:
        plt.axvline(2454037.612+i*2.21857312)
    plt.xlabel('Time [JD]')
    plt.ylabel('Sky')
    
plt.tight_layout()
plt.show()
plt.close()

filelist2 = glob.glob('/data2/talens/3mEast/red_2015060?LPE.hdf5')
filelist2 = np.sort(filelist2)

plt.figure(figsize=(16,8))

for filename, filename2 in zip(filelist, filelist2):
    print filename, filename2
    with h5py.File(filename, 'r') as f, h5py.File(filename2, 'r') as g:
    
        try:
            lc = f['data/807144'].value
        except:
            continue
        rc = g['data2/807144'].value
    
    ax = plt.subplot(211)
    plt.title('ASCC 807144')
    plt.plot(lc['jdmid'], rc['sipc_mag0'], '.')
    plt.ylabel('Flux')
    plt.subplot(212, sharex=ax)
    #plt.plot(lc['jdmid'], lc['sky'], '.')
    plt.plot(lc['jdmid'], -(rc['sipc_mag0']-np.nanmean(rc['sipc_mag0'])), '.')
    bin_flux = np.bincount(lc['lstidx'].astype('int')//50, rc['sipc_mag0'])/np.bincount(lc['lstidx'].astype('int')//50)
    bin_jd = np.bincount(lc['lstidx'].astype('int')//50, lc['jdmid'])/np.bincount(lc['lstidx'].astype('int')//50)
    plt.plot(bin_jd, -(bin_flux-np.nanmean(rc['sipc_mag0'])), 'o')
    plt.plot(lc['jdmid'], TransitModels.MA_Model(lc['jdmid'], 2454037.612, 2.21857312, .0995*(1.138/0.805), 215.1*0.03142/0.805, 0.0041, 85.51*np.pi/180., 90.*np.pi/180., 0.,0.)-1)
    phase = (lc['jdmid']-2454037.612)/2.21857312
    n = np.unique(np.floor(phase))
    for i in n:
        plt.axvline(2454037.612+i*2.21857312)
    plt.xlabel('Time [JD]')
    plt.ylabel('Sky')
    
plt.tight_layout()
plt.show()
plt.close()
