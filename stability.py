#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
    
    
filelist = glob.glob('/data2/talens/3mEast/camip_20150???LPE.hdf5')
filelist = np.sort(filelist)

amplitude = np.full((len(filelist), 272), fill_value=np.nan)
phase = np.full((len(filelist), 272), fill_value=np.nan)
camtrans = np.full((len(filelist), 13502), fill_value=np.nan)

chisq = np.full(len(filelist), fill_value=np.nan)

for i in range(len(filelist)):
    with h5py.File(filelist[i], 'r') as f:
        
        try:
            haidx = f['data/451/haidx_cam'].value
            camtrans[i, haidx] = f['data/451/camtrans']
            
            haidx = f['data/451/haidx_ipx'].value
            a = f['data/451/a'].value
            b = f['data/451/b'].value
            amplitude[i, haidx] = np.sqrt(a**2+b**2)
            phase[i, haidx] = np.arctan2(b, a)
            
            decidx = f['header/decidx'].value
            arg = decidx==451
            tmp = f['header/chisq'].value
            chisq[i] = tmp[arg]
            
        except:
            pass

args, = np.where(chisq > 1)
camtrans[args] = np.nan
amplitude[args] = np.nan
phase[args] = np.nan

camtrans = camtrans - np.nanmedian(camtrans, axis=1, keepdims=True)

#camtrans = camtrans - camtrans[10]
#amplitude = amplitude - amplitude[10]
#phase = phase - phase[10]

plt.figure(figsize=(16,8))
ax = plt.subplot(311)
plt.imshow(camtrans, aspect='auto', vmin=-1, vmax=1, extent=(0,24,0,len(filelist)))
plt.colorbar().set_label('Transmission')
plt.ylabel('Day')

plt.subplot(312, sharex=ax, sharey=ax)
plt.imshow(amplitude, aspect='auto', vmin=0, vmax=.05, extent=(0,24,0,len(filelist)))
plt.colorbar().set_label('Amplitude')
plt.ylabel('Day')

plt.subplot(313, sharex=ax, sharey=ax)
plt.imshow(phase, aspect='auto', vmin=-np.pi, vmax=np.pi, extent=(0,24,0,len(filelist)))
plt.colorbar().set_label('Phase')

plt.xlim(18.5,23)
plt.xlabel('Hour Angle')
plt.ylabel('Day')

plt.tight_layout()
plt.show()

camtrans = camtrans - camtrans[10]
amplitude = amplitude - amplitude[10]
phase = phase - phase[10]

plt.figure(figsize=(16,8))
ax = plt.subplot(311)
plt.imshow(camtrans, aspect='auto', vmin=-.1, vmax=.1, extent=(0,24,0,len(filelist)))
plt.colorbar().set_label('Transmission')
plt.ylabel('Day')

plt.subplot(312, sharex=ax, sharey=ax)
plt.imshow(amplitude, aspect='auto', vmin=0, vmax=.005, extent=(0,24,0,len(filelist)))
plt.colorbar().set_label('Amplitude')
plt.ylabel('Day')

plt.subplot(313, sharex=ax, sharey=ax)
plt.imshow(phase, aspect='auto', vmin=-np.pi/10, vmax=np.pi/10, extent=(0,24,0,len(filelist)))
plt.colorbar().set_label('Phase')

plt.xlim(18.5,23)
plt.xlabel('Hour Angle')
plt.ylabel('Day')

plt.tight_layout()
plt.show()
