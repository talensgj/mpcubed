#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

filelist = glob.glob('/data2/talens/3mEast/camip_20150???LPE.hdf5')
filelist = np.sort(filelist)

ndays = len(filelist)
camtrans = np.full((ndays, 13502), fill_value=np.nan)
amplitude = np.full((ndays, 272), fill_value=np.nan)

decidx = 451
for i in range(ndays):
    with h5py.File(filelist[i]) as f:
        try: tc = f['data/%i'%decidx]
        except: continue

        haidx = tc['haidx_cam'].value
        camtrans[i, haidx] = tc['camtrans'].value
        haidx = tc['haidx_ipx'].value
        a = tc['c'].value
        b = tc['d'].value
        amplitude[i, haidx] = np.sqrt(a**2+b**2)

camtrans = camtrans[:,1:-1]
amplitude = amplitude[:,1:-1]

ylim, xlim = np.where(np.isfinite(camtrans))

camtrans = camtrans - np.nanmean(camtrans, axis=1, keepdims=True)

plt.figure(figsize=(16,8))

plt.subplot(211)
plt.title('Declination idx %i'%decidx)
plt.imshow(camtrans, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar().set_label(r'$\Delta m$')

plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)

plt.subplot(212)
plt.imshow(amplitude, aspect='auto', cmap=viridis, vmin=0, vmax=.1)
plt.colorbar().set_label('Amplitude')

plt.xlim((np.amin(xlim))/50.-.5, (np.amax(xlim))/50.+.5)

plt.xlabel('Hour Angle [idx]')

plt.tight_layout()
plt.show()
