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

filelist = glob.glob('/data2/talens/3mEast/fLC_*LPE.hdf5')
filelist = np.sort(filelist)

filelist2 = glob.glob('/data2/talens/3mEast/red_*LPE.hdf5')
filelist2 = np.sort(filelist2)

ndays = len(filelist)
flux = np.full((ndays, 13500), fill_value=np.nan)

ascc = '822579'
for i in range(ndays):
    with h5py.File(filelist[i]) as f, h5py.File(filelist2[i]) as g:
        try: lc = f['data/'+ascc].value
        except: continue
        rc = g['data/'+ascc].value
        ra = f['header/'+ascc]['ra']
    
    lstidx = lc['lstidx'].astype('int')
    flux[i, lstidx] = rc['intrapix0']

ylim, xlim = np.where(np.isfinite(flux))

#flux = flux - np.nanmean(flux, axis=1, keepdims=True)
#flux = flux - flux[10]

plt.figure(figsize=(16,8))
plt.title('ASCC '+ascc)
plt.imshow(flux, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5)
plt.colorbar().set_label(r'$\Delta m$')
plt.xlabel('LST [idx]')
plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
plt.tight_layout()
plt.show()
