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
    
camtrans1 = np.full((722, 13502), fill_value=np.nan)
amplitude1 = np.full((722, 272), fill_value=np.nan)
with h5py.File('/data1/talens/aperweights/camip_20150601LPE.hdf5') as f:
    
    decidx = f['header/decidx'].value
    for ind in decidx:
        try: trans = f['data/%i'%ind]
        except: pass
        
        haidx = trans['haidx_cam'].value
        camtrans1[ind, haidx] = trans['camtrans'].value
        haidx = trans['haidx_ipx'].value
        amplitude1[ind, haidx] = np.sqrt(trans['a'].value**2 + trans['b'].value**2)

camtrans2 = np.full((722, 13502), fill_value=np.nan)
amplitude2 = np.full((722, 272), fill_value=np.nan)
with h5py.File('/data1/talens/aperweights/camip_20150830LPE.hdf5') as f:
    
    decidx = f['header/decidx'].value
    for ind in decidx:
        try: trans = f['data/%i'%ind]
        except: pass
        
        haidx = trans['haidx_cam'].value
        camtrans2[ind, haidx] = trans['camtrans'].value
        haidx = trans['haidx_ipx'].value
        amplitude2[ind, haidx] = np.sqrt(trans['a'].value**2 + trans['b'].value**2)
        

camtrans1 = camtrans1[1:-1,1:-1]
camtrans2 = camtrans2[1:-1,1:-1]

camtrans1 = camtrans1 - np.nanmean(camtrans1, axis=1, keepdims=True)
camtrans2 = camtrans2 - np.nanmean(camtrans2, axis=1, keepdims=True)

diff1 = camtrans1-camtrans2
diff1 = diff1 - np.nanmean(diff1, axis=1, keepdims=True)

ylim1, xlim1 = np.where(np.isfinite(diff1))

amplitude1 = amplitude1[1:-1,1:-1]
amplitude2 = amplitude2[1:-1,1:-1]

diff2 = amplitude1 - amplitude2

ylim2, xlim2 = np.where(np.isfinite(diff2))


ax = plt.subplot(231)
plt.imshow(camtrans1, aspect='auto', vmin=-.5, vmax=.5, cmap=viridis)
plt.colorbar()

plt.subplot(232, sharex=ax, sharey=ax)
plt.imshow(camtrans2, aspect='auto', vmin=-.5, vmax=.5, cmap=viridis)
plt.colorbar()

plt.subplot(233, sharex=ax, sharey=ax)
plt.imshow(diff1, aspect='auto', vmin=-.1, vmax=.1, cmap=viridis)
plt.colorbar()
plt.xlim(np.amin(xlim1)-.5, np.amax(xlim1)+.5)
plt.ylim(np.amin(ylim1)-.5, np.amax(ylim1)+.5)

ax2 = plt.subplot(234)
plt.imshow(amplitude1, aspect='auto', vmin=0, vmax=.05, cmap=viridis)
plt.colorbar()

plt.subplot(235, sharex=ax2, sharey=ax2)
plt.imshow(amplitude2, aspect='auto', vmin=0, vmax=.05, cmap=viridis)
plt.colorbar()

plt.subplot(236, sharex=ax2, sharey=ax2)
plt.imshow(diff2, aspect='auto', vmin=-.01, vmax=.01, cmap=viridis)
plt.colorbar()
plt.xlim(np.amin(xlim2)-.5, np.amax(xlim2)+.5)
plt.ylim(np.amin(ylim2)-.5, np.amax(ylim2)+.5)

plt.show()

    
