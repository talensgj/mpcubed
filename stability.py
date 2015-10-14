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
    
filelist = glob.glob('/data1/talens/aperweights/camip_20150???LPE.hdf5')
filelist = np.sort(filelist)

chisq = np.full(len(filelist), fill_value=np.nan)
npoints = np.full(len(filelist), fill_value=np.nan)
npars = np.full(len(filelist), fill_value=np.nan)

camtrans = np.full((len(filelist), 13502), fill_value=np.nan)
pointcount = np.full((len(filelist), 13502), fill_value=np.nan)
amplitude = np.full((len(filelist), 272), fill_value=np.nan)

for i in range(len(filelist)):
    
    with h5py.File(filelist[i], 'r') as f:
        try: trans = f['data/451']
        except: pass
        
        haidx = trans['haidx_cam'].value
        camtrans[i, haidx] = trans['camtrans'].value
        pointcount[i, haidx] = trans['pointcount_cam'].value
        haidx = trans['haidx_ipx'].value
        amplitude[i, haidx] = np.sqrt(trans['a'].value**2 + trans['b'].value**2)

        decidx = f['header/decidx'].value
        here = decidx == 451
        chisq[i] = f['header/chisq'][here]

npoints = np.nansum(pointcount, axis=1)

args, = np.where(npoints<75000)
camtrans[args] = np.nan
amplitude[args] = np.nan

#args, = np.where(chisq>1)
#camtrans[args] = np.nan
#amplitude[args] = np.nan

camtrans = camtrans[:,1:-1]
amplitude = amplitude[:,1:-1]

camtrans = camtrans - camtrans[:,11000][:,None]

plt.figure(figsize=(16,8))
ax = plt.subplot(211)
plt.imshow(camtrans, aspect='auto', vmin=-.5, vmax=.5, cmap=viridis, extent=(0,24,0,1))
plt.colorbar().set_label('Transmission')

plt.subplot(212, sharex=ax, sharey=ax)
plt.imshow(pointcount, aspect='auto', cmap=viridis, extent=(0,24,0,1))
plt.colorbar().set_label('Amplitude')

plt.xlim(18.75, 23)
plt.xlabel('Hour Angle')

plt.tight_layout()
plt.show()

camtrans = camtrans - camtrans[10]
amplitude = amplitude - amplitude[10]
camtrans = camtrans - camtrans[:,12000][:,None]

plt.figure(figsize=(16,8))
ax = plt.subplot(211)
plt.imshow(camtrans, aspect='auto', vmin=-.1, vmax=.1, cmap=viridis, extent=(0,24,0,1))
plt.colorbar().set_label('Delta Transmission')

plt.subplot(212, sharex=ax, sharey=ax)
plt.imshow(pointcount, aspect='auto', cmap=viridis, extent=(0,24,0,1))
plt.colorbar().set_label('Delta Transmission')

plt.xlim(18.75, 23)
plt.xlabel('Hour Angle')

plt.tight_layout()
plt.show()


for i in range(len(filelist)):
    
    print filelist[i]
    
    plt.figure(figsize=(16,8))
    ax = plt.subplot(211)
    plt.title(filelist[i])
    plt.plot(camtrans[i], '.')
    plt.ylim(-.1, .1)
    
    plt.subplot(212, sharex=ax)
    plt.imshow(camtrans, aspect='auto', vmin=-.1, vmax=.1, cmap=viridis)
    plt.axhline(i)
    
    plt.xlim(18.75/24.*13500, 23./24.*13500)
    
    plt.show()
    plt.close()
    
    
