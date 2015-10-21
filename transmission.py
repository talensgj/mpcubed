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

filelist = glob.glob('/data2/talens/3mEast/red_aper_xy/camip_201506??LPE.hdf5')
filelist = np.sort(filelist)

ndays = len(filelist)
camtrans = np.full((ndays, 13502), fill_value=np.nan)
Ay = np.full((ndays, 272), fill_value=np.nan)
Ax = np.full((ndays, 272), fill_value=np.nan)

decidx = 451
for i in range(ndays):
    with h5py.File(filelist[i]) as f:
        try: tc = f['data/%i'%decidx]
        except: continue

        haidx = tc['haidx_cam'].value
        camtrans[i, haidx] = tc['camtrans'].value
        haidx = tc['haidx_ipx'].value
        a = tc['a'].value
        b = tc['b'].value
        Ay[i, haidx] = np.sqrt(a**2+b**2)
        c = tc['c'].value
        d = tc['d'].value
        Ax[i, haidx] = np.sqrt(c**2+d**2)


camtrans = camtrans[:,1:-1]
Ax = Ax[:,1:-1]
Ay = Ay[:,1:-1]

ylim, xlim = np.where(np.isfinite(camtrans))

camtrans = camtrans - np.nanmean(camtrans, axis=1, keepdims=True)
camtrans = camtrans - camtrans[10]

plt.figure(figsize=(16,8))

ax = plt.subplot(311)
plt.title('Declination idx %i'%decidx)
plt.imshow(camtrans, aspect='auto', cmap=viridis, vmin=-.1, vmax=.1, extent=(0,24,0,ndays))
plt.colorbar().set_label(r'$\Delta m$')

plt.subplot(312, sharex=ax, sharey=ax)
plt.imshow(Ax, aspect='auto', cmap=viridis, vmin=0, vmax=.1, extent=(0,24,0,ndays))
plt.colorbar().set_label(r'$A_{x}$')

plt.subplot(313, sharex=ax, sharey=ax)
plt.imshow(Ay, aspect='auto', cmap=viridis, vmin=0, vmax=.1, extent=(0,24,0,ndays))
plt.colorbar().set_label(r'$A_{y}$')

plt.xlim(18.5, 23)

plt.xlabel('Hour Angle')

plt.tight_layout()
plt.show()
