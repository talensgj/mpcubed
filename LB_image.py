#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
        
with h5py.File('/data2/talens/3mEast/LBtests/camip_201506LPE.hdf5', 'r') as f:
    
    grp = f['data']
    camtransidx = grp['camtransidx'].value
    z = grp['z'].value
    
    intrapixidx = grp['intrapixidx'].value
    a = grp['a'].value
    b = grp['b'].value
    c = grp['c'].value
    d = grp['d'].value

camtransidx = camtransidx.astype('int')
intrapixidx = intrapixidx.astype('int')
        
pg = PolarGrid(13500, 720)
z = pg.put_values_on_grid(z, camtransidx, np.nan)
        
pg = PolarGrid(270, 720)
a = pg.put_values_on_grid(a, intrapixidx, np.nan)
b = pg.put_values_on_grid(b, intrapixidx, np.nan)
c = pg.put_values_on_grid(c, intrapixidx, np.nan)
d = pg.put_values_on_grid(d, intrapixidx, np.nan)

vmin = np.nanpercentile(z, 1)
vmax = np.nanpercentile(z, 99)

xlim, ylim = np.where(np.isfinite(z))

for i in np.unique(ylim):
    plt.plot(z[:,i], '.')
    plt.show()

plt.figure(figsize=(16,5))

plt.imshow(z.T, aspect='auto', vmin=vmin, vmax=vmax, cmap=viridis)
plt.colorbar().set_label(r'$\Delta m$')
plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')

plt.tight_layout()
plt.show()

Ax = np.sqrt(a**2 + b**2)
phix = np.arctan2(b, a)
xlim, ylim = np.where(np.isfinite(Ax))

plt.figure(figsize=(16,10))

ax = plt.subplot(211)
plt.imshow(Ax.T, aspect='auto', vmin=0, vmax=.1, cmap=viridis)
plt.colorbar().set_label(r'$A_x$')
plt.ylabel('Declination')

plt.subplot(212, sharex=ax, sharey=ax)
plt.imshow(phix.T, aspect='auto', vmin=-np.pi, vmax=np.pi, cmap=viridis)
plt.colorbar().set_label(r'$\phi_x$')
plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')

plt.tight_layout()
plt.show()

Ay = np.sqrt(c**2 + d**2)
phiy = np.arctan2(d, c)
xlim, ylim = np.where(np.isfinite(Ay))

plt.figure(figsize=(16,10))

ax = plt.subplot(211)
plt.imshow(Ay.T, aspect='auto', vmin=0, vmax=.1, cmap=viridis)
plt.colorbar().set_label(r'$A_y$')
plt.ylabel('Declination')

plt.subplot(212, sharex=ax, sharey=ax)
plt.imshow(phiy.T, aspect='auto', vmin=-np.pi, vmax=np.pi, cmap=viridis)
plt.colorbar().set_label(r'$\phi_y$')
plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')

plt.tight_layout()
plt.show()


