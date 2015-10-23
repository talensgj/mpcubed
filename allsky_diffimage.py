#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid, CartesianGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/3mEast/variations/cam_20150611LPE+flags+ipx+weights.hdf5', 'r') as f:
    
    grp = f['data']
    camtransidx = grp['camtransidx'].value
    z = grp['z'].value
    
pg = PolarGrid(13500, 720)
z1 = pg.put_values_on_grid(z, camtransidx, np.nan)

with h5py.File('/data2/talens/3mEast/variations/cam_20150611LPE+flags+ipxcartesian+weights.hdf5', 'r') as f:
    
    grp = f['data']
    camtransidx = grp['camtransidx'].value
    z = grp['z'].value
    
pg = PolarGrid(13500, 720)
z2 = pg.put_values_on_grid(z, camtransidx, np.nan)

delta = z1 - z2

delta = delta - np.nanmean(delta, axis=0)

vmin = np.nanpercentile(delta, 1)
vmax = np.nanpercentile(delta, 99)

xlim, ylim = np.where(np.isfinite(delta))

plt.figure(figsize=(16,8))

plt.subplot(111)
plt.imshow(delta.T, aspect='auto', cmap=viridis, vmin=vmin, vmax=vmax)
plt.colorbar().set_label(r'z')
plt.xlabel('Hour Angle')
plt.ylabel('Declination')
plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)

plt.tight_layout()
plt.show()
