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
import matplotlib.gridspec as gridspec

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

# Read data.
f = fLCfile('/data2/talens/3mEast/fLC_20150611LPE.hdf5')
ascc, ra, dec, nobs, vmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])

with h5py.File('/data2/talens/3mEast/variations/cam_20150611LPE+flags+ipx+weights.hdf5', 'r') as f:
    
    hdr = f['header']
    niter = hdr.attrs['niter']
    chisq = hdr.attrs['chisq']
    npoints = hdr.attrs['npoints']
    npars = hdr.attrs['npars']
    
    grp = f['data']
    staridx = grp['staridx'].value
    m = grp['m'].value
    
    camtransidx = grp['camtransidx'].value
    z = grp['z'].value
    
    #intrapixidx = grp['intrapixidx'].value
    #a = grp['a'].value
    #b = grp['b'].value
    #c = grp['c'].value
    #d = grp['d'].value

pg = PolarGrid(13500, 720)
decidx = pg.find_decidx(dec)
z = pg.put_values_on_grid(z, camtransidx, np.nan)

#pg = PolarGrid(270, 720)
#Ax = pg.put_values_on_grid(np.sqrt(a**2 + b**2), intrapixidx, np.nan)
#Ay = pg.put_values_on_grid(np.sqrt(c**2 + d**2), intrapixidx, np.nan)
#phix = pg.put_values_on_grid(np.arctan2(b, a), intrapixidx, np.nan)
#phiy = pg.put_values_on_grid(np.arctan2(d, c), intrapixidx, np.nan)

#cg = CartesianGrid(250, 167)
#Ax = cg.put_values_on_grid(np.sqrt(a**2 + b**2), intrapixidx, np.nan)
#Ay = cg.put_values_on_grid(np.sqrt(c**2 + d**2), intrapixidx, np.nan)
#phix = cg.put_values_on_grid(np.arctan2(b, a), intrapixidx, np.nan)
#phiy = cg.put_values_on_grid(np.arctan2(d, c), intrapixidx, np.nan)

offsets = np.bincount(decidx[staridx], m - vmag[staridx], minlength=722)/np.bincount(decidx[staridx], minlength=722)

z = z + offsets

vmin = np.nanpercentile(z, 1)
vmax = np.nanpercentile(z, 99)

xlim, ylim = np.where(np.isfinite(z))

plt.figure(figsize=(16,8))

gs = gridspec.GridSpec(1, 3, width_ratios=[15,15,1])

plt.subplot(gs[0])
plt.plot(vmag[staridx], m - offsets[decidx[staridx]], '.', alpha=.2)
plt.plot([2, 8.4], [2, 8.4], c='k')
plt.xlabel('V')
plt.ylabel('m')

plt.subplot(gs[1])
plt.title(r'$\chi^2=%.2f$, niter=%i'%(chisq, niter))
im = plt.imshow(z.T, aspect='auto', cmap=viridis, vmin=vmin, vmax=vmax)
plt.xlabel('Hour Angle')
plt.ylabel('Declination')
plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)

ax = plt.subplot(gs[2])
plt.colorbar(im, cax=ax).set_label(r'z')

plt.tight_layout()
plt.show()

#plt.figure(figsize=(16,8))

#ax = plt.subplot(221)
#plt.imshow(Ax.T, cmap=viridis, vmin=0, vmax=.1)
#plt.colorbar().set_label(r'$A_x$')

#plt.subplot(222, sharex=ax, sharey=ax)
#plt.imshow(Ay.T, cmap=viridis, vmin=0, vmax=.1)
#plt.colorbar().set_label(r'$A_y$')

#plt.subplot(223, sharex=ax, sharey=ax)
#plt.imshow(phix.T, vmin=-np.pi, vmax=np.pi, cmap=viridis)
#plt.colorbar().set_label(r'$\phi_x$')

#plt.subplot(224, sharex=ax, sharey=ax)
#plt.imshow(phiy.T, vmin=-np.pi, vmax=np.pi, cmap=viridis)
#plt.colorbar().set_label(r'$\phi_y$')

#plt.tight_layout()
#plt.show()
