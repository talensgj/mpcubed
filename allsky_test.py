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
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

# Read data.
f = fLCfile('/data2/talens/3mEast/fLC_20150611LPE.hdf5')
ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])
lstidx, lst, x, y, flux, eflux, sky, flag = f.read_data(['lstidx', 'lst', 'x', 'y', 'flux0', 'eflux0', 'sky', 'flag'], ascc, nobs)

lstidx = lstidx.astype('int')
nobs = nobs.astype('int')

# Build indices.
staridx = np.repeat(np.arange(len(ascc)), nobs)

ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
pg = PolarGrid(13500, 720)
haidx, hauni = pg.find_gridpoint(ha, np.repeat(dec, nobs), compact=True)

pg2 = PolarGrid(270, 720)
posidx, posuni = pg2.find_gridpoint(ha, np.repeat(dec, nobs), compact=True)

#cg = CartesianGrid(320, 220, margin=45)
#posidx, posuni = cg.find_gridpoint(x, y, compact=True)

hg = HealpixGrid(8)
skyidx = hg.find_gridpoint(ra, dec)
skyidx = np.repeat(skyidx, nobs)
skyidx = np.ravel_multi_index((skyidx, lstidx), (hg.npix, 13500))
skyidx, skyuni = np.unique(skyidx, return_inverse=True)

# Flag bad data.
here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
lst = lst[here]
staridx = staridx[here]
hauni = hauni[here]
skyuni = skyuni[here]
posuni = posuni[here]
x = x[here]
y = y[here]
flux = flux[here]
eflux = eflux[here]

# Convert flux to magnitudes.
mag = -2.5*np.log10(flux)
emag = 2.5/np.log(10)*eflux/flux

# Calculate a model fit to the data.
m, z, a, b, c, d, niter, chisq, npoints, npars = systematics_dev.trans_ipx(staridx, hauni, posuni, mag, emag, x, y, verbose=True, use_weights=False)
print niter, chisq, npoints, npars
    
z = pg.put_values_on_grid(z, ind=haidx, fill_value=np.nan)
z = z - np.nanmean(z, axis=0, keepdims=True)
z = z - np.nanmean(z)

xlim, ylim = np.where(np.isfinite(z))

plt.imshow(z.T, aspect='auto', vmin=-.5, vmax=1, cmap=viridis)
plt.xlim(np.amin(xlim)-.5, np.amax(xlim)+.5)
plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)
plt.xlabel('Hour Angle [idx]')
plt.ylabel('Declination [idx]')
plt.colorbar()
plt.show()
 
Ax = np.sqrt(a**2 + b**2)
Ax = cg.put_values_on_grid(Ax, posidx, fill_value=np.nan)
Ax = Ax[1:-1,1:-1]

Ay = np.sqrt(c**2 + d**2)
Ay = cg.put_values_on_grid(Ay, posidx, fill_value=np.nan)
Ay = Ay[1:-1,1:-1]

plt.subplot(221)
plt.imshow(Ax.T, vmin=0, vmax=.1, cmap=viridis, extent=(45,4008-45,45,2672-45))
plt.colorbar()
plt.xlim(0,4008)
plt.ylim(0,2672)

plt.subplot(222)
plt.imshow(Ay.T, vmin=0, vmax=.1, cmap=viridis, extent=(45,4008-45,45,2672-45))
plt.colorbar()
plt.xlim(0,4008)
plt.ylim(0,2672)

plt.subplot(223)
plt.imshow((Ax + Ay).T, vmin=0, vmax=.1, cmap=viridis, extent=(45,4008-45,45,2672-45))
plt.colorbar()
plt.xlim(0,4008)
plt.ylim(0,2672)

plt.show()


 
