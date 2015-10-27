#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics

from core import systematics_dev
from time import time
from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

f = fLCfile('/data2/talens/3mEast/fLC_20150716LPE.hdf5')
ascc, ra, dec, nobs, vmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
lstidx, lst, x, y, flux, eflux, sky, flag = f.read_data(['lstidx', 'lst', 'x', 'y', 'flux0', 'eflux0', 'sky', 'flag'], ascc, nobs)

nobs = nobs.astype('int')
lstidx = lstidx.astype('int')

# Build indices
staridx = np.repeat(np.arange(len(ascc)), nobs)

ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
pg = PolarGrid(13500, 720)
camtransidx = pg.find_gridpoint(ha, np.repeat(dec, nobs))

pg2 = PolarGrid(270, 720)
intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec, nobs))

hg = HealpixGrid(8)
skyidx = hg.find_gridpoint(ra, dec)
skyidx = np.repeat(skyidx, nobs)
skytransidx = np.ravel_multi_index((skyidx, lstidx), (hg.npix, 13500))

# Flag bad data.
here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
x = x[here]
y = y[here]
flux = flux[here]
eflux = eflux[here]

staridx = staridx[here]
camtransidx = camtransidx[here]
intrapixidx = intrapixidx[here]
skytransidx = skytransidx[here]

# Convert flux to magnitudes
mag = 25 - 2.5*np.log10(flux)
emag = 2.5/np.log(10)*eflux/flux

# Get unique indices.
staridx, staruni = np.unique(staridx, return_inverse=True)
camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)
intrapixidx, intrapixuni = np.unique(intrapixidx, return_inverse=True)
skytransidx, skytransuni = np.unique(skytransidx, return_inverse=True)

# Calculate a model fit to the data.
sol2 = 0.
for i in range(5):
    m1, z, a, b, c, d, niter, chisq, npoints, npars = systematics_dev.trans_ipx(staruni, camtransuni, intrapixuni, mag-sol2, emag, x, y, verbose=True, use_weights=False, maxiter=20)
    sol1 = z[camtransuni] + a[intrapixuni]*np.sin(2*np.pi*x) + b[intrapixuni]*np.cos(2*np.pi*x) + c[intrapixuni]*np.sin(2*np.pi*y) + d[intrapixuni]*np.cos(2*np.pi*y)
    #m1, z, niter, chisq, npoints, nparts = systematics_dev.trans(staruni, camtransuni, mag-sol2, emag, verbose=True, use_weights=False, maxiter=20)
    #sol1 = z[camtransuni]
     
    m2, s, niter, chisq, npoints, nparts = systematics_dev.trans(staruni, skytransuni, mag-sol1, emag, verbose=True, use_weights=False, maxiter=20)
    sol2 = s[skytransuni]
    
    array = pg.put_values_on_grid(z, camtransidx, np.nan)
    
    plt.imshow(array.T, aspect='auto', cmap=viridis)
    plt.colorbar()
    plt.show()
    plt.close()
    
    plt.plot(vmag[staridx], m1, '.')
    plt.show()
    plt.close()
    
    plt.plot(m1, m1-m2, '.', alpha=.2)
    plt.show()
    plt.close()
    
    if (i > 0):
        tmp = np.abs(m1-m1_old)
        print 'm'
        print np.max(tmp), np.percentile(tmp, 99), np.percentile(tmp, 95)
        
        tmp = np.abs(z-z_old)
        print 'z'
        print np.max(tmp), np.percentile(tmp, 99), np.percentile(tmp, 95)
        
        tmp = np.abs(s-s_old)
        print 's'
        print np.max(tmp), np.percentile(tmp, 99), np.percentile(tmp, 95)
    
        array = pg.put_values_on_grid(z-z_old, camtransidx, np.nan)
    
        plt.imshow(array.T, aspect='auto', cmap=viridis)
        plt.colorbar()
        plt.show()
        plt.close()
    
    m1_old = np.copy(m1)
    z_old = np.copy(z)
    s_old = np.copy(s)
