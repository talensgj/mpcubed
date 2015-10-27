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
nobs = nobs.astype('int')

pg = PolarGrid(13500, 720)
decidx = pg.find_decidx(dec)

pg2 = PolarGrid(270, 720)

hg = HealpixGrid(8)
skyidx = hg.find_gridpoint(ra, dec)

sol2 = 0.


for i in range(5):
    arr1 = 0.
    arr2 = 0.
    arr3 = 0.
    arr4 = 0.
    arr5 = 0.
    print 'Camera'
    for ind in np.unique(decidx):
        print ind
        here = (decidx == ind)
        lstidx, lst, x, y, flux, eflux, sky, flag = f.read_data(['lstidx', 'lst', 'x', 'y', 'flux0', 'eflux0', 'sky', 'flag'], ascc[here], nobs[here])
        lstidx = lstidx.astype('int')

        # Build indices
        staridx = np.repeat(np.arange(len(ascc[here])), nobs[here])

        ha = np.mod(lst*15 - np.repeat(ra[here], nobs[here]), 360.)
        camtransidx = pg.find_gridpoint(ha, np.repeat(dec[here], nobs[here]))

        intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec[here], nobs[here]))

        longsky = np.repeat(skyidx[here], nobs[here])

        # Flag bad data.
        here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
        x = x[here]
        y = y[here]
        flux = flux[here]
        eflux = eflux[here]

        staridx = staridx[here]
        camtransidx = camtransidx[here]
        intrapixidx = intrapixidx[here]
        longsky = longsky[here]
        lstidx = lstidx[here]

        # Convert flux to magnitudes
        mag = 25 - 2.5*np.log10(flux)
        emag = 2.5/np.log(10)*eflux/flux

        if i>0:
            sol2 = s[longsky, lstidx]

        # Get unique indices.
        staridx, staruni = np.unique(staridx, return_inverse=True)
        camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)
        intrapixidx, intrapixuni = np.unique(intrapixidx, return_inverse=True)

        m, z, a, b, c, d, niter, chisq, npoints, npars = systematics_dev.trans_ipx(staruni, camtransuni, intrapixuni, mag - sol2, emag, x, y, verbose=False, use_weights=False, maxiter=20)
        
        arr1 = arr1 + pg.put_values_on_grid(z, camtransidx, 0)
        arr2 = arr2 + pg2.put_values_on_grid(a, intrapixidx, 0)
        arr3 = arr3 + pg2.put_values_on_grid(a, intrapixidx, 0)
        arr4 = arr4 + pg2.put_values_on_grid(a, intrapixidx, 0)
        arr5 = arr5 + pg2.put_values_on_grid(a, intrapixidx, 0)
        
    z = np.ravel(arr1)
    a = np.ravel(arr2)
    b = np.ravel(arr3)
    c = np.ravel(arr4)
    d = np.ravel(arr5)

    arr1 = np.full((hg.npix, 13500), fill_value=np.nan)
    print 'Sky'
    for ind in np.unique(skyidx):
        print ind
        here = (skyidx == ind)
        lstidx, lst, x, y, flux, eflux, sky, flag = f.read_data(['lstidx', 'lst', 'x', 'y', 'flux0', 'eflux0', 'sky', 'flag'], ascc[here], nobs[here])
        lstidx = lstidx.astype('int')

        # Build indices
        staridx = np.repeat(np.arange(len(ascc[here])), nobs[here])

        ha = np.mod(lst*15 - np.repeat(ra[here], nobs[here]), 360.)
        camtransidx = pg.find_gridpoint(ha, np.repeat(dec[here], nobs[here]))

        intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec[here], nobs[here]))
        
        # Flag bad data.
        here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
        x = x[here]
        y = y[here]
        flux = flux[here]
        eflux = eflux[here]

        staridx = staridx[here]
        camtransidx = camtransidx[here]
        intrapixidx = intrapixidx[here]
        lstidx = lstidx[here]
        
        sol1 = z[camtransidx] - a[intrapixidx]*np.sin(2*np.pi*x)  - b[intrapixidx]*np.cos(2*np.pi*x)  - c[intrapixidx]*np.sin(2*np.pi*y)  - d[intrapixidx]*np.cos(2*np.pi*y)
        
        # Convert flux to magnitudes
        mag = 25 - 2.5*np.log10(flux)
        emag = 2.5/np.log(10)*eflux/flux

        # Get unique indices.
        staridx, staruni = np.unique(staridx, return_inverse=True)
        camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)
        intrapixidx, intrapixuni = np.unique(intrapixidx, return_inverse=True)
        lstidx, lstuni = np.unique(lstidx, return_inverse=True)

        m, s, niter, chisq, npoints, nparts = systematics_dev.trans(staruni, lstuni, mag - sol1, emag, verbose=True, use_weights=False, maxiter=20)
        arr1[ind, lstidx] = s

    s = arr1

exit()
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
