#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid
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

with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_test_norm.hdf5', 'r') as f:
    idx1 = f['data/camtrans/idx'].value
    z = f['data/camtrans/z'].value
    idx2 = f['data/intrapix/idx'].value
    a = f['data/intrapix/a'].value
    b = f['data/intrapix/b'].value
    c = f['data/intrapix/c'].value
    d = f['data/intrapix/d'].value

# Initialize reader and coordinate grids.
f = fLCfile('/data2/talens/3mEast/LBtests/15day.hdf5')
pg = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)
hg = HealpixGrid(8)

z = pg.put_values_on_grid(z, idx1, np.nan)
a = pg2.put_values_on_grid(a, idx2, np.nan)
b = pg2.put_values_on_grid(b, idx2, np.nan)
c = pg2.put_values_on_grid(c, idx2, np.nan)
d = pg2.put_values_on_grid(d, idx2, np.nan)

z = np.ravel(z)
a = np.ravel(a)
b = np.ravel(b)
c = np.ravel(c)
d = np.ravel(d)

# Read header data.
Mascc, Mra, Mdec, Mnobs, Mvmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
Mnobs = Mnobs.astype('int')

# Build indices for the stars.
Mstaridx = np.arange(len(Mascc))
skyidx, skyuni = hg.find_gridpoint(Mra, Mdec, compact=True)

# Initialize arrays.
nbins = len(skyidx)      
niter = np.zeros(nbins)
chisq = np.zeros(nbins)
npoints = np.zeros(nbins)
npars = np.zeros(nbins)

m = np.full(len(Mascc), fill_value=np.nan)
sigma1 = np.full(len(Mascc), fill_value=np.nan)
s = np.full((hg.npix, 15*13500), fill_value=np.nan)
sigma2 = np.full((hg.npix, 15*13500), fill_value=np.nan)

for ind in range(nbins):
        
    here = (skyuni == ind)
    ascc = Mascc[here]
    ra = Mra[here]
    dec = Mdec[here]
    nobs = Mnobs[here]
    staridx = Mstaridx[here]

    jdmid, lstidx, lst, flux, eflux, sky, flag, x, y = f.read_data(['jdmid', 'lstidx', 'lst', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y'], ascc, nobs)
    lstidx = lstidx.astype('int')

    # Build indices.    
    staridx = np.repeat(staridx, nobs)
    
    dayidx = np.floor(jdmid).astype('int')
    dayidx = dayidx - 2457175
    
    skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))
        
    ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
    camtransidx = pg.find_gridpoint(ha, np.repeat(dec, nobs))
    intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec, nobs))
        
    # Flag bad data.
    here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    flux = flux[here]
    eflux = eflux[here]
    x = x[here]
    y = y[here]

    staridx = staridx[here]
    skytransidx = skytransidx[here]
    camtransidx = camtransidx[here]
    intrapixidx = intrapixidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    sol = z[camtransidx] + a[intrapixidx]*np.sin(2*np.pi*x) + b[intrapixidx]*np.cos(2*np.pi*x) + c[intrapixidx]*np.sin(2*np.pi*y) + d[intrapixidx]*np.cos(2*np.pi*y)
    mag = mag - sol

    # Get unique indices.
    staridx, staruni = np.unique(staridx, return_inverse=True)
    skytransidx, skytransuni = np.unique(skytransidx, return_inverse=True)

    # Calculate a model fit to the data.
    m[staridx], s[skyidx[ind], skytransidx], niter[ind], chisq[ind], npoints[ind], npars[ind] = systematics_dev.trans(staruni, skytransuni, mag, emag, verbose=True, use_weights=False)
    #m[staridx], s[skyidx[ind], skytransidx], sigma1[staridx], sigma2[skyidx[ind], skytransidx], niter[ind], chisq[ind], npoints[ind], npars[ind] = systematics_dev.trans(staruni, skytransuni, mag, emag, verbose=True, use_weights=True)
    
    offset = np.nanmedian(m[staridx] - Mvmag[staridx])
    m[staridx] = m[staridx] - offset
    s[skyidx[ind], skytransidx] = s[skyidx[ind], skytransidx] + offset
    print offset
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_test_normout.hdf5') as f:
    
    hdr = f.create_group('header')
    hdr.create_dataset('skyidx', data=skyidx)
    hdr.create_dataset('niter', data=niter) 
    hdr.create_dataset('chisq', data=chisq)
    hdr.create_dataset('npoints', data=npoints)
    hdr.create_dataset('npars', data=npars)
    
    grp = f.create_group('data')
    
    grp.create_dataset('magnitudes/ascc', data=Mascc)
    grp.create_dataset('magnitudes/m', data=m)
    grp.create_dataset('magnitudes/sigma', data=sigma1)
    
    idx, lstseq = np.where(~np.isnan(s))
    grp.create_dataset('skytrans/idx', data=idx)
    grp.create_dataset('skytrans/lstseq', data=lstseq)
    grp.create_dataset('skytrans/s', data=s[idx, lstseq])
    grp.create_dataset('skytrans/sigma', data=sigma2[idx, lstseq])
    
    grp['skytrans'].attrs['grid'] = 'healpix'
    grp['skytrans'].attrs['nx'] = 8
    grp['skytrans'].attrs['lstmin'] = 0
    grp['skytrans'].attrs['lstlen'] = 15*13500
