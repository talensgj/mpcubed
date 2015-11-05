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

# Initialize reader and coordinate grids.
f = fLCfile('/data2/talens/3mEast/LBtests/June2.hdf5')
pg = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)
hg = HealpixGrid(8)

# Read skymap.
with h5py.File('/data2/talens/3mEast/LBtests/skyip_June2.hdf5', 'r') as g:

    idx = g['data/skytrans/idx'].value
    lstseq = g['data/skytrans/lstseq'].value
    s = g['data/skytrans/s'].value
    
    #sigma1 = g['data/magnitudes/sigma'].value
    #sigma2 = g['data/skytrans/sigma'].value

tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s
s = tmp

#tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
#tmp[idx, lstseq] = sigma2
#sigma2 = tmp

# Read header data.
Mascc, Mra, Mdec, Mnobs, Mvmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
Mnobs = Mnobs.astype('int')

# Build indices for the stars.
Mstaridx = np.arange(len(Mascc))
decidx, decuni = pg.find_decidx(Mdec, compact=True)
Mskyidx = hg.find_gridpoint(Mra, Mdec)

# Initialize arrays.
nbins = len(decidx) 
niter = np.zeros(nbins)
chisq = np.zeros(nbins)
npoints = np.zeros(nbins)
npars = np.zeros(nbins)

m = np.full(len(Mascc), fill_value=np.nan)
z = np.full((13502*722,), fill_value=np.nan)
a = np.full((272*722,), fill_value=np.nan)
b = np.full((272*722,), fill_value=np.nan)
c = np.full((272*722,), fill_value=np.nan)
d = np.full((272*722,), fill_value=np.nan)

for ind in range(nbins):
        
    here = (decuni == ind)
    ascc = Mascc[here]
    ra = Mra[here]
    dec = Mdec[here]
    nobs = Mnobs[here]
    staridx = Mstaridx[here]
    skyidx = Mskyidx[here]

    jdmid, lstidx, lst, flux, eflux, sky, flag, x, y = f.read_data(['jdmid', 'lstidx', 'lst', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y'], ascc, nobs)
    lstidx = lstidx.astype('int')

    # Build indices.    
    staridx = np.repeat(staridx, nobs)
    
    ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
    camtransidx = pg.find_gridpoint(ha, np.repeat(dec, nobs))
    
    intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec, nobs))
    
    skyidx = np.repeat(skyidx, nobs)
    dayidx = np.floor(jdmid).astype('int')
    #dayidx = dayidx - 2457175 #June1
    dayidx = dayidx - 2457190
 
    skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))
        
    # Flag bad data.
    here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    flux = flux[here]
    eflux = eflux[here]
    x = x[here]
    y = y[here]

    staridx = staridx[here]
    camtransidx = camtransidx[here]
    intrapixidx = intrapixidx[here]
    skyidx = skyidx[here]
    skytransidx = skytransidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    sol = s[skyidx, skytransidx]
    mag = mag - sol
    #emag = np.sqrt(emag**2 + (sigma1[staridx])**2 + (sigma2[skyidx, skytransidx])**2)

    # Get unique indices.
    staridx, staruni = np.unique(staridx, return_inverse=True)
    camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)
    intrapixidx, intrapixuni = np.unique(intrapixidx, return_inverse=True)

    # Calculate a model fit to the data.
    m[staridx], z[camtransidx], a[intrapixidx], b[intrapixidx], c[intrapixidx], d[intrapixidx], niter[ind], chisq[ind], npoints[ind], npars[ind] = systematics_dev.trans_ipx(staruni, camtransuni, intrapixuni, mag, emag, x, y, verbose=True, use_weights=False)
        
    offset = np.nanmedian(m[staridx] - Mvmag[staridx])
    m[staridx] = m[staridx] - offset
    z[camtransidx] = z[camtransidx] + offset
        
with h5py.File('/data2/talens/3mEast/LBtests/camip_June2_iter1.hdf5') as f:
    
    hdr = f.create_group('header')
    hdr.create_dataset('decidx', data=decidx)
    hdr.create_dataset('niter', data=niter)
    hdr.create_dataset('chisq', data=chisq)
    hdr.create_dataset('npoints', data=npoints)
    hdr.create_dataset('npars', data=npars)
    
    grp = f.create_group('data')
    
    grp.create_dataset('magnitudes/ascc', data=Mascc)
    grp.create_dataset('magnitudes/m', data=m)
    
    idx, = np.where(~np.isnan(z))
    grp.create_dataset('camtrans/idx', data=idx)
    grp.create_dataset('camtrans/z', data=z[idx])
    
    grp['camtrans'].attrs['grid'] = 'polar'
    grp['camtrans'].attrs['nx'] = 13500
    grp['camtrans'].attrs['ny'] = 720   
    
    idx, = np.where(~np.isnan(a))
    grp.create_dataset('intrapix/idx', data=idx)
    grp.create_dataset('intrapix/a', data=a[idx])
    grp.create_dataset('intrapix/b', data=b[idx])
    grp.create_dataset('intrapix/c', data=c[idx])
    grp.create_dataset('intrapix/d', data=d[idx])
    
    grp['intrapix'].attrs['grid'] = 'polar'
    grp['intrapix'].attrs['nx'] = 270
    grp['intrapix'].attrs['ny'] = 720
