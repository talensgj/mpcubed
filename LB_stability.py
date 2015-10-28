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

with h5py.File('/data2/talens/3mEast/LBtests/sky_15day.hdf5', 'r') as f:
    s = f['data/s'].value

f = fLCfile('/data2/talens/3mEast/LBtests/15day.hdf5')
pg = PolarGrid(13500, 720)
hg = HealpixGrid(8)

Mascc, Mra, Mdec, Mnobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])
Mnobs = Mnobs.astype('int')

Mstaridx = np.arange(len(Mascc))
decidx, decuni = pg.find_decidx(Mdec, compact=True)
nbins = len(decidx)      

Mskyidx = hg.find_gridpoint(Mra, Mdec)

m = np.zeros(len(Mascc))
z = np.full((13502*722,), fill_value=np.nan)

niter = np.zeros(nbins)
chisq = np.zeros(nbins)
npoints = np.zeros(nbins)
npars = np.zeros(nbins)

for ind in range(nbins):
        
    here = (decuni == ind)
    ascc = Mascc[here]
    ra = Mra[here]
    dec = Mdec[here]
    nobs = Mnobs[here]
    staridx = Mstaridx[here]
    skyidx = Mskyidx[here]

    jdmid, lstidx, lst, flux, eflux, sky, flag = f.read_data(['jdmid', 'lstidx', 'lst', 'flux0', 'eflux0', 'sky', 'flag'], ascc, nobs)
    lstidx = lstidx.astype('int')

    # Build indices.    
    staridx = np.repeat(staridx, nobs)
    
    ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
    camtransidx = pg.find_gridpoint(ha, np.repeat(dec, nobs))
    
    skyidx = np.repeat(skyidx, nobs)
    dayidx = np.floor(jdmid).astype('int')
    dayidx = dayidx - 2457175
    skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))
        
    # Flag bad data.
    here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    flux = flux[here]
    eflux = eflux[here]

    staridx = staridx[here]
    camtransidx = camtransidx[here]
    skyidx = skyidx[here]
    skytransidx = skytransidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    # Get unique indices.
    staridx, staruni = np.unique(staridx, return_inverse=True)
    camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)

    sol = s[skyidx, skytransidx]
    mag = mag - sol

    # Calculate a model fit to the data.
    m[staridx], z[camtransidx], niter[ind], chisq[ind], npoints[ind], npars[ind] = systematics_dev.trans(staruni, camtransuni, mag, emag, verbose=True, use_weights=False)
        
with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter1.hdf5') as f:
    
    hdr = f.create_group('header')
    hdr.create_dataset('decidx', data=decidx)
    hdr.create_dataset('niter', data=niter)
    hdr.create_dataset('chisq', data=chisq)
    hdr.create_dataset('npoints', data=npoints)
    hdr.create_dataset('npars', data=npars)
    
    grp = f.create_group('data')
    grp.create_dataset('ascc', data=Mascc)
    grp.create_dataset('m', data=m)
    dset2 = grp.create_dataset('z', data=z)
    dset2.attrs['grid'] = 'polar'
    dset2.attrs['nx'] = 13500
    dset2.attrs['ny'] = 720    
        
        
        
