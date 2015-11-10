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
f = fLCfile('/data2/talens/3mEast/LBtests/June1.hdf5')
pg = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)
hg = HealpixGrid(8)

# Read transmap.
with h5py.File('/data2/talens/3mEast/LBtests/camip_June1.hdf5', 'r') as g:
    idx1 = g['data/camtrans/idx'].value
    z = g['data/camtrans/z'].value
    idx2 = g['data/intrapix/idx'].value
    a = g['data/intrapix/a'].value
    b = g['data/intrapix/b'].value
    c = g['data/intrapix/c'].value
    d = g['data/intrapix/d'].value

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

# Read skymap.
with h5py.File('/data2/talens/3mEast/LBtests/skyip_June1.hdf5', 'r') as g:
    idx = g['data/skytrans/idx'].value
    lstseq = g['data/skytrans/lstseq'].value
    s = g['data/skytrans/s'].value

tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s
s = tmp

# Read header data.
Mascc, Mra, Mdec, Mnobs, Mvmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
Mnobs = Mnobs.astype('int')

# Build indices for the stars.
Mstaridx = np.arange(len(Mascc))
Mskyidx = hg.find_gridpoint(Mra, Mdec)

for ind in range(len(Mascc)):
        
    here = (Mstaridx == ind)
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
    
    dayidx = np.floor(jdmid).astype('int')
    dayidx = dayidx - 2457175 #June1
    #dayidx = dayidx - 2457190
    
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
    
    lstidx = lstidx[here]
    dayidx = dayidx[here]
    staridx = staridx[here]
    skytransidx = skytransidx[here]
    camtransidx = camtransidx[here]
    intrapixidx = intrapixidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    sol = z[camtransidx] + a[intrapixidx]*np.sin(2*np.pi*x) + b[intrapixidx]*np.cos(2*np.pi*x) + c[intrapixidx]*np.sin(2*np.pi*y) + d[intrapixidx]*np.cos(2*np.pi*y)
    mag = mag - sol

    sol = s[skyidx, skytransidx]
    mag = mag - sol

    # Get unique indices.
    staridx, staruni = np.unique(staridx, return_inverse=True)
    skytransidx, skytransuni = np.unique(skytransidx, return_inverse=True)

    # Calculate a model fit to the data.
    m, t, niter, chisq, npoints, npars = systematics_dev.trans(dayidx, lstidx//50, mag, emag, verbose=True, use_weights=False)
    
    plt.plot(t, '.')
    plt.show()
    plt.close()
    
    print t[~np.isnan(t)]
    
    
#with h5py.File('/data2/talens/3mEast/LBtests/skyip_June2_iter5.hdf5') as f:
    
    #hdr = f.create_group('header')
    #hdr.create_dataset('skyidx', data=skyidx)
    #hdr.create_dataset('niter', data=niter) 
    #hdr.create_dataset('chisq', data=chisq)
    #hdr.create_dataset('npoints', data=npoints)
    #hdr.create_dataset('npars', data=npars)
    
    #grp = f.create_group('data')
    
    #grp.create_dataset('magnitudes/ascc', data=Mascc)
    #grp.create_dataset('magnitudes/m', data=m)
    #grp.create_dataset('magnitudes/sigma', data=sigma1)
    
    #idx, lstseq = np.where(~np.isnan(s))
    #grp.create_dataset('skytrans/idx', data=idx)
    #grp.create_dataset('skytrans/lstseq', data=lstseq)
    #grp.create_dataset('skytrans/s', data=s[idx, lstseq])
    #grp.create_dataset('skytrans/sigma', data=sigma2[idx, lstseq])
    
    #grp['skytrans'].attrs['grid'] = 'healpix'
    #grp['skytrans'].attrs['nx'] = 8
    #grp['skytrans'].attrs['lstmin'] = 0
    #grp['skytrans'].attrs['lstlen'] = 15*13500
