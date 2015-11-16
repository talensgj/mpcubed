#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics

from core import systematics_dev
from core import sysrem

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

pg1 = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)

with h5py.File('/data2/talens/3mEast/LBtests/camip_June1_iter5.hdf5', 'r') as f:
    
    idx1 = f['data/camtrans/idx'].value
    z = f['data/camtrans/z'].value
    
    idx2 = f['data/intrapix/idx'].value
    a = f['data/intrapix/a'].value
    b = f['data/intrapix/b'].value
    c = f['data/intrapix/c'].value
    d = f['data/intrapix/d'].value
    
z = pg1.put_values_on_grid(z, idx1, np.nan).ravel()
a = pg2.put_values_on_grid(a, idx2, np.nan).ravel()
b = pg2.put_values_on_grid(b, idx2, np.nan).ravel()
c = pg2.put_values_on_grid(c, idx2, np.nan).ravel()
d = pg2.put_values_on_grid(d, idx2, np.nan).ravel()
    
hg = HealpixGrid(8)
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_June1_iter5.hdf5', 'r') as f:
    
    m = f['data/magnitudes/m'].value
    sigma1 = f['data/magnitudes/sigma'].value
    
    idx = f['data/skytrans/idx'].value
    lstseq = f['data/skytrans/lstseq'].value
    s = f['data/skytrans/s'].value
    sigma2 = f['data/skytrans/sigma'].value
    
tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = s
s = tmp
    
tmp = np.full((hg.npix, 15*13500), fill_value=np.nan)
tmp[idx, lstseq] = sigma2
sigma2 = tmp    
    
f = fLCfile('/data2/talens/3mEast/LBtests/June1.hdf5')
    
ascc, ra, dec, nobs, vmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
nobs = nobs.astype('int')

skyidx = hg.find_gridpoint(ra, dec)

for i in range(len(ascc)):

    if ascc[i] not in ['807144']: continue

    jdmid, lstidx, lst, flux, eflux, sky, flag, x, y = f.read_data(['jdmid', 'lstidx', 'lst', 'flux0', 'eflux0', 'sky', 'flag', 'x', 'y'], [ascc[i]], [nobs[i]])
    lstidx = lstidx.astype('int')

    # Build indices.    
    staridx = np.repeat(i, nobs[i])
    
    ha = np.mod(lst*15 - ra[i], 360.)
    camtransidx = pg1.find_gridpoint(ha, np.repeat(dec[i], nobs[i]))
    intrapixidx = pg2.find_gridpoint(ha, np.repeat(dec[i], nobs[i]))
    
    dayidx = np.floor(jdmid).astype('int')
    dayidx = dayidx - 2457175
    
    skytransidx = np.ravel_multi_index((dayidx, lstidx), (15, 13500))
    
    binidx = np.ravel_multi_index((dayidx, lstidx//50), (15, 270))
    
    # Flag bad data.
    here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    jdmid = jdmid[here]
    flux = flux[here]
    eflux = eflux[here]
    x = x[here]
    y = y[here]

    staridx = staridx[here]
    camtransidx = camtransidx[here]
    intrapixidx = intrapixidx[here]
    skytransidx = skytransidx[here]

    dayidx = dayidx[here]
    lstidx = lstidx[here]
    binidx = binidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    # Find the corrections to the magnitude and the error.
    sol = m[i]
    sol = sol + z[camtransidx]
    sol = sol + a[intrapixidx]*np.sin(2*np.pi*x) + b[intrapixidx]*np.cos(2*np.pi*x) + c[intrapixidx]*np.sin(2*np.pi*y) + d[intrapixidx]*np.cos(2*np.pi*y)
    sol = sol + s[skyidx[i], skytransidx]

    mag = mag - sol
    emag = np.sqrt(emag**2 + sigma1[i]**2 + sigma2[skyidx[i], skytransidx]**2)

    # Apply one final correction for the lst dependence.
    weights = 1/emag**2
    t = np.bincount(lstidx, weights*mag)/np.bincount(lstidx, weights)

    mag = mag - t[lstidx]

    # Bin the data.
    count = index_statistics(binidx, None, statistic='count')
    bin_jdmid = index_statistics(binidx, jdmid, statistic='mean')
    bin_mag = index_statistics(binidx, mag, statistic='mean')
    bin_emag = index_statistics(binidx, mag, statistic='std')/np.sqrt(count)
    
    here = (count == 50)
    bin_jdmid = bin_jdmid[here]
    bin_mag = bin_mag[here]
    bin_emag = bin_emag[here]
    
    plt.plot(jdmid, mag, '.')
    plt.errorbar(bin_jdmid, bin_mag, yerr=bin_emag, fmt='o')
    plt.show()
    
    

