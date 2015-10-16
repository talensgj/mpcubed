#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid
from core.index_functions import index_statistics

from core import intrapix_dev
from core import coarse_decor

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

with h5py.File('/data2/talens/3mEast/fLC_20150818LPE.hdf5', 'r') as f:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    nobs = f['header_table/nobs'].value.astype('int')
    vmag = f['header_table/vmag'].value
    
    decidx = pg1.find_decidx(dec)
    print decidx[ascc=='807144']
    
    here = (decidx == 451)
    ascc = ascc[here]
    ra = ra[here]
    dec = dec[here]
    nobs = nobs[here]
    vmag = vmag[here]

    lst = np.array([])
    y = np.array([])
    x = np.array([])
    flux = np.array([])
    eflux = np.array([])
    sflux = np.array([])
    for i in range(len(ascc)):
        lc = f['data/'+ascc[i]]
        
        lst = np.append(lst, lc['lst'])
        y = np.append(y, lc['y'])
        x = np.append(x, lc['x'])
        flux = np.append(flux, lc['flux0'])
        eflux = np.append(eflux, lc['eflux0'])
        tmp = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std', keeplength=True)
        sflux = np.append(sflux, tmp)
    
eflux = sflux
    
ha = np.mod(lst*15.-np.repeat(ra, nobs), 360)

haidx1 = pg1.find_raidx(ha)
haidx2 = pg2.find_raidx(ha)
staridx = np.repeat(np.arange(len(ascc)), nobs)

here = (flux > 0) & (eflux > 0)
staridx = staridx[here]
haidx1 = haidx1[here]
haidx2 = haidx2[here]
y = y[here]
x = x[here]
flux = flux[here]
eflux = eflux[here]

mag = -2.5*np.log10(flux)
emag = 2.5/np.log(10.)*eflux/flux

m, z, a, b, c, d, niter, chisq, npoints, npars = intrapix_dev.coarse_positions(staridx, haidx1, haidx2, x, y, mag, emag, verbose=True, use_weights=False)

plt.subplot(211)
plt.plot(np.sqrt(a**2+b**2), '.')
plt.plot(np.sqrt(c**2+d**2), '.')
plt.subplot(212)
plt.plot(np.arctan2(b, a), '.')
plt.plot(np.arctan2(d, c), '.')

m, z, a, b, c, d, niter, chisq, npoints, npars = coarse_decor.coarse_positions(staridx, haidx1, haidx2, x, y, mag, emag, verbose=True, use_weights=False)

plt.subplot(211)
plt.plot(np.sqrt(a**2+b**2), '.')
plt.plot(np.sqrt(c**2+d**2), '.')
plt.subplot(212)
plt.plot(np.arctan2(b, a), '.')
plt.plot(np.arctan2(d, c), '.')
plt.show()

plt.subplot(211)
plt.plot(mag, '.')
plt.plot(m[staridx] + z[haidx1], '.')
plt.subplot(212)
plt.plot(mag - m[staridx] - z[haidx1], '.')
plt.plot(a[haidx2]*np.sin(2*np.pi*x)+b[haidx2]*np.cos(2*np.pi*x)+c[haidx2]*np.sin(2*np.pi*y)+d[haidx2]*np.cos(2*np.pi*y), '.')
plt.show()

