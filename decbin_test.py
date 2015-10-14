#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid
from core.index_functions import index_statistics

from core import sysrem
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
    print ascc
    print ra[ascc=='822579']
    exit()
    lst = np.array([])
    y = np.array([])
    flux = np.array([])
    eflux = np.array([])
    sflux = np.array([])
    for i in range(len(ascc)):
        lc = f['data/'+ascc[i]]
        
        lst = np.append(lst, lc['lst'])
        y = np.append(y, lc['y'])
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
flux = flux[here]
eflux = eflux[here]

# Sysrem.
a1, a2, niter, chisq, npoints, npars = sysrem_dev.sysrem(staridx, haidx1, flux, eflux)
print niter, chisq, npoints, npars
plt.plot(a2, '.')

a1, a2, a, b, niter, chisq, npoints, npars = sysrem_dev.intrarem(staridx, haidx1, haidx2, y, flux, eflux)
print niter, chisq, npoints, npars
plt.plot(a2, '.')

plt.show()
    
# Coarse decor.
mag = -2.5*np.log10(flux)
emag = 2.5/np.log(10.)*eflux/flux

a1, a2, sigam2, niter, chisq, npoints, npars = coarse_dev.coarse_decorrelation(staridx, haidx1, mag, emag, use_weights=False)
print niter, chisq, npoints, npars
plt.plot(a2, '.')

a1, a2, a, b, niter, chisq, npoints, npars = coarse_dev.coarse_positions(staridx, haidx1, haidx2, y, mag, emag, use_weights=False)
print niter, chisq, npoints, npars
plt.plot(a2, '.')

plt.show()

a1, a2, sigam2, niter, chisq, npoints, npars = coarse_dev.coarse_decorrelation(staridx, haidx1, mag, emag, use_weights=True)
print niter, chisq, npoints, npars
plt.plot(a2, '.')

a1, a2, a, b, niter, chisq, npoints, npars = coarse_dev.coarse_positions(staridx, haidx1, haidx2, y,  mag, emag, use_weights=True)
print niter, chisq, npoints, npars
plt.plot(a2, '.')

plt.show()
