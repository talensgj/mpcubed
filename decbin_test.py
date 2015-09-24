#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

import intrarem
import sysrem
from coarse_decor import coarse_positions

from coordinate_grids import PolarGrid

from index_functions import index_statistics

from scipy.optimize import minimize

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

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', 'r') as f:
    
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
    
    ha = np.mod(lst*15.-np.repeat(ra, nobs), 360)
    
    haidx1 = pg1.find_raidx(ha)
    haidx2 = pg2.find_raidx(ha)
    staridx = np.repeat(np.arange(len(ascc)), nobs)
    
    here = (flux>0)&(sflux>0)
    staridx = staridx[here]
    haidx1 = haidx1[here]
    haidx2 = haidx2[here]
    y = y[here]
    flux = flux[here]
    eflux = eflux[here]
    sflux = sflux[here]
    
    mag = -2.5*np.log10(flux)
    smag = 2.5/np.log(10)*sflux/flux
    emag = 2.5/np.log(10)*eflux/flux
    
    F, T, a, b = intrarem.trans_intrapixel(staridx, haidx1, haidx2, y, flux, eflux)[:4]
    fit0 = F[staridx]*T[haidx1]*(a[haidx2]*np.sin(2*np.pi*y)+b[haidx2]*np.cos(2*np.pi*y)+1)
    
    m, z, a, b = coarse_positions(staridx, haidx1, haidx2, y, mag, emag, verbose=True, weights=False)
    fit1 = m[staridx]+z[haidx1]+a[haidx2]*np.sin(2*np.pi*y)+b[haidx2]*np.cos(2*np.pi*y)
    
    m, z, a, b, sigma1 = coarse_positions(staridx, haidx1, haidx2, y, mag, emag, verbose=True)
    fit2 = m[staridx]+z[haidx1]+a[haidx2]*np.sin(2*np.pi*y)+b[haidx2]*np.cos(2*np.pi*y)
    
    for arg in np.unique(staridx):
        
        here = staridx == arg
        
        ax = plt.subplot(211)
        plt.plot(mag[here], '.', c='k')
        plt.plot(-2.5*np.log10(fit0[here]), '.', c='grey')
        plt.plot(fit1[here], '.', c='r')
        plt.plot(fit2[here], '.', c='g')
        plt.subplot(212, sharex=ax)
        plt.plot(-2.5*np.log10(flux[here]/fit0[here]), '.', c='grey')
        plt.plot(mag[here]-fit1[here], '.', c='r')
        plt.plot(mag[here]-fit2[here], '.', c='g')
        plt.show()
        plt.close()
        
        
        
