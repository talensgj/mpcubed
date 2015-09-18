#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

import intrarem
from coordinate_grids import PolarGrid

from index_functions import index_statistics

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
    
    decidx = pg1.find_decidx(dec)
    print decidx[ascc=='807144']
    
    here = (decidx == 451)
    ascc = ascc[here]
    ra = ra[here]
    dec = dec[here]
    nobs = nobs[here]
    
    lst = np.array([])
    y = np.array([])
    flux = np.array([])
    eflux = np.array([])
    for i in range(len(ascc)):
        lc = f['data/'+ascc[i]]
        
        lst = np.append(lst, lc['lst'])
        y = np.append(y, lc['y'])
        flux = np.append(flux, lc['flux0'])
        sflux = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std', keeplength=True)
        #eflux = np.append(eflux, lc['eflux0'])
        eflux = np.append(eflux, sflux)
    
    ha = np.mod(lst*15.-np.repeat(ra, nobs), 360)
    
    haidx1 = pg1.find_raidx(ha)
    haidx2 = pg2.find_raidx(ha)
    staridx = np.repeat(np.arange(len(ascc)), nobs)
    
    here = (flux>0)&(eflux>0)
    staridx = staridx[here]
    haidx1 = haidx1[here]
    haidx2 = haidx2[here]
    y = y[here]
    flux = flux[here]
    eflux = eflux[here]
    
    F1, T1, a1, b1, niter, chisq, npoints, npars = intrarem.trans_intrapixel(staridx, haidx1, haidx2, y, flux, eflux, verbose=True, eps=1e-6)
    fit1 = F1[staridx]*T1[haidx1]*(a1[haidx2]*np.sin(2*np.pi*y)+b1[haidx2]*np.cos(2*np.pi*y)+1)
    
    F1 *= np.nanmax(T1)
    T1 /= np.nanmax(T1)
   
    array1 = np.full((len(ascc), 13502), fill_value=np.nan)
    array1[staridx, haidx1] = flux/(F1[staridx])
    #array1 = array1[:,1:-1]
    array1 = array1[np.argsort(ra)]
    array1 = np.roll(array1, 13500/2, axis=1)
    
    array2 = np.full((len(ascc), 13502), fill_value=np.nan)
    array2[staridx, haidx1] = (a1[haidx2]*np.sin(2*np.pi*y)+b1[haidx2]*np.cos(2*np.pi*y)+1)
    #array2 = array2[:,1:-1]
    array2 = array2[np.argsort(ra)]
    array2 = np.roll(array2, 13500/2, axis=1)
    
    array3 = np.full((len(ascc), 13502), fill_value=np.nan)
    array3[staridx, haidx1] = flux/(F1[staridx]*T1[haidx1])
    #array3 = array3[:,1:-1]
    array3 = array3[np.argsort(ra)]
    array3 = np.roll(array3, 13500/2, axis=1)
    
    count = np.bincount(haidx1)
    here = count>0
    ha1 = pg1.find_ra(np.arange(len(T1)))/360*24
    ha1 = np.mod(ha1-12, 24)-12
    
    array4 = array3/array2
    
    ind = np.arange(len(T1))
    ind = np.mod(ind-13502/2, 13502)
    
    haidx1 = np.mod(haidx1-13502/2, 13502)
    
    plt.figure(figsize=(16,12))
    
    ax = plt.subplot(321)
    plt.imshow(array1, aspect='auto', cmap=viridis, vmin=0, vmax=1)
    cb = plt.colorbar()
    plt.xlabel('HA [idx]')
    cb.set_label('flux0/F')
    
    plt.subplot(322, sharex=ax)
    plt.plot(ind, T1, '.', c='k')
    plt.ylim(0, 1)
    cb = plt.colorbar()
    plt.xlabel('HA [idx]')
    plt.ylabel('T')
    
    plt.subplot(323, sharex=ax, sharey=ax)
    plt.imshow(array3, aspect='auto', cmap=viridis, vmin=.9, vmax=1.1)
    cb = plt.colorbar()
    plt.xlabel('HA [idx]')
    cb.set_label('flux0/(F*T)')
    
    plt.subplot(324, sharex=ax, sharey=ax)
    plt.imshow(array2, aspect='auto', cmap=viridis, vmin=.9, vmax=1.1)
    cb = plt.colorbar()
    plt.xlabel('HA [idx]')
    cb.set_label('IP')
    
    plt.subplot(325, sharex=ax, sharey=ax)
    plt.imshow(array3/array2, aspect='auto', cmap=viridis, vmin=.9, vmax=1.1)
    cb = plt.colorbar()
    plt.xlim(np.amin(haidx1)-.5, np.amax(haidx1)+.5)
    plt.xlabel('HA [idx]')
    cb.set_label('flux0/(F*T*IP)')
    
    plt.subplot(326)
    plt.hist(array4[np.isfinite(array4)], bins=np.linspace(.5, 1.5, 200), histtype='step', lw=2, normed=True)
    plt.xlabel('flux0/(F*T*IP)')
    plt.ylabel('PDF')
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    
    
    
    

        
        
        
