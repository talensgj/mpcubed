#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from sysrem import sysrem
from coordinate_grids import HealpixGrid

from index_functions import index_statistics

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

hg = HealpixGrid(8)

with h5py.File('/data2/talens/Jul2015/fLC_20150714LPC.hdf5', 'r') as f, h5py.File('/data2/talens/Jul2015/red_20150714LPC.hdf5', 'r') as g:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    nobs = f['header_table/nobs'].value.astype('int')
    
    skyidx = hg.find_gridpoint(ra, dec)
    print skyidx[ascc=='807144']

    here = (skyidx == 266)
    ascc = ascc[here]
    ra = ra[here]
    dec = dec[here]
    nobs = nobs[here]
    
    lstidx = np.array([])
    flux = np.array([])
    eflux = np.array([])
    for i in range(len(ascc)):
        lc = f['data/'+ascc[i]]
        rc = g['data/'+ascc[i]]
        lstidx = np.append(lstidx, lc['lstidx'])
        flux = np.append(flux, rc['ipcflux0'])
        eflux = np.append(eflux, rc['sipcflux0'])
        #sflux = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std', keeplength=True)
        #eflux = np.append(eflux, sflux)
    
    lstidx = lstidx.astype('int')
    staridx = np.repeat(np.arange(len(ascc)), nobs)
    
    here = (flux>0)&(eflux>0)
    staridx = staridx[here]
    lstidx = lstidx[here]
    flux = flux[here]
    eflux = eflux[here]
    
    F, S = sysrem(staridx, lstidx, flux, eflux, verbose=True, eps=1e-6)[:2]
    
    F *= np.nanmedian(S)
    S /= np.nanmedian(S)
   
    array1 = np.full((len(ascc), 13500), fill_value=np.nan)
    array1[staridx, lstidx] = flux/(F[staridx])
    array1 = array1[np.argsort(ra)]
    
    array2 = np.full((len(ascc), 13500), fill_value=np.nan)
    array2[staridx, lstidx] = flux/(F[staridx]*S[lstidx])
    array2 = array2[np.argsort(ra)]
    
    plt.figure(figsize=(16,8))
    ax = plt.subplot(221)
    plt.imshow(array1, aspect='auto', cmap=viridis, vmin=.5, vmax=1.5)
    cb = plt.colorbar()
    cb.set_label('flux0/F')
    plt.xlabel('LST [idx]')
    plt.subplot(222, sharex=ax)
    plt.plot(S, '.', c='k')
    plt.ylim(0.5,1.5)
    plt.xlabel('LST [idx]')
    plt.ylabel('S')
    plt.subplot(223, sharex=ax, sharey=ax)
    plt.imshow(array2, aspect='auto', cmap=viridis, vmin=.5, vmax=1.5)
    cb = plt.colorbar()
    plt.xlabel('LST [idx]')
    cb.set_label('flux0/(F*T)')
    plt.xlim(np.amin(lstidx)-.5, np.amax(lstidx)+.5)
    plt.subplot(224)
    plt.hist(array2[np.isfinite(array2)], bins=np.linspace(0,2,200), histtype='step', lw=2, normed=True)
    plt.xlabel('flux0/(F*T)')
    plt.ylabel('PDF')
    plt.tight_layout()
    plt.show()
    plt.close()
    
    
    
    

        
        
        
