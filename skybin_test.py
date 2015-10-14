#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core import sysrem
from core import coarse_decor
from core.coordinate_grids import HealpixGrid

from core.index_functions import index_statistics

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

hg = HealpixGrid(8)


with h5py.File('/data2/talens/3mEast/fLC_20150611LPE.hdf5', 'r') as f, h5py.File('/data2/talens/3mEast/red_20150611LPE.hdf5', 'r') as g:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    nobs = f['header_table/nobs'].value.astype('int')
    
    skyidx = hg.find_gridpoint(ra, dec)
    print skyidx[ascc=='807144']
    print np.amin(dec), np.amax(dec)
    here = (skyidx == 266)
    ascc = ascc[here]
    ra = ra[here]
    dec = dec[here]
    nobs = nobs[here]
    
    print np.amin(dec), np.amax(dec)
    exit()
    jdmid = np.array([])
    lstidx = np.array([])
    mag = np.array([])
    emag = np.array([])
    for i in range(len(ascc)):
        
        lc = f['data/'+ascc[i]]
        rc = g['data/'+ascc[i]]
        
        jdmid = np.append(jdmid, lc['jdmid'])
        lstidx = np.append(lstidx, lc['lstidx'])
        mag = np.append(mag, rc['ipcflux0'])
        emag = np.append(emag, 2.5/np.log(10)*lc['eflux0']/lc['flux0'])
    
    lstidx = lstidx.astype('int')
    staridx = np.repeat(np.arange(len(ascc)), nobs)
    
    here = np.isfinite(mag) & (emag > 0)
    staridx = staridx[here]
    jdmid = jdmid[here]
    lstidx = lstidx[here]
    mag = mag[here]
    emag = emag[here]
    
    lstuni, lstidx = np.unique(lstidx, return_inverse=True)
    staruni, staridx = np.unique(staridx, return_inverse=True)
    
    m, z, sigma1, sigma2 = coarse_decorrelation(staridx, lstidx, mag, emag, verbose=True, maxiter=25)
    
    print sigma1
    print sigma2
    
    plt.hist(sigma2, bins=np.linspace(0,2,100))
    plt.show()
    
    fit = m[staridx]+z[lstidx]
    
    for i in range(len(ascc)):
        here = staridx == i
        
        plt.subplot(211)
        plt.plot(jdmid[here], mag[here], '.')
        plt.plot(jdmid[here], fit[here], '.')
        plt.subplot(212)
        
        bin_jd = index_statistics(lstidx[here]//50, jdmid[here], statistic='mean')
        bin_mag = index_statistics(lstidx[here]//50, mag[here]-fit[here], statistic='mean')
        
        #plt.plot(jdmid[here], mag[here]-fit[here], '.')
        plt.plot(bin_jd, bin_mag, '.')
        plt.show()
