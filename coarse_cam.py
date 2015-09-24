#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from coordinate_grids import PolarGrid
from coarse_decor import coarse_positions

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

def read_data(filename, ascc, nobs):
    
    nstars = len(ascc)
    ndata = np.sum(nobs)
    select = np.append(0, np.cumsum(nobs))
    
    lst = np.zeros(ndata)
    y = np.zeros(ndata)
    flux0 = np.zeros(ndata)
    eflux0 = np.zeros(ndata)
    sky = np.zeros(ndata)
    flags = np.zeros(ndata)
    
    with h5py.File(filename, 'r') as f:
    
        lc = f['data']
    
        for i in range(nstars):
            
            lst[select[i]:select[i+1]] = lc[ascc[i]]['lst']
            y[select[i]:select[i+1]] = lc[ascc[i]]['y']
            flux0[select[i]:select[i+1]] = lc[ascc[i]]['flux0']
            eflux0[select[i]:select[i+1]] = lc[ascc[i]]['eflux0']
            sky[select[i]:select[i+1]] = lc[ascc[i]]['sky']
            flags[select[i]:select[i+1]] = lc[ascc[i]]['flag']

    return lst, y, flux0, eflux0, sky, flags

filename = '/data2/talens/Jul2015/fLC_20150710LPC.hdf5'

with h5py.File(filename, 'r') as f:
    
    ascc1 = f['header_table/ascc'].value
    ra1 = f['header_table/ra'].value
    dec1 = f['header_table/dec'].value
    nobs1 = f['header_table/nobs'].value.astype('int')
    
pg1 = PolarGrid(13500, 720)   
pg2 = PolarGrid(270, 720)

decidx, decuni = pg1.find_decidx(dec1, compact=True)

for ind in range(len(decidx)):
    here = (decuni == ind)
    
    ascc = ascc1[here]
    ra = ra1[here]
    dec = dec1[here]
    nobs = nobs1[here]
    
    lst, y, flux, eflux, sky, flags = read_data(filename, ascc, nobs)
    
    ha = np.mod(lst*15.-np.repeat(ra, nobs), 360)
    
    haidx1 = pg1.find_raidx(ha)
    haidx2 = pg2.find_raidx(ha)
    staridx = np.repeat(np.arange(len(ascc)), nobs)
    
    here = (flux > 0) & (eflux > 0) & (sky > 0) & (flags < 1)
    staridx = staridx[here]
    haidx1 = haidx1[here]
    haidx2 = haidx2[here]
    y = y[here]
    flux = flux[here]
    eflux = eflux[here]
    
    haidx_cam, hauni_cam = np.unique(haidx1, return_inverse=True)
    pointcount_cam = np.bincount(hauni_cam)
    
    haidx_ipx, hauni_ipx = np.unique(haidx2, return_inverse=True)
    pointcount_ipx = np.bincount(hauni_ipx)
    
    mag = -2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux
    
    m, z, a, b = coarse_positions(staridx, hauni_cam, hauni_ipx, y, mag, emag, verbose=True)
    
    with h5py.File('/data2/talens/Jul2015/coarsecam.hdf5') as f:
        grp = f.create_group('data/%i'%decidx[ind])
        grp.create_dataset('haidx_cam', data=haidx_cam)
        grp.create_dataset('pointcount_cam', data=pointcount_cam)
        grp.create_dataset('camtrans', data=z)
        grp.create_dataset('haidx_ipx', data=haidx_ipx)
        grp.create_dataset('pointcount_ipx', data=pointcount_ipx)
        grp.create_dataset('a', data=a)
        grp.create_dataset('b', data=b)
        
        
        
