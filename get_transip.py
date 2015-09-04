#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

from coordinate_grids import PolarGrid, CartesianGrid
from index_functions import index_statistics

from sysrem import sysrem
from intrarem import trans_intrapixel, bin_intrapixel
    
fLC = '/data2/talens/Jul2015/fLC_20150716LPN.hdf5'

with h5py.File(fLC, 'r') as f:
    
    ascc1 = f['table_header/ascc'].value
    ra1 = f['table_header/ra'].value
    dec1 = f['table_header/dec'].value
    nobs1 = f['table_header/nobs'].value.astype('int')
    vmag1 = f['table_header/vmag'].value

pg = PolarGrid(13500, 720)
decidx, decuni = pg.find_decidx(dec1, compact=True)

pg1 = PolarGrid(270, 720)

niter = np.zeros(len(decidx), dtype='int')
chisq = np.zeros(len(decidx), dtype='float')
npoints = np.zeros(len(decidx), dtype='int')
npars = np.zeros(len(decidx), dtype='int')

for ind in range(len(decidx)):

    ascc = ascc1[decuni==ind]
    ra = ra1[decuni==ind]
    dec = dec1[decuni==ind]
    nobs = nobs1[decuni==ind]
    vmag = vmag1[decuni==ind]

    ndata = np.sum(nobs)
    select = np.append(0, np.cumsum(nobs))

    flux = np.zeros(ndata)
    eflux = np.zeros(ndata)
    y = np.zeros(ndata)
    lst = np.zeros(ndata)
    flags = np.zeros(ndata)
    sflux = np.zeros(ndata)
    sky = np.zeros(ndata)

    with h5py.File(fLC, 'r') as f:
        
        for i in range(len(ascc)):
            
            lc = f['data/'+ascc[i]].value
            
            flux[select[i]:select[i+1]] = lc['flux0']
            eflux[select[i]:select[i+1]] = lc['eflux0']
            y[select[i]:select[i+1]] = lc['y']
            lst[select[i]:select[i+1]] = lc['lst']
            flags[select[i]:select[i+1]] = lc['flag']
            sflux[select[i]:select[i+1]] = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std', keeplength=True)
            sky[select[i]:select[i+1]] = lc['sky']
            
    ha = np.mod(lst*15.-np.repeat(ra, nobs), 360.)
    
    haidx1 = pg.find_raidx(ha)
    haidx2 = pg1.find_raidx(ha)
    staridx = np.repeat(np.arange(len(ascc)), nobs)
    
    here = (flux>0)&(sflux>0)&(sky>0)&(flags<1)
    flux = flux[here]
    eflux = eflux[here]
    sflux = sflux[here]
    y = y[here]
    haidx1 = haidx1[here]
    haidx2 = haidx2[here]
    staridx = staridx[here]
    
    if len(flux) == 0 : continue
    
    haidx1, hauni1 = np.unique(haidx1, return_inverse=True)
    haidx2, hauni2 = np.unique(haidx2, return_inverse=True)
    
    F, T, a, b, niter[ind], chisq[ind], npoints[ind], npars[ind] = trans_intrapixel(staridx, hauni1, hauni2, y, flux, sflux, verbose=True)
    
    with h5py.File('/data2/talens/Jul2015/transipN.hdf5') as f:
        
        grp = f.create_group('data/%i'%decidx[ind])
        grp.create_dataset('haidx_cam', data=haidx1)
        grp.create_dataset('camtrans', data=T)
        grp.create_dataset('haidx_ip', data=haidx2)
        grp.create_dataset('a', data=a)
        grp.create_dataset('b', data=b)

with h5py.File('/data2/talens/Jul2015/transipN.hdf5') as f:

    grp = f.create_group('header')

    grp.create_dataset('decidx', data = decidx)
    grp.create_dataset('niter', data = niter)
    grp.create_dataset('npoints', data = npoints)
    grp.create_dataset('npars', data = npars)
    grp.create_dataset('chisq', data = chisq)
