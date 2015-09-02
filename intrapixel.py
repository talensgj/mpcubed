#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt

from coordinate_grids import PolarGrid
from index_functions import index_statistics

fLC = '/data2/talens/Jul2015/fLC_20150716LPC.hdf5'
red = '/data2/talens/Jul2015/red_20150716LPC.hdf5'

from scipy.optimize import minimize
from scipy.ndimage.filters import gaussian_filter1d

from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

def solver(lstidx, y, flux, eflux):
    
    ind = ((lstidx+25)//50).astype('int')
    
    count = np.bincount(ind)
    idx, = np.where(count>0)
    
    cs = np.cos(2*np.pi*y)
    sn = np.sin(2*np.pi*y)
    w = 1/eflux**2
     
    C = np.bincount(ind, flux*w)/np.bincount(ind, w)
    a = C/10*np.cos(.5)*np.ones(len(count))
    b = C/10*np.sin(.5)*np.zeros(len(count))
    
    for i in range(50):

        C = np.bincount(ind, (flux-a[ind]*cs-b[ind]*sn)*w)/np.bincount(ind, w)
        a = np.bincount(ind, (flux-b[ind]*sn-C[ind])*cs*w)/np.bincount(ind, cs**2*w)
        b = np.bincount(ind, (flux-a[ind]*cs-C[ind])*sn*w)/np.bincount(ind, sn**2*w)
        
    #plt.subplot(211)
    #plt.plot(flux, '.')
    #plt.plot(a[ind]*cs+b[ind]*sn+C[ind], '.')
    #plt.subplot(212)
    #plt.plot(flux/(a[ind]*cs+b[ind]*sn+C[ind]), '.')
    #plt.show()
        
    return a[idx], b[idx], C[idx]
    
pg = PolarGrid(13500, 720)

with h5py.File(fLC, 'r') as f, h5py.File(red, 'r') as g:
    
    ascc = f['table_header/ascc'].value
    ra = f['table_header/ra'].value
    dec = f['table_header/dec'].value
    vmag = f['table_header/vmag'].value

    decidx, decuni = pg.find_decidx(dec, compact=True)

    ascc = ascc[decuni==82]
    ra = ra[decuni==82]
    dec = dec[decuni==82]
    vmag = vmag[decuni==82]

    plt.figure(figsize=(16,8))
    for j in range(len(ascc)):
        
        sid = ascc[j]
        
        lc = f['data/'+sid]
        rc = g['data/'+sid]
        
        idx = ((lc['lstidx']+25)//50).astype('int')
        ha = np.mod(lc['lst']*15.-ra[j], 360.)
        ha = np.mod(ha-180, 360)
        ha = np.bincount(idx, ha)/np.bincount(idx)
        ha = ha[np.unique(idx)]
        
        a, b, C = solver(lc['lstidx'], lc['y'], lc['flux0'], lc['eflux0'])
        
        plt.subplot(211)
        plt.plot(ha, np.arctan(b/a), '.', c='k')
        plt.ylabel('Phase')
        plt.ylim(0, np.pi)
        
        plt.subplot(212)
        plt.plot(ha, np.sqrt(a**2+b**2)/C, '.', c='k')
        plt.ylabel('Amplitude')
        
        plt.xlabel('HA [deg]')
    
    plt.tight_layout()
    plt.show()
        
      
        
