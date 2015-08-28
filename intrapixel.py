#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt

from index_functions import index_statistics

fLC = '/data2/talens/Jul2015/fLC_20150716LPC.hdf5'
red = '/data2/talens/Jul2015/red_20150716LPC.hdf5'

from scipy.optimize import minimize
from scipy.ndimage.filters import gaussian_filter1d

def intrapixel1(pars, x, y=None, yerr=None):
    
    tmp = np.cos(2*np.pi*(x+pars[0]))
    
    fit = (pars[1]*tmp+1)
    
    if y is None:
        return fit
    elif yerr is None:
        return np.sum(np.abs((y-fit)))
    else:
        return np.sum((y-fit)**2/yerr**2)

def intrapixel2(pars, x, y=None, yerr=None):
    
    tmp = x+pars[0]
    tmp = np.abs(tmp-np.floor(tmp)-.5)
    tmp = tmp*4-1
    
    fit = (pars[1]*tmp+1)
    
    if y is None:
        return fit
    elif yerr is None:
        return np.sum(np.abs((y-fit)))
    else:
        return np.sum((y-fit)**2/yerr**2)


with h5py.File(fLC, 'r') as f, h5py.File(red, 'r') as g:
    
    ascc = f['table_header/ascc'].value
    dec = f['table_header/dec'].value

    pixa = np.array([])
    deca = np.array([])
    pars1 = np.array([])
    pars2 = np.array([])

    for j in range(len(ascc)):
        
        sid = ascc[j]
        
        lc = f['data/'+sid]
        rc = g['data/'+sid]
        
        x = lc['y']
        gauss = gaussian_filter1d(rc['cflux0'], 100)
        y = rc['cflux0']/gauss
        
        pix = np.unique(np.floor(x))
        
        for i in range(len(pix)):
            here = np.floor(x) == pix[i]
        
            if np.sum(here) < 25: continue
        
            res = minimize(intrapixel1, [0.5, 0.5], args=(x[here],y[here]))
            pars1 = np.append(pars1, res.x[0])
            pars2 = np.append(pars2, res.x[0])
            pixa = np.append(pixa, pix[i])
            deca = np.append(deca, dec[j])

np.savez('intrapixel.npz', pix=pixa, dec=dec, pars1=pars1, pars2=pars2)

plt.scatter(pixa, dec, c=pars1)
plt.colorbar()
plt.show()
