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
    
    decidx = np.searchsorted(np.linspace(-90,90,721), dec)
    
    ascc = ascc[decidx==451]

    for sid in ascc:
        print sid
        lc = f['data/'+sid]
        rc = g['data/'+sid]
        
        jd = lc['jdmid']-np.floor(lc['jdmid'])
        x = lc['y']
        gauss = gaussian_filter1d(rc['cflux0'], 100)
        y = rc['cflux0']/gauss
        
        pix = np.unique(np.floor(x))
        
        fit1 = np.ones(len(y))
        fit2 = np.ones(len(y))
        pars1 = np.zeros((len(pix), 2))
        pars2 = np.zeros((len(pix), 2))
        
        for i in range(len(pix)):
            here = np.floor(x) == pix[i]
        
            if np.sum(here) < 25: continue
        
            res1 = minimize(intrapixel1, [0.5, 0.5], args=(x[here],y[here]))
            pars1[i] = res1.x
            fit1[here] = intrapixel1(res1.x, x[here])
            
            res2 = minimize(intrapixel2, [0.5, 0.5], args=(x[here],y[here]))
            pars2[i] = res2.x
            fit2[here] = intrapixel2(res2.x, x[here])
        
        #plt.plot(jd, y, '.')
        #plt.plot(jd, fit1, '.')
        #plt.plot(jd, fit2, '.')
        #plt.show()
        
        ax1 = plt.subplot(221)
        plt.title('Cosine')
        plt.plot(pix, pars1[:,0]%.5, '.', c='k')
        plt.ylabel('phase')
        ax2 = plt.subplot(222, sharex=ax1, sharey=ax1)
        plt.title('Sawtooth')
        plt.plot(pix, pars2[:,0]%.5, '.', c='k')
        plt.ylim(0,.4)
        ax3 = plt.subplot(223, sharex=ax1)
        plt.plot(pix, np.abs(pars1[:,1]), '.', c='k')
        plt.ylabel('amplitude')
        plt.xlabel('y position')
        ax4 = plt.subplot(224, sharex=ax1, sharey=ax3)
        plt.plot(pix, np.abs(pars2[:,1]), '.', c='k')
        plt.xlabel('y position')
        plt.xlim(1630, 1670)
        plt.ylim(0.,.1)
    
    plt.show()
