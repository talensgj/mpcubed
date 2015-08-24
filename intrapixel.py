#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt

from index_functions import index_statistics

fLC = '/data2/talens/Jul2015/fLC_20150716LPC.hdf5'
red = '/data2/talens/Jul2015/red_20150716LPC.hdf5'

from scipy.optimize import minimize

def intrapixel(pars, x, y):
    
    fit = pars[0]*(pars[2]*np.cos(2*np.pi*x+pars[1])+1)
    
    return np.sum(np.abs(y-fit))

def intrapixel_eval(pars, x):
    fit = pars[0]*(pars[2]*np.cos(2*np.pi*x+pars[1])+1)
    return fit
    
    
def intrapixel1(pars, x, y):
    
    fit = pars[0]*np.sin(2*np.pi*x+pars[1])+1
    
    return np.sum(np.abs(y-fit))
    
def intrapixel1_eval(pars, x):
    
    fit = pars[0]*np.sin(2*np.pi*x+pars[1])+1
    
    return fit
    
    
with h5py.File(fLC, 'r') as f, h5py.File(red, 'r') as g:
    
    lc = f['data/807144']
    rc = g['data/807144']
    
    jd = lc['jdmid']
    x = lc['y']
    y = rc['cflux0']
    
    mean = index_statistics(lc['lstidx']//50, lc['peak']/lc['flux0'], statistic='mean')
    
    slow_idx = lc['lstidx']//50
    ind, slow_idx = np.unique(slow_idx, return_inverse=True)
    
    pars = np.zeros((len(ind),3))
    fit = np.ones(len(y))
    for i in range(len(ind)):
        here = slow_idx == i

        res = minimize(intrapixel, [np.mean(y[here]), np.pi/2, .1], args=(x[here], y[here]))
        fit[here] = intrapixel_eval(res.x, x[here])
        pars[i] = res.x

    
    pars1 = np.zeros((len(ind), 2))
    fit1 = np.ones(len(y))
    niter = 0
    while niter < 5:
    
        F = np.mean(y/fit1)
        print F
        for i in range(len(ind)):
            here = slow_idx == i
    
            res = minimize(intrapixel1, [.1, np.pi/2], args=(x[here], y[here]/F))
            fit1[here] = intrapixel1_eval(res.x, x[here])
            pars1[i] = res.x
            
        niter += 1
    
    fit1 = F*fit1
    fit2 = fit/np.repeat(pars[:,0], np.bincount(slow_idx))
        
    ax = plt.subplot(311)
    plt.plot(jd, y, '.', c='k')
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(jd, fit, '.', c='b')
    plt.plot(jd, fit1, '.', c='r')
    plt.plot(jd, fit2*(np.median(y/fit2)), c='c')
    plt.subplot(313, sharex=ax)
    plt.plot(jd, y/fit, '.', c='b')
    plt.plot(jd, y/fit1, '.', c='r')
    plt.plot(jd, y/fit2/(np.median(y/fit2)), '.', c='c')
    plt.show()
    
    plt.subplot(311)
    plt.plot(pars[:,0], '.', c='b')
    plt.axhline(F, c='r')
    plt.subplot(312)
    plt.plot(np.abs(pars[:,2]), '.', c='b')
    plt.plot(np.abs(pars1[:,0]), '.', c='r')
    plt.subplot(313)
    plt.plot(pars[:,1]%np.pi, '.', c='b')
    plt.plot(pars1[:,1]%np.pi, '.', c='r')
    plt.show()
    
    plt.plot(np.abs(pars[:,2]), mean, '.', c='b')
    plt.plot(np.abs(pars1[:,0]), mean, '.', c='r')
    plt.show()
