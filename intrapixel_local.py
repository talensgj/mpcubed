#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib import rcParams

from index_functions import index_statistics

fLC = '/data2/talens/Jul2015/fLC_20150716LPC.hdf5'
red = '/data2/talens/Jul2015/red_20150716LPC.hdf5'

from scipy.optimize import minimize
from scipy.ndimage.filters import gaussian_filter1d

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'

def intrapixel1(pars, x, y=None, yerr=None):
    
    fit = pars[0]*(pars[1]*np.cos(2*np.pi*(x+pars[2]))+1)
    
    if y is None:
        return fit
    elif yerr is None:
        return np.sum(np.abs((y-fit)))
    else:
        return np.sum((y-fit)**2/yerr**2)

def intrapixel2(pars, x, y=None, yerr=None):
    
    fit = (pars[0]*np.cos(2*np.pi*x+pars[1])+1)
    
    if y is None:
        return fit
    elif yerr is None:
        return np.sum(np.abs((y-fit)))
    else:
        return np.sum((y-fit)**2/yerr**2)
        
def intrapixel3(pars, x1, x2, y=None, yerr=None):
    
    fit = (np.polyval(pars[:2], x2)*np.cos(2*np.pi*x1+pars[2])+1)
    
    if y is None:
        return fit
    elif yerr is None:
        return np.sum(np.abs((y-fit)))
    else:
        return np.sum((y-fit)**2/yerr**2)
        
def intrapixel4(pars, x1, x2, count, y=None, yerr=None):
    
    fit = (np.polyval(pars[:2], x2)*np.cos(2*np.pi*x1+np.repeat(pars[2:], count))+1)
    
    if y is None:
        return fit
    elif yerr is None:
        return np.sum(np.abs((y-fit)))
    else:
        return np.sum((y-fit)**2/yerr**2)

def sawtooth(pars, x, y=None, yerr=None):
    
    xtmp = x-pars[0]
    
    fit = pars[2]*(pars[1]*(np.abs(xtmp-np.floor(xtmp)-0.5)*4-1)+1)
    
    if y is None:
        return fit
    elif yerr is None:
        return np.sum(np.abs((y-fit)))
    else:
        return np.sum((y-fit)**2/yerr**2)
        

with h5py.File(fLC, 'r') as f, h5py.File(red, 'r') as g:
    
    lc = f['data/807144']
    rc = g['data/807144']
    
    jd = lc['jdmid']-np.floor(lc['jdmid'])
    x1 = lc['y']
    x2 = lc['peak']/lc['flux0']
    gauss = gaussian_filter1d(rc['cflux0'], 100)
    y = rc['cflux0']/gauss
    
    x3 = x1-np.floor(x1)
    
    ax = plt.subplot(311)
    plt.plot(jd, y, '.')
    plt.subplot(312, sharex=ax)
    plt.plot(jd, np.abs(x3-.5), '.')
    plt.subplot(313, sharex=ax)
    plt.plot(jd, np.cos(2*np.pi*x1), '.')
    plt.xlim(.53,.57)
    plt.tight_layout()
    plt.show()
    
    
    pix = np.unique(np.floor(x1))
    fit = np.ones(len(y))
    fit1 = np.ones(len(y))
    pars = np.full((len(pix), 3), fill_value=np.nan)
    pars1 = np.full((len(pix), 3), fill_value=np.nan)
    for i in range(len(pix)):
    
        here = np.floor(x1) == pix[i]
        if np.sum(here)<25: continue
        
        res = minimize(sawtooth, [0., .1, 1.], args=(x1[here], y[here]))
        pars[i] = res.x
        fit[here] = sawtooth(res.x, x1[here])
   
        res = minimize(intrapixel1, [1., .1, 0.], args=(x1[here], y[here]))
        pars1[i] = res.x
        fit1[here] = intrapixel1(res.x, x1[here])
   
    plt.subplot(311)
    plt.plot(pix, pars[:,0]%.5, '.')
    plt.plot(pix, np.abs(pars1[:,2]%.5-.5), '.')
    plt.subplot(312)
    plt.plot(pix, np.abs(pars[:,1]), '.')
    plt.plot(pix, np.abs(pars1[:,1]), '.')
    plt.subplot(313)
    plt.plot(pix, pars[:,2], '.')
    plt.plot(pix, pars1[:,0], '.')
    plt.show()
   
    ax = plt.subplot(211)
    plt.plot(jd, y, '.', c='k')
    plt.plot(jd, fit, '.', c='r')
    plt.plot(jd, fit1, '.', c='c')
    
    plt.subplot(212, sharex=ax, sharey=ax)
    plt.plot(jd, y/fit, '.', c='r')
    plt.plot(jd, y/fit1, '.', c='c')
    plt.xlim(.53,.57)
    plt.show()
    exit()
    plt.figure(figsize=(16,8))
    ax = plt.subplot(311)
    plt.plot(jd, rc['cflux0'], '.', c='k')
    plt.plot(jd, gauss, c='r')
    plt.ylabel('flux0')
    plt.subplot(312, sharex=ax)
    plt.plot(jd, np.abs(x1-np.floor(x1)-.5), '.', c='k')
    plt.ylabel(r'$cos(2\pi y)$')
    plt.subplot(313, sharex=ax)
    plt.plot(jd, x2, '.', c='k')
    plt.xlabel('Time [JD-2457220]')
    plt.ylabel('peak/flux0')
    plt.xlim(np.amin(jd), np.amax(jd))
    plt.tight_layout()
    plt.show()

    res = minimize(intrapixel1, [1., 0., 0.], args=(x1, y))
    fit = intrapixel1(res.x, x1)
    
    plt.figure(figsize=(16,8))
    ax = plt.subplot(311)
    plt.plot(jd, y, '.', c='k')
    plt.plot(jd, fit, '.', c='r')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data, Fit')
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(jd, y/fit, '.', c='k')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data/Fit')
    plt.subplot(313, sharex=ax)
    plt.plot(jd, y-fit, '.', c='k')
    plt.xlim(np.amin(jd), np.amax(jd))
    plt.xlabel('Time [JD-2457220]')
    plt.ylabel('Data-Fit')
    plt.tight_layout()
    plt.show()
    
    slow_ind = lc['lstidx']//50
    fit = np.zeros(len(y))
    ind, slow_ind = np.unique(slow_ind, return_inverse=True)
    pars = np.zeros((len(ind), 3))
    for i in range(len(ind)):
        here = slow_ind == i
        res = minimize(intrapixel1, [1., 0., 0.], args=(x1[here], y[here]))
        fit[here] = intrapixel1(res.x, x1[here])
        pars[i] = res.x
    
    plt.figure(figsize=(16,8))
    ax = plt.subplot(311)
    plt.plot(jd, y, '.', c='k')
    plt.plot(jd, fit, '.', c='r')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data, Fit')
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(jd, y/fit, '.', c='k')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data/Fit')
    plt.subplot(313, sharex=ax)
    plt.plot(jd, y-fit, '.', c='k')
    plt.xlim(np.amin(jd), np.amax(jd))
    plt.xlabel('Time [JD-2457220]')
    plt.ylabel('Data-Fit')
    plt.tight_layout()
    plt.show()
    
    slow_ind = lc['lstidx']//50
    fit = np.zeros(len(y))
    ind, slow_ind = np.unique(slow_ind, return_inverse=True)
    pars = np.zeros((len(ind), 2))
    for i in range(len(ind)):
        here = slow_ind == i
        res = minimize(intrapixel2, [0., 0.], args=(x1[here], y[here]))
        fit[here] = intrapixel2(res.x, x1[here])
        pars[i] = res.x
    
    plt.figure(figsize=(16,8))
    ax = plt.subplot(311)
    plt.plot(jd, y, '.', c='k')
    plt.plot(jd, fit, '.', c='r')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data, Fit')
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(jd, y/fit, '.', c='k')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data/Fit')
    plt.subplot(313, sharex=ax)
    plt.plot(jd, y-fit, '.', c='k')
    plt.xlim(np.amin(jd), np.amax(jd))
    plt.xlabel('Time [JD-2457220]')
    plt.ylabel('Data-Fit')
    plt.tight_layout()
    plt.show()
    
    res = minimize(intrapixel3, [0., 0., 0.], args=(x1, x2, y))
    fit = intrapixel3(res.x, x1, x2)
    
    plt.figure(figsize=(16,8))
    ax = plt.subplot(311)
    plt.plot(jd, y, '.', c='k')
    plt.plot(jd, fit, '.', c='r')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data, Fit')
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(jd, y/fit, '.', c='k')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data/Fit')
    plt.subplot(313, sharex=ax)
    plt.plot(jd, y-fit, '.', c='k')
    plt.xlim(np.amin(jd), np.amax(jd))
    plt.xlabel('Time [JD-2457220]')
    plt.ylabel('Data-Fit')
    plt.tight_layout()
    plt.show()
    
    slow_ind = lc['lstidx']//50
    ind, slow_ind = np.unique(slow_ind, return_inverse=True)
    count = np.bincount(slow_ind)
    res = minimize(intrapixel4, [0., 0.]+[0.]*len(ind), args=(x1, x2, count, y))
    fit = intrapixel4(res.x, x1, x2, count)
    
    plt.figure(figsize=(16,8))
    ax = plt.subplot(311)
    plt.plot(jd, y, '.', c='k')
    plt.plot(jd, fit, '.', c='r')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data, Fit')
    plt.subplot(312, sharex=ax, sharey=ax)
    plt.plot(jd, y/fit, '.', c='k')
    plt.axhline(np.mean(y)+np.std(y))
    plt.axhline(np.mean(y)-np.std(y))
    plt.ylabel('Data/Fit')
    plt.subplot(313, sharex=ax)
    plt.plot(jd, y-fit, '.', c='k')
    plt.xlim(np.amin(jd), np.amax(jd))
    plt.xlabel('Time [JD-2457220]')
    plt.ylabel('Data-Fit')
    plt.tight_layout()
    plt.show()
    
