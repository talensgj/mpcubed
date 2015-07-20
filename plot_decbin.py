#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt
from time import time

import os
import glob

from reduction_methods import sysrem
from scipy.optimize import minimize

def window_std(lst_idx, array):
    
    slow_id = lst_idx//50
    result = np.zeros(array.shape)
    for i in range(270):
        args, = np.where(slow_id==i)        
        result[args] = np.nanstd(array[args])
    
    args, = np.where(result == 0)
    result[args] = np.nan
    
    return result

def make_idx(variable, range, nbins):
    
    bins = np.linspace(range[0], range[1], nbins+1)
    idx = np.searchsorted(bins, variable)
    
    if np.any(idx == 0) | np.any(idx == len(bins)):
        print 'Warning: there where values out of range.'
        exit()
        
    idx -= 1
    idx = idx.astype('int')
    
    offset = np.amin(idx)
    length = np.ptp(idx) + 1
    
    return idx, length, offset

def fitsin(pars, x, y, yerr):
    
    fit = pars[0]*np.sin(2*np.pi*x+pars[1])+1
    
    return np.nansum(((y-fit)/yerr)**2)

with h5py.File('/data2/mascara/LaPalma/20150203LPE/fLC/fLC_20150203LPE.hdf5') as f:
    
    hdr = f['table_header']
    ascc = hdr['ascc'].value
    ra = hdr['ra'].value
    dec = hdr['dec'].value
    Nobs = hdr['nobs'].value.astype('int')
    
    ra_idx, length, offset = make_idx(ra, (0,360), 13500) # hardcoded...
    dec_idx, length, offset_dec = make_idx(dec, (-90,90), 1440) # hardcoded...
    
    args, = np.where((dec<28.25)&(dec>28))
    
    data = f['data']
    
    lst_idx = np.array([])
    eflux0 = np.array([])
    flux0 = np.array([])
    x = np.array([])
    flags = np.array([])
    for i in args:
        lst_idx = np.hstack([lst_idx, data[ascc[i]]['lstidx']])
        flux0 = np.hstack([flux0, data[ascc[i]]['flux0']])
        eflux0 = np.hstack([eflux0, window_std(data[ascc[i]]['lstidx'], data[ascc[i]]['flux0'])])
        #eflux0 = np.hstack([eflux0, data[ascc[i]]['eflux0']])
        flags = np.hstack([flags, data[ascc[i]]['flag']])
        x = np.hstack([x, data[ascc[i]]['x']])

eflux0[flags>0] = np.nan
flux0[np.isnan(eflux0)] = np.nan

        
lst_idx = lst_idx.astype('int')
ha_idx = lst_idx - np.repeat(ra_idx[args], Nobs[args])
ha_idx = np.mod(ha_idx, 13500) # hardcoded...
offset_ha = np.amin(ha_idx)

data = np.full((len(args), np.ptp(ha_idx)+1), fill_value=np.nan)
error = np.full((len(args), np.ptp(ha_idx)+1), fill_value=np.nan)
pos = np.full((len(args), np.ptp(ha_idx)+1), fill_value=np.nan)

select = np.append(0, np.cumsum(Nobs[args]))
for i in range(len(args)):
    data[i, ha_idx[select[i]:select[i+1]]-offset_ha] = flux0[select[i]:select[i+1]]
    error[i, ha_idx[select[i]:select[i+1]]-offset_ha] = eflux0[select[i]:select[i+1]]
    pos[i, ha_idx[select[i]:select[i+1]]-offset_ha] = x[select[i]:select[i+1]]

c_i, a_j, niter, chi2_new = sysrem(data, error)

print chi2_new/(np.sum(np.isfinite(data))-np.sum(np.isfinite(c_i))-np.sum(np.isfinite(a_j)))

ax = plt.subplot(311)
plt.imshow(data/c_i[:,np.newaxis], aspect='auto', interpolation='None', vmin=0.5, vmax=1.5)
plt.colorbar()

plt.subplot(312, sharex=ax)
plt.scatter(np.arange(len(a_j)), a_j, c=a_j, vmin=0.5, vmax=1.5)
plt.colorbar()
plt.ylim(0.5, 1.5)
plt.xlim(-0.5, len(a_j)-.5)

plt.subplot(313, sharex=ax)
plt.imshow(data/np.outer(c_i, a_j), aspect='auto', interpolation='None', vmin=0.8, vmax=1.2)
plt.colorbar()

plt.show()

exit()
for j in range(len(args)):
    print ascc[args[j]]
    lims, = np.where(np.isfinite(data[j]))
    
    if len(lims) == 0: continue
    
    ax = plt.subplot(411)
    plt.errorbar(np.arange(len(data[j])), data[j]/c_i[j], yerr=error[j]/c_i[j], fmt='.')
    plt.ylim(0.5,1.5)
    
    plt.subplot(412, sharex=ax)
    plt.plot(a_j, '.')
    plt.ylim(0.5,1.5)
    
    plt.subplot(413, sharex=ax)
    plt.errorbar(np.arange(len(data[j])), data[j]/(a_j*c_i[j]), yerr=error[j]/(a_j*c_i[j]), fmt='.')
    plt.xlim(np.amin(lims), np.amax(lims))
    plt.ylim(0.9,1.1)
    
    plt.subplot(414)
    plt.errorbar(pos[j]%1, data[j]/(a_j*c_i[j]), yerr=error[j]/(a_j*c_i[j]), fmt='.')
    plt.ylim(0.9,1.1)
    
    plt.show()
    plt.close()


    
    
    
