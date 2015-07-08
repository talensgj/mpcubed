#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np
from reduction_functions import make_array, make_transmission_map

import matplotlib.pyplot as plt
from time import time

import os
import glob

from scipy.stats import binned_statistic_2d

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
    
def window_std(lst_idx, array):
    
    if np.isscalar(array):
        return np.nan
    
    slow_id = lst_idx//50
    result = np.zeros(array.shape)
    for i in range(270):
        args, = np.where(slow_id==i)        
        result[args] = np.std(array[args])
    
    args, = np.where(result == 0)
    result[args] = np.nan
    
    return result
    
def transmission(filename):

    with h5py.File(filename) as f:
        
        hdr = f['table_header']
        ascc = hdr['ascc'].value
        vmag = hdr['vmag'].value
        ra = hdr['ra'].value
        dec = hdr['dec'].value
        Nobs = hdr['nobs'].value.astype('int')
        
        select = np.append(0, np.cumsum(Nobs))
    
        lst_idx = np.zeros(np.sum(Nobs)).astype('int')
        flux0 = np.zeros(np.sum(Nobs))
        eflux0 = np.zeros(np.sum(Nobs))
        x = np.zeros(np.sum(Nobs))
        y = np.zeros(np.sum(Nobs))
        flags = np.zeros(np.sum(Nobs))
        
        data = f['data']
        for i in range(len(ascc)):
            lst_idx[select[i]:select[i+1]] = data[ascc[i]]['lstidx']
            flux0[select[i]:select[i+1]] = data[ascc[i]]['flux0']
            #eflux0[select[i]:select[i+1]] = data[ascc[i]]['eflux0']
            eflux0[select[i]:select[i+1]] = window_std(data[ascc[i]]['lstidx'], data[ascc[i]]['flux0'])
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']
            x[select[i]:select[i+1]] = data[ascc[i]]['x']
            y[select[i]:select[i+1]] = data[ascc[i]]['y']

    #Create a position index.
    stat, xedg, yedg, binnum = binned_statistic_2d(x, y, flux0, statistic='count', bins=(np.linspace(0,4008,101),np.linspace(0,2672,76)))
    stat0, xedg, yedg, binnum0 = binned_statistic_2d(x, y, flux0, statistic='count', bins=(np.linspace(0,4008,201),np.linspace(0,2672,151)))
    stat1, xedg, yedg, binnum1 = binned_statistic_2d(x, y, flux0, statistic='count', bins=(np.linspace(0,4008,401),np.linspace(0,2672,301)))
    stat2, xedg, yedg, binnum2 = binned_statistic_2d(x, y, flux0, statistic='count', bins=(np.linspace(0,4008,801),np.linspace(0,2672,601)))
    
    #Create a star index.
    star_id = np.repeat(np.arange(len(ascc)), Nobs)
    
    # Set bad data to NaN.
    here, = np.where((flags < 1)&np.isfinite(eflux0))
    flux0 = flux0[here]
    eflux0 = eflux0[here]
    binnum = binnum[here]
    binnum0 = binnum0[here]
    binnum1 = binnum1[here]
    binnum2 = binnum2[here]
    star_id = star_id[here]

    a_j = np.ones(802*602)
    for i in range(20):
        c_i = np.bincount(star_id, flux0*a_j[binnum2]/eflux0**2)/np.bincount(star_id, (a_j**2)[binnum2]/eflux0**2)
        a_j = np.bincount(binnum, flux0*c_i[star_id]/eflux0**2)/np.bincount(binnum, (c_i**2)[star_id]/eflux0**2)
        c_i = np.bincount(star_id, flux0*a_j[binnum]/eflux0**2)/np.bincount(star_id, (a_j**2)[binnum]/eflux0**2)
        a_j = np.bincount(binnum0, flux0*c_i[star_id]/eflux0**2)/np.bincount(binnum0, (c_i**2)[star_id]/eflux0**2)
        c_i = np.bincount(star_id, flux0*a_j[binnum0]/eflux0**2)/np.bincount(star_id, (a_j**2)[binnum0]/eflux0**2)
        a_j = np.bincount(binnum1, flux0*c_i[star_id]/eflux0**2)/np.bincount(binnum1, (c_i**2)[star_id]/eflux0**2)
        c_i = np.bincount(star_id, flux0*a_j[binnum1]/eflux0**2)/np.bincount(star_id, (a_j**2)[binnum1]/eflux0**2)
        a_j = np.bincount(binnum2, flux0*c_i[star_id]/eflux0**2)/np.bincount(binnum2, (c_i**2)[star_id]/eflux0**2)
        
    res = np.full((802*602), fill_value=np.nan)
    res[np.arange(len(a_j))] = a_j
    res = res.reshape((602,802))
    
    plt.semilogy(vmag[np.unique(star_id)], c_i[np.unique(star_id)], '.')
    plt.show()

    #stat[stat<1] = np.nan
    #plt.imshow(stat.T, interpolation='None', vmax = np.nanmax(stat), vmin = np.nanmin(stat[stat>0]))
    #plt.colorbar()
    #plt.show()
    
    plt.imshow(res, interpolation='None', vmin=0, vmax=2)
    plt.show()

    return 0

if __name__ == '__main__':
    #import argparse
    
    #parser = argparse.ArgumentParser(description='Compute the transmission map for a particular night and camera.')
    
    #parser.add_argument('path', type=str, help='the global path for search')
    #parser.add_argument('-n', '--night', default='', type=str, help='the night to reduce')
    #parser.add_argument('-c', '--camera',default='', type=str, help='the camera to reduce')
    
    #args = parser.parse_args()
    
    filelist = glob.glob('/data2/mascara/LaPalma/20150203LPW/fLC/fLC_*.hdf5')
    filelist = np.sort(filelist)

    filelist = ['/data2/talens/fLC_20150203LPC.hdf5']

    for filename in filelist:
        print 'Data:', filename
        transmission(filename)
