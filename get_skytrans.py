#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob

import h5py
import numpy as np

import healpy

import matplotlib.pyplot as plt
from viridis import viridis

from index_functions import index_statistics
from coordinate_grids import HealpixGrid, PolarGrid, CartesianGrid
from sysrem import sysrem

def transmission(filename):

    with h5py.File(filename) as f, h5py.File('/data2/talens/Jul2015/red_20150714LPC.hdf5') as g:
        
        hdr = f['table_header']
        ascc = hdr['ascc'].value
        vmag = hdr['vmag'].value
        ra = hdr['ra'].value
        dec = hdr['dec'].value
        Nobs = hdr['nobs'].value.astype('int')
        
        select = np.append(0, np.cumsum(Nobs))
    
        lstidx = np.zeros(np.sum(Nobs))
        cflux0 = np.zeros(np.sum(Nobs))
        eflux0 = np.zeros(np.sum(Nobs))
        sky = np.zeros(np.sum(Nobs))
        flags = np.zeros(np.sum(Nobs))
        
        data = f['data']
        red = g['data']
        for i in range(len(ascc)):
            lstidx[select[i]:select[i+1]] = data[ascc[i]]['lstidx']
            cflux0[select[i]:select[i+1]] = red[ascc[i]+'/cflux0'].value
            eflux0[select[i]:select[i+1]] = index_statistics(data[ascc[i]]['lstidx']//50, red[ascc[i]+'/cflux0'].value, statistic='std', keeplength=True)
            sky[select[i]:select[i+1]] = data[ascc[i]]['sky']
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']
           
        print 'Done reading.'
        
        # Create the indices.
        hg = HealpixGrid(32)
        bins, binnum = hg.find_gridpoint(ra, dec, compact=True)
            
        lstidx = lstidx.astype('int')
        times, lstidx = np.unique(lstidx, return_inverse=True)
        print len(times)
        sky_id = np.ravel_multi_index([np.repeat(binnum,Nobs),lstidx], (len(bins), len(times)))
        
        star_id = np.repeat(np.arange(len(ascc)), Nobs)
        
        # Remove bad data.
        here = np.where((flags < 1)&(cflux0>0)&(sky>0)&(eflux0>0))
        sky_id = sky_id[here]
        star_id = star_id[here]
        cflux0 = cflux0[here]
        eflux0 = eflux0[here]
            
        a2, a1, niter, chisq, chisq_pbin2, chisq_pbin1, npoints, npars = sysrem(sky_id, star_id, cflux0, np.ones(len(cflux0)), a2=1e7*10**(vmag/-2.5))
        
        array = np.full((len(bins), len(times)), fill_value=np.nan)
        array[np.unravel_index(sky_id, (len(bins), len(times)))] = a2[sky_id]
        
        plt.imshow(array, interpolation='None', origin='lower', aspect='auto', vmin=0, vmax=1.5, cmap=viridis)
        plt.show()
        
        return 0

if __name__ == '__main__':
    #import argparse
    
    #parser = argparse.ArgumentParser(description='Compute the transmission map for a particular night and camera.')
    
    #parser.add_argument('path', type=str, help='the global path for search')
    #parser.add_argument('-n', '--night', default='', type=str, help='the night to reduce')
    #parser.add_argument('-c', '--camera',default='', type=str, help='the camera to reduce')
    
    #args = parser.parse_args()
    
    filelist = glob.glob('/data2/talens/Jul2015/fLC_20150714LPC.hdf5')
    #filelist = np.sort(filelist)

    #filelist = ['/data2/talens/fLC_20150203LPC.hdf5']

    for filename in filelist:
        print 'Data:', filename
        transmission(filename)
