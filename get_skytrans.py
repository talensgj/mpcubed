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
        
        
        # Create the indices.
        hg = HealpixGrid(64)
        bins, binnum = hg.find_gridpoint(ra, dec, compact=True)
        count = index_statistics(binnum, vmag, statistic='count', keeplength=False)
        plt.hist(count, bins = np.linspace(0.5, 50.5, 51), normed=True, histtype='step')
        
        hg = HealpixGrid(32)
        bins, binnum = hg.find_gridpoint(ra, dec, compact=True)
        count = index_statistics(binnum, vmag, statistic='count', keeplength=False)
        plt.hist(count, bins = np.linspace(0.5, 50.5, 51), normed=True, histtype='step')
        
        hg = HealpixGrid(16)
        bins, binnum = hg.find_gridpoint(ra, dec, compact=True)
        count = index_statistics(binnum, vmag, statistic='count', keeplength=False)
        plt.hist(count, bins = np.linspace(0.5, 50.5, 51), normed=True, histtype='step')
        
        plt.xlabel('# Stars')
        plt.show()
        exit()
        
        ind = np.where(ascc=='807144')
        here = binnum == binnum[ind]
        ascc = ascc[here]
        vmag = vmag[here]
        ra = ra[here]
        dec = dec[here]
        Nobs = Nobs[here]
        binnum = binnum[here]
        
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
            
        lstidx = lstidx.astype('int')
        times, lstidx = np.unique(lstidx, return_inverse=True)
        
        star_id = np.repeat(np.arange(len(ascc)), Nobs)
        
        # Remove bad data.
        here = np.where((flags < 1)&(cflux0>0)&(sky>0)&(eflux0>0))
        lstidx = lstidx[here]
        star_id = star_id[here]
        cflux0 = cflux0[here]
        eflux0 = eflux0[here]
        
        _, lstidx = np.unique(lstidx, return_inverse=True)
            
        a2, a1, niter, chisq, chisq_pbin2, chisq_pbin1, npoints, npars = sysrem(lstidx, star_id, cflux0, eflux0, a2=1e7*10**(vmag/-2.5))
        
        for i in np.unique(star_id):
            
            here = star_id==i
        
            ax = plt.subplot(311)
            plt.plot(lstidx[here], cflux0[here]/np.median(cflux0[here]), '.')
            
            plt.subplot(312, sharex=ax, sharey=ax)
            plt.plot(np.unique(lstidx), a2/np.median(a2), '.')
            
            plt.subplot(313, sharex=ax, sharey=ax)
            plt.plot(lstidx[here], (cflux0/a2[lstidx])[here]/np.median((cflux0/a2[lstidx])[here]), '.')
            plt.ylim(.8,1.2)
            
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
