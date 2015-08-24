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

    with h5py.File(filename) as f:
        
        hdr = f['table_header']
        ascc = hdr['ascc'].value
        vmag = hdr['vmag'].value
        ra = hdr['ra'].value
        dec = hdr['dec'].value
        Nobs = hdr['nobs'].value.astype('int')
        
        select = np.append(0, np.cumsum(Nobs))
    
        lst = np.zeros(np.sum(Nobs))
        flux0 = np.zeros(np.sum(Nobs))
        eflux0 = np.zeros(np.sum(Nobs))
        sky = np.zeros(np.sum(Nobs))
        flags = np.zeros(np.sum(Nobs))
        x = np.zeros(np.sum(Nobs))
        y = np.zeros(np.sum(Nobs))
        
        data = f['data']
        for i in range(len(ascc)):
            lst[select[i]:select[i+1]] = data[ascc[i]]['lst']
            flux0[select[i]:select[i+1]] = data[ascc[i]]['flux0']
            #eflux0[select[i]:select[i+1]] = data[ascc[i]]['eflux0']
            eflux0[select[i]:select[i+1]] = index_statistics(data[ascc[i]]['lstidx']//50, data[ascc[i]]['flux0'], statistic='std', keeplength=True)
            sky[select[i]:select[i+1]] = data[ascc[i]]['sky']
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']
            x[select[i]:select[i+1]] = data[ascc[i]]['x']
            y[select[i]:select[i+1]] = data[ascc[i]]['y']
    
    dec_1 = np.copy(dec)
    
    print 'Done reading.'

    # Sky coordinates of data.
    ha = np.mod(lst*15.-np.repeat(ra,Nobs), 360.)
    dec = np.repeat(dec, Nobs)

    # Create the indices.
    star_id = np.repeat(np.arange(len(ascc)), Nobs)
    
    # Remove bad data.
    here, = np.where((flags < 1)&(flux0>0)&(sky>0)&(eflux0>0))
    flux0 = flux0[here]
    eflux0 = eflux0[here]
    ha = ha[here]
    dec = dec[here]
    x = x[here]
    y = y[here]
    star_id = star_id[here]

    hg = PolarGrid(2700, 720)
    bins, binnum = hg.find_gridpoint(ha, dec, compact=True)
    
    #hg = CartesianGrid(400, 300, margin=50)
    #bins, binnum = hg.find_gridpoint(x, y, compact=True)
    
    #hg = HealpixGrid(256)
    #bins, binnum = hg.find_gridpoint(ha, dec, compact=True)
    
    count = index_statistics(binnum, flux0, statistic='count', keeplength=False)
    print np.percentile(count, 5), np.percentile(count, 95)
    
    # Compute the transmission using sysrem.
    a2, a1, niter, chisq, chisq_pbin2, chisq_pbin1, npoints, npars = sysrem(binnum, star_id, flux0, eflux0, a2=1e7*10**(vmag/-2.5))
    
    with h5py.File('/data2/talens/Jul2015/cam_20150716LPC_pg2700x720.hdf5') as f:
        
        grp = f.create_group('header')
        grp.attrs['grid'] = 'polar'
        grp.attrs['nx'] = 2700
        grp.attrs['ny'] = 720
        grp.attrs['margin'] = 0
        grp.attrs['niter'] = niter
        grp.attrs['npoints'] = npoints
        grp.attrs['npars'] = npars
        grp.attrs['chisq'] = chisq
        
        grp = f.create_group('data')
        grp.create_dataset('binnum', data = bins)
        grp.create_dataset('count', data = count)
        grp.create_dataset('trans', data = a2)
        grp.create_dataset('chisq_trans', data = chisq_pbin2)
        
        stars = np.unique(star_id)
        grp.create_dataset('ascc', data = ascc[stars])
        grp.create_dataset('vmag', data = vmag[stars])
        grp.create_dataset('flux', data = a1[stars])
        grp.create_dataset('dec', data = dec_1[stars])
        grp.create_dataset('chisq_flux', data = chisq_pbin1[stars])
        
    return 0

if __name__ == '__main__':
    #import argparse
    
    #parser = argparse.ArgumentParser(description='Compute the transmission map for a particular night and camera.')
    
    #parser.add_argument('path', type=str, help='the global path for search')
    #parser.add_argument('-n', '--night', default='', type=str, help='the night to reduce')
    #parser.add_argument('-c', '--camera',default='', type=str, help='the camera to reduce')
    
    #args = parser.parse_args()
    
    filelist = glob.glob('/data2/talens/Jul2015/fLC_20150716LPC.hdf5')
    #filelist = np.sort(filelist)

    #filelist = ['/data2/talens/fLC_20150203LPS.hdf5']

    for filename in filelist:
        print 'Data:', filename
        transmission(filename)
