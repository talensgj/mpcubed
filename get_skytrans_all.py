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

from time import time
    
def skytrans(fLC, red):

    ascc = np.array([])
    vmag = np.array([])
    ra = np.array([])
    dec = np.array([])
    Nobs = np.array([]).astype('int')
    cardinal = np.array([])

    for camera in ['LPN', 'LPE', 'LPS', 'LPW', 'LPC']:

        with h5py.File('/data2/talens/Jul2015/fLC_20150714%s.hdf5'%camera, 'r') as f:
            
            hdr = f['table_header']
            ascc = np.append(ascc, hdr['ascc'].value)
            vmag = np.append(vmag, hdr['vmag'].value)
            ra = np.append(ra, hdr['ra'].value)
            dec = np.append(dec, hdr['dec'].value)
            Nobs = np.append(Nobs, hdr['nobs'].value.astype('int'))
            cardinal = np.append(cardinal, [camera]*len(hdr['ascc'].value))
    
    grid = HealpixGrid(8)
    bins, binnum = grid.find_gridpoint(ra, dec, compact=True)
    starcount = index_statistics(binnum, ra, statistic='count')
    #hmap = grid.put_values_on_grid(count, bins)
    
    #healpy.mollview(hmap)
    #plt.show()
    print binnum[np.where(ascc=='807144')]
    for ind in range(len(bins)):
        
        here = binnum == ind
        ascc_tmp = ascc[here]
        vmag_tmp = vmag[here]
        Nobs_tmp = Nobs[here]
        cardinal_tmp = cardinal[here]
        ra_tmp = ra[here]
        
        select = np.append(0, np.cumsum(Nobs_tmp))

        lstidx = np.zeros(np.sum(Nobs_tmp))
        sky = np.zeros(np.sum(Nobs_tmp))
        flags1 = np.zeros(np.sum(Nobs_tmp))
        
        cflux0 = np.zeros(np.sum(Nobs_tmp))
        ecflux0 = np.zeros(np.sum(Nobs_tmp))
        flags2 = np.zeros(np.sum(Nobs_tmp))
        
        for i in range(len(ascc_tmp)):
        
            with h5py.File('/data2/talens/Jul2015/fLC_20150714%s.hdf5'%cardinal_tmp[i], 'r') as f, h5py.File('/data2/talens/Jul2015/red_20150714%s.hdf5'%cardinal_tmp[i], 'r') as g:
            
                lc = f['data']
                rc = g['data']
            
                lstidx[select[i]:select[i+1]] = lc[ascc_tmp[i]]['lstidx']
                sky[select[i]:select[i+1]] = lc[ascc_tmp[i]]['sky']
                flags1[select[i]:select[i+1]] = lc[ascc_tmp[i]]['flag']
    
                cflux0[select[i]:select[i+1]] = rc[ascc_tmp[i]]['cflux0']
                ecflux0[select[i]:select[i+1]] = rc[ascc_tmp[i]]['ecflux0']
                flags2[select[i]:select[i+1]] = rc[ascc_tmp[i]]['flags']
           
        print 'Done reading.'
        
        #lstidx = lstidx.astype('int')
        #times, lstidx = np.unique(lstidx, return_inverse=True)
        #star_id = np.repeat(np.arange(len(ascc_tmp)), Nobs_tmp)

        #array = np.full((len(times), len(ascc_tmp)), fill_value=np.nan)
        #array[lstidx, star_id] = cflux0
        #array = array[:, np.argsort(ra_tmp)]
        
        #plt.imshow((array/np.nanmedian(array, axis=0)).T, interpolation='none', aspect='auto', cmap=viridis, vmin=.5, vmax=1.5)
        #plt.colorbar()
        #plt.show()

        lstidx = lstidx.astype('int')
        count = np.bincount(lstidx)
        
        star_id = np.repeat(np.arange(len(ascc_tmp)), Nobs_tmp)
        
        here = (flags1<1)&(flags2<1)&(cflux0>0)&(sky>0)&(ecflux0>0)
        cflux0 = cflux0[here]
        ecflux0 = ecflux0[here]
        lstidx = lstidx[here]
        star_id = star_id[here]
        
        if len(cflux0) == 0:
            continue
            
        times, lstidx = np.unique(lstidx, return_inverse=True)
        a2, a1, niter, chisq, chisq_pbin2, chisq_pbin1, npoints, npars = sysrem(lstidx, star_id, cflux0, ecflux0, a2=1e7*10**(vmag_tmp/-2.5))
        
        sky = np.full(13500, fill_value = np.nan)
        sky[times[np.arange(len(a2))]] = a2
    
        for i in range(len(ascc_tmp)):
        
            with h5py.File('/data2/talens/Jul2015/fLC_20150714%s.hdf5'%cardinal_tmp[i], 'r') as f, h5py.File('/data2/talens/Jul2015/red_20150714%s.hdf5'%cardinal_tmp[i], 'r+') as g:
                print cardinal_tmp[i], ascc_tmp[i]
                lc = f['data']
                rc = g['data']
                
                try:
                    grp = g['data3']
                except:
                    grp = g.create_group('data3')
                else:
                    if ind == 0 & i==0:
                        del g['data3']
                        grp = g.create_group('data3')
               
                sky0 = sky[lc[ascc_tmp[i]]['lstidx'].astype('int')]
                scflux0 = rc[ascc_tmp[i]]['cflux0']/sky0
                
                flags = np.where(np.isnan(scflux0), 1, 0)
                flags = flags + np.where(count[lc[ascc_tmp[i]]['lstidx'].astype('int')]<=5, 2, 0)
                if starcount[ind] <= 5:
                    flags = flags + 4
                    
                record = np.rec.fromarrays([sky0, scflux0, flags], names=['sky0', 'scflux0', 'flags'])
                grp.create_dataset(ascc_tmp[i], data=record)
                
    return 0

if __name__ == '__main__':
    skytrans('/data2/talens/Jul2015/fLC_20150714LPW.hdf5', '/data2/talens/Jul2015/red_20150714LPW.hdf5')
