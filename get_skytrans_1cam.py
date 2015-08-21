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

def smooth_func(x, x_eval, y, sigma):
    delta_x = x_eval[:,None] - x
    weights = np.exp(-delta_x*delta_x / (2*sigma*sigma))
    weights /= np.sum(weights, axis=1, keepdims=True)
    y_eval = np.dot(weights, y)
    return y_eval

def smooth(ind1, ind2, par1, values, errors, window, a2=None, maxiter=50, eps=1e-3):
    
    if a2 is None:
        a2 = np.ones(np.amax(ind2)+1)
    
    s = values/errors**2
    r = 1./errors**2
    
    niter = 0
    end_crit = True
    
    #delta_x = par1[:,None] - par1[ind1]
    #weights = np.exp(-delta_x*delta_x / (2*window**2))
    #weights = weights/errors**2
    #weights /= np.sum(weights, axis=1, keepdims=True)
    
    while (niter < maxiter) & (end_crit):
        
        a1 = np.zeros(len(par1))
        for i in range(len(par1)):
            weights = np.exp(-(par1-par1[i])**2/(2.*window**2))
            a1[i] = np.sum(s*a2[ind2]*weights[ind1])/np.sum(r*(a2**2)[ind2]*weights[ind1])

        #temp = weights*(a2**2)[ind2]
        #temp /= np.sum(temp, axis=1, keepdims=True)
        #a1 = np.dot(temp, values/a2[ind2])

        a2 = np.bincount(ind2, s*a1[ind1])/np.bincount(ind2, r*(a1**2)[ind1])
        
        sol = a1[ind1]*a2[ind2]
        
        if niter == 0:
            end_crit = True
        else:
            crit = np.nanmax(np.abs((solo-sol)/solo))
            end_crit = (crit > eps)
        
        solo = np.copy(sol)
        
        niter += 1
        
    return a1, a2
    

def skytrans(fLC, red):

    with h5py.File(fLC, 'r') as f:
        
        hdr = f['table_header']
        ascc = hdr['ascc'].value
        vmag = hdr['vmag'].value
        ra = hdr['ra'].value
        dec = hdr['dec'].value
        Nobs = hdr['nobs'].value.astype('int')
        
    grid = HealpixGrid(8)
    bins, binnum = grid.find_gridpoint(ra, dec, compact=True)
    starcount = np.bincount(binnum)
    
    print 'Processing %i positions on the sky.'%len(bins)
    
    for ind in range(len(bins)):
        
        here = binnum == ind
        ascc_tmp = ascc[here]
        vmag_tmp = vmag[here]
        Nobs_tmp = Nobs[here]
        
        select = np.append(0, np.cumsum(Nobs_tmp))

        lstidx = np.zeros(np.sum(Nobs_tmp))
        sky = np.zeros(np.sum(Nobs_tmp))
        flags1 = np.zeros(np.sum(Nobs_tmp))
        
        cflux0 = np.zeros(np.sum(Nobs_tmp))
        ecflux0 = np.zeros(np.sum(Nobs_tmp))
        flags2 = np.zeros(np.sum(Nobs_tmp))
            
        with h5py.File(fLC, 'r') as f, h5py.File(red, 'r') as g:
            
            lc = f['data']
            rc = g['data']
            
            for i in range(len(ascc_tmp)):
                lstidx[select[i]:select[i+1]] = lc[ascc_tmp[i]]['lstidx']
                sky[select[i]:select[i+1]] = lc[ascc_tmp[i]]['sky']
                flags1[select[i]:select[i+1]] = lc[ascc_tmp[i]]['flag']
                
                cflux0[select[i]:select[i+1]] = rc[ascc_tmp[i]]['cflux0']
                ecflux0[select[i]:select[i+1]] = rc[ascc_tmp[i]]['ecflux0']
                flags2[select[i]:select[i+1]] = rc[ascc_tmp[i]]['flags']
           
        print 'Done reading.'
        
        lstidx = lstidx.astype('int')
        count = np.bincount(lstidx)
        
        star_id = np.repeat(np.arange(len(ascc_tmp)), Nobs_tmp)
        
        here = (flags1<1)#&(flags2<1)&(cflux0>0)&(sky>0)&(ecflux0>0)
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
        
        with h5py.File(fLC, 'r') as f, h5py.File(red, 'r+') as g:
            
            lc = f['data']
            rc = g['data']
            
            try:
                grp = g['data2']
            except:
                grp = g.create_group('data2')
            else:
                if ind == 0:
                    del g['data2']
                    grp = g.create_group('data2')
           
            for i in range(len(ascc_tmp)):
        
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
