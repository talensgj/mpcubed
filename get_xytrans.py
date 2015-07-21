#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob

import h5py
import numpy as np

import healpy
from scipy.stats import binned_statistic_2d

from time import time

import matplotlib.pyplot as plt

from coordinate_grids import HealpixGrid, PolarGrid

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

def sigma_clip(values, niter=5):
    
    arr = np.sort(values)
    n = len(values)
    
    return np.mean(arr[n*.05:n*.95])
    

def index_statistic(indices, values, statistic='mean', keeplength=True):
    
    indices = indices.astype('int')
    
    known_statistics = ['mean', 'std', 'count', 'sum', 'median']
    if statistic not in known_statistics:
        raise ValueError('invalid statistic %r' % (statistic,))
    
    flatcount = np.bincount(indices, None)
    a = flatcount.nonzero()
    
    if statistic == 'mean':
        flatsum = np.bincount(indices, values)
        result = flatsum[a] / flatcount[a]
        
    elif statistic == 'std':
        flatsum = np.bincount(indices, values)
        flatsum2 = np.bincount(indices, values ** 2)
        result = np.sqrt(flatsum2[a] / flatcount[a] - (flatsum[a] / flatcount[a]) ** 2)
        
    elif statistic == 'count':
        result = flatcount[a]
        
    elif statistic == 'sum':
        flatsum = np.bincount(indices, values)
        result = flatsum[a]
        
    elif statistic == 'median':
        result = np.zeros(len(flatcount))
        for i in np.unique(indices):
            result[i] = sigma_clip(values[indices == i])
        #result = result[a]
        
    if keeplength:
        return np.repeat(result, flatcount[a])
    else:
        return result

def sysrem(ind1, ind2, values, errors, a2=None, maxiter=50, eps=1e-9, idx=None):
    
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    if a2 is None:
        a2 = np.ones(npars2) # WRONG!!!
    
    if idx is not None:
        ref = np.copy(a2)
        stars = np.unique(ind2)
    
    niter = 0
    chisq = np.zeros(maxiter)
    chisq_crit = True
    
    s = values/errors**2
    r = 1./errors**2
    
    while (niter < maxiter) & (chisq_crit):
        
        #a1 = np.bincount(ind1, s*a2[ind2])/np.bincount(ind1, r*(a2**2)[ind2])
        #a2 = np.bincount(ind2, s*a1[ind1])/np.bincount(ind2, r*(a1**2)[ind1])
        
        a1 = index_statistic(ind1, values/a2[ind2], statistic='median', keeplength=False)
        a2 = index_statistic(ind2, values/a1[ind1], statistic='median', keeplength=False)
        
        chisq[niter] = np.sum(r*(values-a1[ind1]*a2[ind2])**2)/(npoints-npars)
        
        if niter == 0:
            chisq_crit = True
        else:
            chisq_crit = np.abs(chisq[niter]-chisq[niter-1]) > eps

        print niter, chisq[niter]
        
        if idx is not None:
            ratio = a2[stars]/ref[stars]
            median = index_statistic(idx[stars], ratio, statistic='median', keeplength=False)
            a2[stars] = a2[stars]/median[idx[stars]]

            plt.plot(np.unique(idx), median[np.unique(idx)], '.', alpha=.1)
            plt.show()

        niter += 1
    
    chisq_pbin1 = np.bincount(ind1, r*(values-a1[ind1]*a2[ind2])**2)
    chisq_pbin2 = np.bincount(ind2, r*(values-a1[ind1]*a2[ind2])**2)
    
    return a1, a2, niter, chisq[niter-1], chisq_pbin1, chisq_pbin2, npoints, npars
    
    
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
        wflux0 = np.zeros(np.sum(Nobs))
        sky = np.zeros(np.sum(Nobs))
        flags = np.zeros(np.sum(Nobs))
        
        data = f['data']
        for i in range(len(ascc)):
            lst[select[i]:select[i+1]] = data[ascc[i]]['lst']
            flux0[select[i]:select[i+1]] = data[ascc[i]]['flux0']
            eflux0[select[i]:select[i+1]] = data[ascc[i]]['eflux0']
            wflux0[select[i]:select[i+1]] = index_statistic(data[ascc[i]]['lstidx']//50, data[ascc[i]]['flux0'], statistic='std')
            sky[select[i]:select[i+1]] = data[ascc[i]]['sky']
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']
    
    dec_1 = np.copy(dec)
    
    eflux0 = np.where(eflux0>wflux0, eflux0, wflux0)
    
    print 'Done reading.'

    # Sky coordinates of data.
    ha = np.mod(lst*15.-np.repeat(ra,Nobs), 360.)
    dec = np.repeat(dec, Nobs)

    # Create the indices.
    star_id = np.repeat(np.arange(len(ascc)), Nobs)
    #hg = HealpixGrid(128)
    #binnum, count = hg.find_gridpoint(ha, dec)
    
    hg = PolarGrid(1350, 720)
    binnum, count = hg.find_gridpoint(ha, dec)
    
    print np.sum(count>0)
    print np.percentile(count[count>0], 5), np.percentile(count[count>0], 95)

    Nok = np.bincount(star_id, flags==0)
    frac = Nok/Nobs
    
    frac = np.repeat(frac, Nobs)
    Nok = np.repeat(Nok, Nobs)

    # Remove bad data.
    here, = np.where((flags < 1)&(eflux0>0)&(flux0>0)&(frac>.95)&(Nok>=50)&(sky>0))
    flux0 = flux0[here]
    eflux0 = eflux0[here]
    binnum = binnum[here]
    star_id = star_id[here]

    # Compute the transmission using sysrem.
    idx = np.digitize(dec_1, bins=np.linspace(-90,90,721))
    a1, a2, niter, chisq, chisq_pbin1, chisq_pbin2, npoints, npars = sysrem(binnum, star_id, flux0, eflux0, a2=1e7*10**(vmag/-2.5), idx=idx)
    
    trans = hg.put_values_on_grid(a1)
    plt.imshow(trans.T, interpolation='None')
    plt.show()
    
    #healpy.mollview(trans, xsize=8000)
    #plt.show()

    
    return 0

if __name__ == '__main__':
    #import argparse
    
    #parser = argparse.ArgumentParser(description='Compute the transmission map for a particular night and camera.')
    
    #parser.add_argument('path', type=str, help='the global path for search')
    #parser.add_argument('-n', '--night', default='', type=str, help='the night to reduce')
    #parser.add_argument('-c', '--camera',default='', type=str, help='the camera to reduce')
    
    #args = parser.parse_args()
    
    filelist = glob.glob('/data2/mascara/LaPalma/20150203LPE/fLC/fLC_*.hdf5')
    #filelist = np.sort(filelist)

    #filelist = ['/data2/talens/fLC_20150203LPC.hdf5']

    for filename in filelist:
        print 'Data:', filename
        transmission(filename)
