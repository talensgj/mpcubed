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
    
    select = range(len(values))
    for i in range(niter):
        mean = np.mean(values[select])
        std = np.std(values[select])
        select = np.argwhere(np.abs(values-mean) < 2*std)
        
    return mean

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
            result[i] = np.median(values[indices == i])
        #result = result[a]
        
    if keeplength:
        return np.repeat(result, flatcount[a])
    else:
        return result

def sysrem(ind1, ind2, values, errors, a2=None, maxiter=50, eps=1e-3, idx=None):
    
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    if a2 is None:
        a2 = np.ones(np.amax(ind2)+1) # WRONG!!!
    
    s = values/errors**2
    r = 1./errors**2
    
    niter = 0
    end_crit = True
    while (niter < maxiter) & (end_crit):
        
        print niter
        
        a1 = np.bincount(ind1, s*a2[ind2])/np.bincount(ind1, r*(a2**2)[ind2])
        a2 = np.bincount(ind2, s*a1[ind1])/np.bincount(ind2, r*(a1**2)[ind1])
        
        if niter == 0:
            end_crit = True
        else:
            end_crit = (np.nanmax(np.abs((a1o-a1)/a1o)) > eps) | (np.nanmax(np.abs((a2o-a2)/a2o)) > eps)
            print np.nanmax(np.abs((a1o-a1)/a1o)), np.nanmax(np.abs((a2o-a2)/a2o))
        
        a1o = np.copy(a1)
        a2o = np.copy(a2)
            
        niter += 1
    
    chisq = np.sum(r*(values-a1[ind1]*a2[ind2])**2)/(npoints-npars)
    chisq_pbin1 = np.bincount(ind1, r*(values-a1[ind1]*a2[ind2])**2)
    chisq_pbin2 = np.bincount(ind2, r*(values-a1[ind1]*a2[ind2])**2)
    
    return a1, a2, niter, chisq, chisq_pbin1, chisq_pbin2, npoints, npars
    
    
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
        
        data = f['data']
        for i in range(len(ascc)):
            lst[select[i]:select[i+1]] = data[ascc[i]]['lst']
            flux0[select[i]:select[i+1]] = data[ascc[i]]['flux0']
            eflux0[select[i]:select[i+1]] = data[ascc[i]]['eflux0']
            sky[select[i]:select[i+1]] = data[ascc[i]]['sky']
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']
    
    dec_1 = np.copy(dec)
    
    print 'Done reading.'

    # Sky coordinates of data.
    ha = np.mod(lst*15.-np.repeat(ra,Nobs), 360.)
    dec = np.repeat(dec, Nobs)

    # Create the indices.
    hg = PolarGrid(1350, 720)
    binnum, count = hg.find_gridpoint(ha, dec)
    star_id = np.repeat(np.arange(len(ascc)), Nobs)
    
    print np.percentile(count[count>0], 5), np.percentile(count[count>0], 95)
    
    # Remove bad data.
    here, = np.where((flags < 1)&(flux0>0)&(sky>0))
    flux0 = flux0[here]
    eflux0 = eflux0[here]
    binnum = binnum[here]
    star_id = star_id[here]

    # Compute the transmission using sysrem.
    idx = np.digitize(dec_1, np.linspace(-90,90,721))
    a2, a1, niter, chisq, chisq_pbin2, chisq_pbin1, npoints, npars = sysrem(binnum, star_id, flux0, eflux0, a2=1e7*10**(vmag/-2.5))
    
    trans = hg.put_values_on_grid(a2)
    chisq_map = hg.put_values_on_grid(chisq_pbin2, fill_value=np.nan)
    stars = np.unique(star_id)
    
    median = index_statistic(idx[stars], a1[stars]/(1e7*10**(vmag[stars]/-2.5)), statistic='median', keeplength=False)
    
    plt.subplot(211)
    plt.semilogy(dec_1[stars], a1[stars]/(1e7*10**(vmag[stars]/-2.5)), '.', alpha=.1)
    plt.xlim(np.amin(dec), np.amax(dec))
    plt.xlabel('Dec [deg]')
    plt.ylabel('F/V')
    
    a1[stars] = a1[stars]/median[idx[stars]]
    
    plt.subplot(212)
    plt.semilogy(dec_1[stars], a1[stars]/(1e7*10**(vmag[stars]/-2.5)), '.', alpha=.1)
    plt.xlim(np.amin(dec), np.amax(dec))
    plt.xlabel('Dec [deg]')
    plt.ylabel('F/V')
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    print '1'
    trans[:,np.unique(idx)] = trans[:,np.unique(idx)]*median[np.unique(idx)]
    print '2'
    ax = plt.subplot(221)
    plt.imshow(trans[1:-1,1:-1].T, interpolation='None', aspect='auto', origin='lower', vmin=0, vmax=1.5, extent=(0,360,-90,90))
    plt.colorbar().set_label('Transmission')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec), np.amax(dec))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    print '3'
    plt.subplot(222, sharex=ax, sharey=ax)
    plt.imshow(count[1:-1,1:-1].T, interpolation='None', aspect='auto', origin='lower', extent=(0,360,-90,90))
    plt.colorbar().set_label('# points')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec), np.amax(dec))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    print '4'
    plt.subplot(223, sharex=ax, sharey=ax)
    plt.imshow(chisq_map[1:-1,1:-1].T, interpolation='None', aspect='auto', origin='lower', vmin=0, vmax=np.percentile(chisq_map[np.isfinite(chisq_map)], 95), extent=(0,360,-90,90))
    plt.colorbar().set_label('chisq')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec), np.amax(dec))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    print '5'
    plt.subplot(224)
    plt.semilogy(vmag[stars], chisq_pbin1[stars], '.', alpha=.1)
    plt.xlim(8.4, 2)
    plt.ylim(1e-2, 1e5)
    plt.xlabel('chisq')
    plt.ylabel('F [ADU]')
    print '6'
    plt.tight_layout()
    plt.savefig('0203LPE_maps.png')
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
