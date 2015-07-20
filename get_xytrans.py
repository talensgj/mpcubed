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

from coordinate_grids import HealpixGrid

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

def index_statistic(indices, values, statistic='mean'):
    
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
        result = result[a]
        
    return np.repeat(result, flatcount[a])

def sysrem(ind1, ind2, values, errors, a2=None, maxiter=50, eps=1e-3):
    
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    if a2 is None:
        a2 = np.ones(npars2) # WRONG!!!
    
    niter = 0
    chisq = np.zeros(maxiter)
    chisq_crit = True
    
    s = values/errors**2
    r = 1./errors**2
    
    while (niter < maxiter) & (chisq_crit):
        
        a1 = np.bincount(ind1, s*a2[ind2])/np.bincount(ind1, r*(a2**2)[ind2])
        a2 = np.bincount(ind2, s*a1[ind1])/np.bincount(ind2, r*(a1**2)[ind1])
        
        # CLIPPING?
        #res = values-a1[ind1]*a2[ind2]
        #mask = np.abs(res - np.mean(res)) < 3*np.std(res)
        
        chisq[niter] = np.sum(r*(values-a1[ind1]*a2[ind2])**2)/(npoints-npars)
        
        if niter == 0:
            chisq_crit = True
        else:
            chisq_crit = np.abs(chisq[niter]-chisq[niter-1]) > eps

        print niter, chisq[niter]

        niter += 1
    
    chisq_pbin1 = np.bincount(ind1, r*(values-a1[ind1]*a2[ind2])**2)
    chisq_pbin2 = np.bincount(ind2, r*(values-a1[ind1]*a2[ind2])**2)
    
    return a1, a2, niter, chisq[niter-1], chisq_pbin1, chisq_pbin2, npoints, npars
    
    
def transmission(filename):

    with h5py.File(filename) as f:
        
        hdr = f['table_header']
        ascc = hdr['ascc'].value
        vmag = hdr['vmag'].value
        bmag = hdr['bmag'].value
        ra = hdr['ra'].value
        dec = hdr['dec'].value
        Nobs = hdr['nobs'].value.astype('int')
        blendval = hdr['blendvalue'].value
        
        select = np.append(0, np.cumsum(Nobs))
    
        jd = np.zeros(np.sum(Nobs))
        lst = np.zeros(np.sum(Nobs))
        flux0 = np.zeros(np.sum(Nobs))
        eflux0 = np.zeros(np.sum(Nobs))
        sky = np.zeros(np.sum(Nobs))
        esky = np.zeros(np.sum(Nobs))
        x = np.zeros(np.sum(Nobs))
        y = np.zeros(np.sum(Nobs))
        flags = np.zeros(np.sum(Nobs))
        
        ratio = np.zeros(len(Nobs))
        
        data = f['data']
        for i in range(len(ascc)):
            jd[select[i]:select[i+1]] = data[ascc[i]]['jdmid']
            lst[select[i]:select[i+1]] = data[ascc[i]]['lst']
            flux0[select[i]:select[i+1]] = data[ascc[i]]['flux0']
            #eflux0[select[i]:select[i+1]] = data[ascc[i]]['eflux0']
            eflux0[select[i]:select[i+1]] = index_statistic(data[ascc[i]]['lstidx']//50, data[ascc[i]]['flux0'], statistic='std')
            sky[select[i]:select[i+1]] = data[ascc[i]]['sky']
            esky[select[i]:select[i+1]] = data[ascc[i]]['esky']
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']
            x[select[i]:select[i+1]] = data[ascc[i]]['x']
            y[select[i]:select[i+1]] = data[ascc[i]]['y']
            ratio[i] = np.mean(data[ascc[i]]['flux1']/data[ascc[i]]['flux0'])
    
    print 'Done reading.'

    dec_1 = np.copy(dec)

    # Sky coordinates of data.
    ha = np.mod(lst*15.-np.repeat(ra,Nobs), 360.)
    dec = np.repeat(dec, Nobs)

    # Create the indices.
    star_id = np.repeat(np.arange(len(ascc)), Nobs)
    hg = HealpixGrid(512)
    binnum, count = hg.find_gridpoint(ha, dec)
    
    print np.percentile(count[count>0], 5), np.percentile(count[count>0], 95)
    
    Nok = np.bincount(star_id, flags==0)
    print Nok/Nobs
    frac = np.repeat(Nok/Nobs, Nobs)
    Nok = np.repeat(Nok, Nobs)
    print len(flux0)
    # Remove bad data.
    here, = np.where((flags < 1)&(eflux0>0)&(esky/sky<.1))
    jd = jd[here]
    flux0 = flux0[here]
    eflux0 = eflux0[here]
    x = x[here]
    y = y[here]
    binnum = binnum[here]
    star_id = star_id[here]

    print len(flux0)

    # Compute the transmission using sysrem.
    a1, a2, niter, chisq, chisq_pbin1, chisq_pbin2, npoints, npars = sysrem(binnum, star_id, flux0, eflux0, a2=1e7*10**(vmag/-2.5))

    trans = hg.put_values_on_grid(a1)
    
    healpy.mollview(trans)
    healpy.graticule(dpar=5)
    plt.show()

    plt.plot(ratio[np.unique(star_id)], a2[np.unique(star_id)]/(1e7*10**(vmag[np.unique(star_id)]/-2.5)), '.', alpha=.1)
    plt.show()

    plt.plot(dec_1[np.unique(star_id)], a2[np.unique(star_id)]/(1e7*10**(vmag[np.unique(star_id)]/-2.5)), '.', alpha=.1)
    plt.show()

    fit = a1[binnum]*a2[star_id]
    
    stars = np.unique(star_id)
    Nok = np.bincount(star_id)
    
    
    ascc = ascc[stars]
    Nobs = Nobs[stars]
    Nok = Nok[stars]
    vmag = vmag[stars]
    a2 = a2[stars]
    chisq_pbin2 = chisq_pbin2[stars]
    
    stat1 = a2/(1e7*10**(vmag/-2.5))
    stat2 = chisq_pbin2/Nok
    
    args = np.argsort(stat1)
    args = args[::-1]
    
    for i in args:
        print i
        if stat1[i] < 1.5: break
        print 'OK'
        here = star_id==stars[i]
        print 'OK'
        plt.subplot(211)
        plt.title('ASCC %s, Nobs = %i, Nok = %i, V = %.2f, F/V = %.2f, chisq/Nok= %.2f'%(ascc[i], Nobs[i], Nok[i], vmag[i], stat1[i], stat2[i])) 
        plt.errorbar(jd[here], flux0[here], yerr=eflux0[here], fmt='.')
        plt.plot(jd[here], fit[here])
        print 'OK'
        plt.subplot(212)
        plt.errorbar(jd[here], flux0[here]-fit[here], yerr=eflux0[here], fmt='.')
        print 'OK'
        #plt.tight_layout()
        plt.savefig('0203LPE/wsigmasky/ascc%s.png'%ascc[i])
        plt.close()
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
