#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob

import h5py
import numpy as np

import healpy
from scipy.stats import binned_statistic_2d

import matplotlib.pyplot as plt

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

def sysrem(ind1, ind2, values, errors, a1=None, a2=None):
    
    return a1, a2

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
    
        lst = np.zeros(np.sum(Nobs))
        flux0 = np.zeros(np.sum(Nobs))
        eflux0 = np.zeros(np.sum(Nobs))
        x = np.zeros(np.sum(Nobs))
        y = np.zeros(np.sum(Nobs))
        flags = np.zeros(np.sum(Nobs))
        
        data = f['data']
        for i in range(len(ascc)):
            lst[select[i]:select[i+1]] = data[ascc[i]]['lst']
            flux0[select[i]:select[i+1]] = data[ascc[i]]['flux0']
            eflux0[select[i]:select[i+1]] = data[ascc[i]]['eflux0']
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']
            x[select[i]:select[i+1]] = data[ascc[i]]['x']
            y[select[i]:select[i+1]] = data[ascc[i]]['y']

    dec_1 = np.copy(dec)

    # 
    ha = np.mod(lst*15.-np.repeat(ra,Nobs), 360.)
    dec = np.repeat(dec, Nobs)

    nside = 1024
    filename = 'Tvmag_20150203LPE_ns%i.hdf5'%nside

    # Create the indices.
    npix = healpy.nside2npix(nside)
    binnum = healpy.ang2pix(nside, (dec+90)*np.pi/180., ha*np.pi/180.)
    star_id = np.repeat(np.arange(len(ascc)), Nobs)
    
    Ngood = np.bincount(star_id, flags<1)
    plt.plot(Nobs, Ngood, '.')
    plt.show()
    
    Ngood = np.repeat(Ngood, Nobs)
    
    # Remove bad data.
    here, = np.where((flags < 1)&np.isfinite(eflux0)&(Ngood>150))
    flux0 = flux0[here]
    eflux0 = eflux0[here]
    binnum = binnum[here]
    star_id = star_id[here]
    x = x[here]
    y = y[here]
    
    bins = np.unique(binnum)
    stars = np.unique(star_id)
    
    npoints = len(flux0)
    npars = len(bins)+len(stars)
    
    maxiter = 50
    eps = 1e-3
    niter = 0
    chisq_crit = True
    a_j = np.ones(npix)
    chisq = np.zeros(maxiter)
    c_i = 1e7*10**(vmag/-2.5)
    while (niter < maxiter) & (chisq_crit):
        
        a_j = np.bincount(binnum, flux0*c_i[star_id]/eflux0**2)/np.bincount(binnum, (c_i**2)[star_id]/eflux0**2)
        c_i = np.bincount(star_id, flux0*a_j[binnum]/eflux0**2)/np.bincount(star_id, (a_j**2)[binnum]/eflux0**2)
        
        chisq[niter] = np.nansum((flux0-c_i[star_id]*a_j[binnum])**2/eflux0**2)/(npoints-npars)
        
        if niter == 0:
            chisq_crit = True
        else:
            chisq_crit = np.abs(chisq[niter]-chisq[niter-1]) > eps

        print niter, chisq[niter]

        niter += 1
    
    chisq_pbin = np.bincount(binnum, (flux0-c_i[star_id]*a_j[binnum])**2/eflux0**2)
        
    trans = np.full(npix, fill_value=np.nan)
    trans[bins] = a_j[bins]

    chisq = np.full(npix, fill_value=np.nan)
    chisq[bins] = chisq_pbin[bins]
        
    count = np.full(npix, fill_value=np.nan)
    count[bins] = np.bincount(binnum)[bins]
    
    print np.nanmean(count)
    
    f = h5py.File(filename)
    f.create_dataset('ra', data=ra[stars])
    f.create_dataset('ascc', data=ascc[stars])
    f.create_dataset('dec', data=dec_1[stars])
    f.create_dataset('nobs', data=Nobs[stars])
    f.create_dataset('vmag', data=vmag[stars])
    f.create_dataset('bmag', data=bmag[stars])
    f.create_dataset('blendvalue', data=blendval[stars])
    f.create_dataset('flux', data=c_i[stars])
    f.create_dataset('count', data=count)
    f.create_dataset('trans', data=trans)
    f.create_dataset('chisq', data=chisq)
    f.close()
    
    #  
    nx = (200, 400, 800)
    ny = (150, 300, 600)  
    for nx, ny in zip(nx, ny):
        print nx, ny
    
        count_xy, xedges, yedges, binnum = binned_statistic_2d(x, y, flux0, statistic='count', bins = [np.linspace(50,3958,nx+1), np.linspace(50,2622,ny+1)])
        a_j = np.bincount(binnum, flux0*c_i[star_id]/eflux0**2)/np.bincount(binnum, (c_i**2)[star_id]/eflux0**2)
        chisq_pbin = np.bincount(binnum, (flux0-c_i[star_id]*a_j[binnum])**2/eflux0**2)
        
        bins = np.unique(binnum)
        npars = len(bins)+len(stars)
        
        trans_xy = np.full((nx+2)*(ny+2), fill_value=np.nan)
        trans_xy[bins] = a_j[bins]
        trans_xy = trans_xy.reshape((ny+2,nx+2))
        
        chisq_xy = np.full((nx+2)*(ny+2), fill_value=np.nan)
        chisq_xy[bins] = chisq_pbin[bins]
        chisq_xy = chisq_xy.reshape((ny+2,nx+2))
        
        count_xy = count_xy.T
            
        print np.nansum(chisq_xy)/(npoints-npars)
        
        f = h5py.File(filename)
        grp = f.create_group('nx%iny%i'%(nx,ny))
        grp.create_dataset('count', data=count_xy)
        grp.create_dataset('trans', data=trans_xy)
        grp.create_dataset('chisq', data=chisq_xy)
        f.close()

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
