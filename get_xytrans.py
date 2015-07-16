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

def regrid(ind1, ind2, values, errors, a2):
    
    npoints = len(values)
    npars = len(np.unique(ind1)) + len(np.unique(ind2))
        
    a1 = np.bincount(ind1, values*a2[ind2]/errors**2)/np.bincount(ind1, (a2**2)[ind2]/errors**2)
        
    chisq = np.nansum((values-a1[ind1]*a2[ind2])**2/errors**2)/(npoints-npars)
    
    chisq_pbin1 = np.bincount(ind1, (values-a1[ind1]*a2[ind2])**2/errors**2)
    chisq_pbin2 = np.bincount(ind2, (values-a1[ind1]*a2[ind2])**2/errors**2)
    
    return a1, a2, chisq, chisq_pbin1, chisq_pbin2, npoints, npars
    
def polar_grid(ra, dec, nx, ny):
    
    ind1 = np.searchsorted(np.linspace(0, 360, nx+1), ra, 'right')
    ind2 = np.searchsorted(np.linspace(-90, 90, ny+1), dec, 'right')
    
    ind = np.ravel_multi_index([ind1,ind2], (nx+2, ny+2))
    
    count = np.bincount(ind, minlength=(nx+2)*(ny+2))
    count = count.reshape((nx+2, ny+2))
    
    return ind, count
    
def test_grid(values, nx, ny):
    
    array = np.zeros((nx+2)*(ny+2))
    array[:len(values)] = values
    
    return array
    
class PolarGrid():
    
    def __init__(self, nx, ny):
        
        self.nx = nx
        self.ny = ny
        
        self.bins1 = np.linspace(0, 360, self.nx+1)
        self.bins2 = np.linspace(-90, 90, self.ny+1)
        
    def find_gridpoint(self, ra, dec):
        
        ind1 = np.searchsorted(self.bins1, ra, 'right')
        ind2 = np.searchsorted(self.bins2, dec, 'right')
        
        ind = np.ravel_multi_index([ind1,ind2], (self.nx+2, self.ny+2))
        
        count = np.bincount(ind)
        count = self.put_values_on_grid(count)
        
        return ind, count
        
    def put_values_on_grid(self, values):
        
        array = np.zeros((self.nx+2)*(self.ny+2))
        array[:len(values)] = values
        array = array.reshape((self.nx+2, self.ny+2))
    
        return array
    
class CartesianGrid():
    
    def __init__(self, nx, ny, margin=0, Lx=4008, Ly=2672):
        
        self.nx = nx
        self.ny = ny
        
        self.bins1 = np.linspace(margin, Lx-margin, self.nx+1)
        self.bins2 = np.linspace(margin, Ly-margin, self.ny+1)
        
    def find_gridpoint(self, x, y):
        
        ind1 = np.searchsorted(self.bins1, x, 'right')
        ind2 = np.searchsorted(self.bins2, y, 'right')
        ind = np.ravel_multi_index([ind1,ind2], (self.nx+2, self.ny+2))
    
        count = np.bincount(ind)
        count = self.put_values_on_grid(count)
        
        return ind, count
        
    def put_values_on_grid(self, values):
        
        array = np.zeros((self.nx+2)*(self.ny+2))
        array[:len(values)] = values
        array = array.reshape((self.nx+2, self.ny+2))
        
        return array
    
class HealpixGrid():
    
    def __init__(self, nside):
        
        self.nside = nside
        self.npix = healpy.nside2npix(self.nside)
    
    def find_gridpoint(self, ra, dec):
        
        ind = healpy.ang2pix(self.nside, (dec+90)*np.pi/180., ra*np.pi/180.)
        
        count = np.bincount(ind)
        count = self.put_values_on_grid(count)
        
        return ind, count
        
    def put_values_on_grid(self, values):
        
        array = np.zeros(self.npix)
        array[:len(values)] = values
        
        return array
        
    
def healpy_grid(ra, dec, nside):
    
    npix = healpy.nside2npix(nside)
    ind = healpy.ang2pix(nside, (dec+90)*np.pi/180., ra*np.pi/180.)
    
    count = np.bincount(ind, minlength=npix)
    
    return ind, count
    
def cartesian_grid(x, y, nx, ny, Lx=4008, Ly=2672, margin=0):
    
    ind1 = np.searchsorted(np.linspace(margin, Lx-margin, nx+1), x, 'right')
    ind2 = np.searchsorted(np.linspace(margin, Ly-margin, ny+1), y, 'right')
    
    ind = np.ravel_multi_index([ind1,ind2], (nx+2, ny+2))
    
    count = np.bincount(ind, minlength=(nx+2)*(ny+2))
    count = count.reshape((nx+2, ny+2))
    
    return ind, count
    
    
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
    
    print 'Done reading.'

    dec_1 = np.copy(dec)

    # Sky coordinates of data.
    ha = np.mod(lst*15.-np.repeat(ra,Nobs), 360.)
    dec = np.repeat(dec, Nobs)

    # Create the indices.
    star_id = np.repeat(np.arange(len(ascc)), Nobs)
    binnum, count = healpy_grid(ha, dec, 512)
    
    hg = HealpixGrid(512)
    binnum1, count1 = hg.find_gridpoint(ha, dec)
    
    print np.percentile(count[count>0], 5), np.percentile(count[count>0], 95)
    
    # Remove bad data.
    here, = np.where(flags < 1)
    flux0 = flux0[here]
    eflux0 = eflux0[here]
    x = x[here]
    y = y[here]
    binnum = binnum[here]
    star_id = star_id[here]

    # Compute the transmission using sysrem.
    a1, a2, niter, chisq, chisq_pbin1, chisq_pbin2, npoints, npars = sysrem(binnum, star_id, flux0, eflux0, a2=1e7*10**(vmag/-2.5))
    
    
    array = np.zeros(len(count))
    array[:len(a1)] = a1
    
    healpy.mollview(array)
    plt.show()
    
    stars = np.unique(star_id)
    ratio = a2[stars]/(1e7*10**(vmag[stars]/-2.5))
    
    plt.plot(dec_1[stars], ratio, '.', alpha=.1)
    plt.show()
    
    # Recompute the transmission on a new grid.
    binnum, count = cartesian_grid(x, y, 800, 600, margin=50)
    a1, a2, chisq, chisq_pbin1, chisq_pbin2, npoints, npars = regrid(binnum, star_id, flux0, eflux0, a2)

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
