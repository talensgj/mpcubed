#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import healpy

class PolarGrid():
    
    # Does not deal with out of range bins correctly.
    
    def __init__(self, nx, ny):
        
        self.nx = nx
        self.ny = ny
        
        self.bins1 = np.linspace(0, 360, self.nx+1)
        self.bins2 = np.linspace(-90, 90, self.ny+1)
    
    def find_raidx(self, ra, compact=False):
        
        ind = np.searchsorted(self.bins1, ra, 'right')
        
        if not compact:
            return ind
        else: 
            return np.unique(ind, return_inverse=True)
        
    def find_decidx(self, dec, compact=False):
        
        ind = np.searchsorted(self.bins2, dec, 'right')
        
        if not compact:
            return ind
        else: 
            return np.unique(ind, return_inverse=True)
    
    def find_gridpoint(self, ra, dec, compact=False):
        
        ind1 = self.find_raidx(ra)
        ind2 = self.find_decidx(dec)
        
        if np.any(ind1==0) | np.any(ind1==self.nx+1):
            print 'Warning 1'
        if np.any(ind2==0) | np.any(ind2==self.ny+1):
            print 'Warning 2'
        
        ind = np.ravel_multi_index([ind1,ind2], (self.nx+2, self.ny+2))
        
        if not compact:
            return ind
        else:
            return np.unique(ind, return_inverse=True)
        
    def find_ra(self, raidx):
        
        ra = np.full(self.nx+2, fill_value=np.nan)
        ra[1:-1] = (self.bins1[:-1]+self.bins1[1:])/2.
        
        return ra[raidx]
        
    def find_dec(self, decidx):
        
        dec = np.full(self.ny+2, fill_value=np.nan)
        dec[1:-1] = (self.bins2[:-1]+self.bins2[1:])/2.
        
        return dec[decidx]
        
    def grid_coordinates(self, ind=None):
        
        ra = (self.bins1[:-1]+self.bins1[1:])/2.
        dec = (self.bins2[:-1]+self.bins2[1:])/2.
        
        ra, dec = np.meshgrid(ra, dec)
        ra = np.ravel(ra)
        dec = np.ravel(dec)
        
        if ind is None:
            return ra, dec
        else:
            return ra[ind], dec[ind]
            
    def put_values_on_grid(self, values, ind=None, fill_value=0):
        
        array = np.full((self.nx+2)*(self.ny+2), fill_value=fill_value)
        
        if ind is None:
            array[:len(values)] = values
        else:
            array[ind] = values
        
        array = array.reshape((self.nx+2, self.ny+2))
        
        return array

    
class CartesianGrid():
    
    def __init__(self, nx, ny, margin=0, Lx=4008, Ly=2672):
        
        self.nx = nx
        self.ny = ny
        
        self.bins1 = np.linspace(margin, Lx-margin, self.nx+1)
        self.bins2 = np.linspace(margin, Ly-margin, self.ny+1)
        
    def find_xidx(self, x, compact=False):
        
        ind = np.searchsorted(self.bins1, x, 'right')
        
        if not compact:
            return ind
        else: 
            return np.unique(ind, return_inverse=True)
            
    def find_yidx(self, y, compact=False):
        
        ind = np.searchsorted(self.bins2, y, 'right')
        
        if not compact:
            return ind
        else: 
            return np.unique(ind, return_inverse=True)
        
    def find_gridpoint(self, x, y, compact=False):
        
        ind1 = self.find_xidx(x)
        ind2 = self.find_yidx(y)
        
        if np.any(ind1==0) | np.any(ind1==self.nx+1):
            print 'Warning 1'
        if np.any(ind2==0) | np.any(ind2==self.ny+1):
            print 'Warning 2'
        
        ind = np.ravel_multi_index([ind1,ind2], (self.nx+2, self.ny+2))
    
        if not compact:
            return ind
        else:
            return np.unique(ind, return_inverse=True)
        
    def xidx2x(self, xidx):
        
        x = np.full(self.nx+2, fill_value=np.nan)
        x[1:-1] = (self.bins1[:-1]+self.bins1[1:])/2.
        
        return x[xidx]
        
    def yidx2y(self, yidx):
        
        y = np.full(self.ny+2, fill_value=np.nan)
        y[1:-1] = (self.bins2[:-1]+self.bins2[1:])/2.
        
        return y[yidx]
        
        
    def put_values_on_grid(self, values, ind=None, fill_value=0):
        
        array = np.full((self.nx+2)*(self.ny+2), fill_value=fill_value)
        
        if ind is None:
            array[:len(values)] = values
        else:
            array[ind] = values
        
        array = array.reshape((self.nx+2, self.ny+2))
        
        return array

    
class HealpixGrid():
    
    def __init__(self, nside):
        
        self.nside = nside
        self.npix = healpy.nside2npix(self.nside)
    
    def find_gridpoint(self, ra, dec, compact=False):
        
        ind = healpy.ang2pix(self.nside, (90-dec)*np.pi/180., ra*np.pi/180.)
        
        if not compact:
            return ind
        else:
            return np.unique(ind, return_inverse=True)
            
    def grid_coordinates(self, ind=None):
    
        if ind is None:
            theta, phi = healpy.pix2ang(self.nside, np.arange(self.npix))
        else:
            theta, phi = healpy.pix2ang(self.nside, ind)
            
        dec = 90-theta*180./np.pi
        ra = phi*180./np.pi
            
        return ra, dec
        
    def put_values_on_grid(self, values, ind=None, fill_value=0):
        
        array = np.full(self.npix, fill_value=fill_value)
        
        if ind is None:
            array[:len(values)] = values
        else:
            array[ind] = values
        
        return array
