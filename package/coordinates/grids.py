#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import healpy

class PolarGrid(object):
    """ Return a polar coordinate grid."""
    
    def __init__(self, nx, ny):
        
        self.nx = nx
        self.ny = ny
        
        self.xedges = np.linspace(0, 360, self.nx+1)
        self.yedges = np.linspace(-90, 90, self.ny+1)
        
        self.npix = (nx + 2)*(ny + 2)
    
    def _ra2idx(self, ra):
        
        idx = np.searchsorted(self.xedges, ra, 'right')
        
        return idx
        
    def _dec2idx(self, dec):
        
        idx = np.searchsorted(self.yedges, dec, 'right')
        
        return idx
    
    def radec2idx(self, ra, dec):
        
        idx1 = self._ra2idx(ra)
        idx2 = self._dec2idx(dec)
        
        return idx1, idx2
            
    def _idx2ra(self, idx):
        
        ra = np.full(self.nx+2, fill_value=np.nan)
        ra[1:-1] = (self.xedges[:-1] + self.yedges[1:])/2.
        
        return ra[idx]
        
    def _idx2dec(self, idx):
        
        dec = np.full(self.ny+2, fill_value=np.nan)
        dec[1:-1] = (self.yedges[:-1] + self.yedges[1:])/2.
        
        return dec[idx]
        
    def idx2radec(self, idx1, idx2):
        
        ra = self._idx2ra(idx1)
        dec = self._idx2dec(idx2)
        
        return ra, dec
        
    def values2grid(self, idx1, idx2, values, fill_value=0):
        
        array = np.full((self.nx+2, self.ny+2), fill_value=fill_value)
        array[idx1, idx2] = values
        
        return array


class HealpixGrid(object):
    """Return a healpix coordinate grid."""
    
    def __init__(self, nside):
        
        self.nside = nside
        self.npix = healpy.nside2npix(self.nside)
    
    def radec2idx(self, ra, dec):
        
        idx = healpy.ang2pix(self.nside, (90. - dec)*np.pi/180., ra*np.pi/180.)
        
        return idx
            
    def idx2radec(self, idx):
    
        theta, phi = healpy.pix2ang(self.nside, idx)
        dec = 90. - theta*180./np.pi
        ra = phi*180./np.pi
            
        return ra, dec
        
    def values2grid(self, idx, values, fill_value=0):
        
        array = np.full(self.npix, fill_value=fill_value)
        array[idx] = values
        
        return array

    
class CartesianGrid(object):
    """Return a cartesian coordinate grid."""
    
    def __init__(self, nx, ny, margin=0, Lx=4008, Ly=2672):
        
        self.nx = nx
        self.ny = ny
        
        self.xedges = np.linspace(margin, Lx-margin, self.nx+1)
        self.yedges = np.linspace(margin, Ly-margin, self.ny+1)
        
        self.npix = (nx + 2)*(ny + 2)
        
    def _x2idx(self, x):
        
        idx = np.searchsorted(self.xedges, x, 'right')
        
        return idx

    def _y2idx(self, y):
        
        idx = np.searchsorted(self.yedges, y, 'right')
        
        return idx
        
    def xy2idx(self, x, y):
        
        idx1 = self._x2idx(x)
        idx2 = self._y2idx(y)
        
        return idx1, idx2
        
    def _idx2x(self, idx):
        
        x = np.full(self.nx+2, fill_value=np.nan)
        x[1:-1] = (self.bins1[:-1] + self.bins1[1:])/2.
        
        return x[idx]
        
    def _idx2y(self, idx):
        
        y = np.full(self.ny+2, fill_value=np.nan)
        y[1:-1] = (self.bins2[:-1] + self.bins2[1:])/2.
        
        return y[idx]
    
    def idx2xy(self, idx1, idx2):
        
        x = self._idx2x(idx1)
        y = self._idx2y(idx2)
        
        return x, y    
    
    def values2grid(self, idx1, idx2, values, fill_value=0):
        
        array = np.full((self.nx+2, self.ny+2), fill_value=fill_value)
        array[idx1, idx2] = values
        
        return array
