#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import healpy

import matplotlib.pyplot as plt

class PolarGrid():
    
    def __init__(self, nx, ny):
        
        self.nx = nx
        self.ny = ny
        
        self.bins1 = np.linspace(0, 360, self.nx+1)
        self.bins2 = np.linspace(-90, 90, self.ny+1)
        
    def find_gridpoint(self, ra, dec):
        
        ind1 = np.searchsorted(self.bins1, ra, 'right')
        ind2 = np.searchsorted(self.bins2, dec, 'right')
        
        if np.any(ind1==0) | np.any(ind1==self.nx+1):
            print 'Warning 1'
        if np.any(ind2==0) | np.any(ind2==self.ny+1):
            print 'Warning 2'
        
        ind = np.ravel_multi_index([ind1,ind2], (self.nx+2, self.ny+2))
        
        count = np.bincount(ind)
        count = self.put_values_on_grid(count)
        
        return ind, count
        
    def put_values_on_grid(self, values, fill_value=0):
        
        array = np.full((self.nx+2)*(self.ny+2), fill_value=fill_value)
        args = np.where(values>0)
        array[args] = values[args]
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
        
        if np.any(ind1==0) | np.any(ind1==self.nx+1):
            print 'Warning 1'
        if np.any(ind2==0) | np.any(ind2==self.ny+1):
            print 'Warning 2'
        
        ind = np.ravel_multi_index([ind1,ind2], (self.nx+2, self.ny+2))
    
        count = np.bincount(ind)
        count = self.put_values_on_grid(count)
        
        return ind, count
        
    def put_values_on_grid(self, values):
        
        array = np.zeros((self.nx+2)*(self.ny+2))
        array[len(values)] = values
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
