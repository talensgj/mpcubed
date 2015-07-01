#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from reduction_methods import sysrem

def make_array(data, points_per_row, column_index, fill_value=np.nan):
    
    # Get the dimensions of the array.
    nrows = len(points_per_row)
    ncols = np.amax(column_index) + 1
    
    # Get the indexing for making the array.
    select = np.append(0, np.cumsum(points_per_row))
    
    # Make the array.
    array = np.full((nrows, ncols), fill_value=fill_value)
    for i in range(nrows):
        array[i, column_index[select[i]:select[i+1]]] = data[select[i]:select[i+1]] # Maybe use np.s_ here?
    
    return array
    
def make_transmission_map(data, error, dec_id, fill_value=np.nan): # Need also optional parameters for sysrem...
    
    # Create Declination bins.
    stars_per_bin = np.bincount(dec_id)
    here, = np.where(stars_per_bin > 0)
    stars_per_bin = stars_per_bin[here]  
    
    # Get the indexing for making the transmission map.
    select = np.append(0, np.cumsum(stars_per_bin))
    
    # Partition the array according to the declination bins.
    sort = np.argpartition(dec_id, select[1:-1])
    isort = np.arange(len(dec_id), dtype='int')[sort]
    data = data[sort]
    error = error[sort]
    
    nbins = len(stars_per_bin)
    
    transmission = np.full((nbins, data.shape[1]), fill_value=fill_value)
    flags = np.full((nbins, data.shape[1]), fill_value=2)
    chi2 = np.zeros(nbins)
    npoints = np.zeros(nbins)
    npars = np.zeros(nbins)
    niter = np.zeros(nbins)
    
    for i in range(nbins):
        
        # Select stars in the current declination bin.
        tmp = data[select[i]:select[i+1]] # Maybe use np.s_ here?
        etmp = error[select[i]:select[i+1]]
        
        npoints[i] = np.sum(np.isfinite(tmp))
        
        # Flag data where there are <=5 data points in a column.
        count = np.sum(np.isfinite(tmp), axis=0)
        flags[i] = np.where((count <= 5)&(count > 0), 1, flags[i])
        flags[i] = np.where((count > 5), 0, flags[i])
        
        # Perform of a single iteration of sysrem.
        c_i, a_j, niter[i], chi2[i] = sysrem(tmp, etmp)
        transmission[i] = a_j
        
        npars[i] = np.sum(np.isfinite(c_i)) + np.sum(np.isfinite(a_j))
        
    return transmission, flags, chi2, npoints, npars, niter
