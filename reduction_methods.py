#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def sysrem(r_ij, s_ij, a_j=None, maxiter=50, eps=1e-3):
    """
    Computes the best linear fit of the from r_ij = c_i*a_j to an array of data
    r_ij with errors s_ij.
    """
    
    # If no initial estimate of a_j is given try a_j = 1.
    if a_j is None:
        a_j = np.ones(r_ij.shape[1])
    
    # Compute the initial chi2 value.
    chi2_new = np.nansum(r_ij**2./s_ij**2.)
    chi2_old = 1.1*chi2_new+2.*eps
    
    # Iterate until the improvement of the chi2 falls below the threshold or the maximum number of iterations is reached.
    niter = 0
    while (np.abs(chi2_new-chi2_old) > eps) & (niter < maxiter):
        
        chi2_old = chi2_new
        
        # Compute the best linear model for the data.
        c_i = np.nansum(r_ij*a_j/s_ij**2., axis=1)/np.nansum(a_j**2./s_ij**2., axis=1)
        a_j = np.nansum(r_ij*c_i[:,np.newaxis]/s_ij**2., axis=0)/np.nansum(c_i[:,np.newaxis]**2./s_ij**2., axis=0)
        
        # Compute the current chi2 value.
        res = r_ij-np.outer(c_i,a_j)
        chi2_new = np.nansum(res**2./s_ij**2.)
    
        niter += 1
    
    return c_i, a_j, niter, chi2_new
