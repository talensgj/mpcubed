#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from collections import namedtuple

Quality = namedtuple('Quality', 'niter chisq npoints npars') 

def sysrem(value, weights, maxiter=100, dtol=1e-3, verbose=True):
    """ Perform a coarse decorrelation.

    Args:
        value (float): A 2d-array of values.
        weights (float): A 2d-array of weights corresponding to value.
        maxiter (int): The maximum number of iterations to perform. Default
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the change falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.

    Returns:
        par1 (float): The parameters corresponding to the first axis of value.
        par2 (float): The parameters corresponding to the second axis of value.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.

    """
    
    print 'Warning: this function has not been extensively tested.'
    
    # Determine the number of datapoints and parameters to fit.
    npoints = value.size
    npars1 = value.shape[0]
    npars2 = value.shape[1]
    npars = npars1 + npars2

    # Create arrays.
    par2 = np.ones((value.shape[0], 1))
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        par1 = np.nansum(weights*value*par2, axis=1, keepdims=True)/np.nansum(weights*par2**2, axis=1, keepdims=True)
        par2 = np.nansum(weights*value*par1, axis=0, keepdims=True)/np.nansum(weights*par1**2, axis=0, keepdims=True)
    
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit1 = np.nanmax(np.abs(par1 - par1_old))
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            
            if (dcrit1 < dtol) & (dcrit2 < dtol):
                break
        
        par1_old = np.copy(par1)
        par2_old = np.copy(par2)
            
    # Compute the chi-square of the fit.
    fit = np.outer(par1, par2)
    chisq = weights*(value - fit)**2        
    chisq = np.nansum(chisq)
    
    return par1, par2, Quality(niter, chisq, npoints, npars)
