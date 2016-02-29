#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from collections import namedtuple

Quality = namedtuple('Quality', 'niter chisq npoints npars') 

def sysrem(value, weights, maxiter=100, dtol=1e-3, verbose=True):

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
    
        print par1
        print par2
    
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

def main():
    return

if __name__ == '__main__':
    main()
