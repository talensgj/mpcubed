#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def sysrem(ind1, ind2, values, errors, a2=None, maxiter=500, eps=1e-3, verbose=False):
    
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    weights = 1./errors**2
    
    if a2 is None:
        a2 = np.ones(np.amax(ind2)+1)
    
    for niter in range(maxiter):
        
        a1 = np.bincount(ind1, values*a2[ind2]*weights)/np.bincount(ind1, (a2**2)[ind2]*weights)
        a2 = np.bincount(ind2, values*a1[ind1]*weights)/np.bincount(ind2, (a1**2)[ind1]*weights)
        
        if (niter > 0):
        
            crit1 = np.nanmax(np.abs((a1_old-a1)/a1_old))
            crit2 = np.nanmax(np.abs((a2_old-a2)/a2_old))
            
            if (crit1 < eps) & (crit2 < eps):
                break
            
        if verbose:
            print 'niter = %i'%niter
        
        a1_old = np.copy(a1)
        a2_old = np.copy(a2)
    
    chi_tmp = (values-a1[ind1]*a2[ind2])**2*weights
    chisq = np.sum(chi_tmp)/(npoints-npars)
    chisq_pbin1 = np.bincount(ind1, chi_tmp)
    chisq_pbin2 = np.bincount(ind2, chi_tmp)
    
    return a1, a2, niter, chisq, chisq_pbin1, chisq_pbin2, npoints, npars
