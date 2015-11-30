#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from coarse_decor import find_sigma


def spatial_decor(ind1, ind2, ind3, mag, emag, x, y, use_weights=False, maxiter=100, eps=1e-3, verbose=True):
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars3 = 4*len(np.unique(ind3))
    npars = npars1 + npars2 + npars3
    
    # Create the necessary arrays.
    weights = 1/emag**2
    snx = np.sin(2*np.pi*x)
    csx = np.cos(2*np.pi*x)
    sny = np.sin(2*np.pi*y)
    csy = np.cos(2*np.pi*y)
    z = np.zeros(np.amax(ind2) + 1)
    b = np.zeros(np.amax(ind3) + 1)
    d = np.zeros(np.amax(ind3) + 1)
    sol2 = 0.
    sol3 = 0.
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Computation of parameters.
        m = np.bincount(ind1, weights*(mag - z[ind2] - sol2 - sol3))/np.bincount(ind1, weights)
        z = np.bincount(ind2, weights*(mag - m[ind1] - sol2 - sol3))/np.bincount(ind2, weights)
    
        sol1 = m[ind1] + z[ind2]
    
        a = np.bincount(ind3, weights*(mag - sol1 - b[ind3]*csx - sol3)*snx)/np.bincount(ind3, weights*snx**2)
        b = np.bincount(ind3, weights*(mag - sol1 - a[ind3]*snx - sol3)*csx)/np.bincount(ind3, weights*csx**2)
        
        sol2 = a[ind3]*snx + b[ind3]*csx
    
        c = np.bincount(ind3, weights*(mag - sol1 - sol2 - d[ind3]*csy)*sny)/np.bincount(ind3, weights*sny**2)
        d = np.bincount(ind3, weights*(mag - sol1 - sol2 - c[ind3]*sny)*csy)/np.bincount(ind3, weights*csy**2)
    
        sol3 = c[ind3]*sny + d[ind3]*csy
    
        if use_weights:
            res = mag - sol1 - sol2 - sol3
            sigma = find_sigma(ind1, res, emag)
            weights = 1./(emag**2 + (sigma**2)[ind1])
    
        if (niter > 0):
            
            crit1 = np.nanmax(np.abs(sol1 - sol1_old))
            crit2 = np.nanmax(np.abs(sol2 - sol2_old))
            crit3 = np.nanmax(np.abs(sol3 - sol3_old))
            
            if verbose:
                print ' crit1 = %g, crit2 = %g, crit3 = %g'%(crit1, crit2, crit3)
            
            if (crit1 < eps) & (crit2 < eps) & (crit3 < eps):
                break
        
        sol1_old = np.copy(sol1)
        sol2_old = np.copy(sol2)
        sol3_old = np.copy(sol3)
    
    # Compute the chi-square value of the solution.
    chi_tmp = weights*(mag - sol1 - sol2 - sol3)**2
    chisq = np.sum(chi_tmp)/(npoints - npars)
    
    A = np.vstack([a, b, c, d]).T
    
    if use_weights:
        return m, z, sigma, A, niter, chisq, npoints, npars
    return m, z, A, niter, chisq, npoints, npars


def temporal_decor(ind1, ind2, mag, emag, use_weights=False, maxiter=100, eps=1e-3, verbose=True):
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    # Create the necessary arrays.
    weights = 1/emag**2
    z = np.zeros(np.amax(ind2) + 1)
    sigma2 = np.zeros(np.amax(ind2) + 1)
    
    crit1 = np.zeros(maxiter)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Computation of parameters.
        m = np.bincount(ind1, weights*(mag - z[ind2]))/np.bincount(ind1, weights)
        z = np.bincount(ind2, weights*(mag - m[ind1]))/np.bincount(ind2, weights)
    
        sol1 = m[ind1] + z[ind2]
    
        if use_weights:
            res = mag - sol1
            sigma1 = find_sigma(ind1, res, np.sqrt(emag**2 + (sigma2**2)[ind2]))
            sigma2 = find_sigma(ind2, res, np.sqrt(emag**2 + (sigma1**2)[ind1]))
            weights = 1./(emag**2 + (sigma1**2)[ind1] + (sigma2**2)[ind2])
    
        if (niter > 0):
            
            crit1[niter] = np.nanmax(np.abs(sol1 - sol1_old))
            
            if verbose:
                print ' crit1 = %g'%(crit1[niter])
            
            if (crit1[niter] < eps):
                break
                
        if (niter > 1):
            
            dcrit = crit1[niter] - crit1[niter-1]
            
            if verbose:
                print ' dcrit1 = %g'%(dcrit)
            
            if (abs(dcrit) < eps**2):
                break
        
        sol1_old = np.copy(sol1)
    
    # Compute the chi-square value of the solution.
    chi_tmp = weights*(mag - sol1)**2
    chisq = np.sum(chi_tmp)/(npoints - npars)
    
    if use_weights:
        return m, z, sigma1, sigma2, niter, chisq, npoints, npars
    return m, z, niter, chisq, npoints, npars


