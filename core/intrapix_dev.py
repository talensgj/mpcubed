#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import minimize

def coarse_positions(ind1, ind2, ind3, x, y, mag, emag, maxiter=100, eps=1e-3, verbose=False, use_weights=True):
    """ 
    Given data this code will fit a model of the form:
    mag = m[ind1] + z[ind2] + a[ind3]*np.sin(2*np.pi*y) + b[ind3]*np.cos(2*np.pi*y)
    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars3 = 2*len(np.unique(ind3))
    npars = npars1 + npars2 + npars3
    
    # Create the necessary arrays.
    weights = 1./emag**2
    sny = np.sin(2*np.pi*y)
    csy = np.cos(2*np.pi*y)
    snx = np.sin(2*np.pi*x)
    csx = np.cos(2*np.pi*x)
    z = np.zeros(np.amax(ind2)+1)
    b = np.zeros(np.amax(ind3)+1)
    d = np.zeros(np.amax(ind3)+1)
    sol2 = np.zeros(len(mag))
    sol3 = np.zeros(len(mag))
    
    # Start the iterative solution.
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Computation of parameters.
        m = np.bincount(ind1, weights*(mag-z[ind2]-sol2-sol3))/np.bincount(ind1, weights)
        z = np.bincount(ind2, weights*(mag-m[ind1]-sol2-sol3))/np.bincount(ind2, weights)
        
        sol1 = m[ind1] + z[ind2]
        
        a = np.bincount(ind3, weights*(mag-sol1-b[ind3]*csy-sol3)*sny)/np.bincount(ind3, weights*sny**2)
        b = np.bincount(ind3, weights*(mag-sol1-a[ind3]*sny-sol3)*csy)/np.bincount(ind3, weights*csy**2)
        
        sol2 = a[ind3]*sny + b[ind3]*csy
        
        c = np.bincount(ind3, weights*(mag-sol1-sol2-d[ind3]*csx)*snx)/np.bincount(ind3, weights*snx**2)
        d = np.bincount(ind3, weights*(mag-sol1-sol2-c[ind3]*snx)*csx)/np.bincount(ind3, weights*csx**2)
        
        sol3 = c[ind3]*snx + d[ind3]*csx
        
        if use_weights:
            res = mag - sol1 - sol2 - sol3
            sigma = find_sigma(ind1, res, emag)
            weights = 1./(emag**2+(sigma**2)[ind1])
        
        if (niter > 0):
            
            # Check if the solution has converged.
            critm = np.nanmax(np.abs(m-m_old))
            critz = np.nanmax(np.abs(z-z_old))
            crita = np.nanmax(np.abs(a-a_old))
            critb = np.nanmax(np.abs(b-b_old))
            critc = np.nanmax(np.abs(c-c_old))
            critd = np.nanmax(np.abs(d-d_old))
            
            if verbose:
                print ' critm = %g, critz = %g, crita = %g, critb = %g'%(critm, critz, crita, critb)
            
            if (critm < eps) & (critz < eps) & (crita < eps) & (critb < eps) & (critc < eps) & (critd < eps):
                break
        
        m_old = np.copy(m)
        z_old = np.copy(z)
        a_old = np.copy(a)
        b_old = np.copy(b)
        c_old = np.copy(c)
        d_old = np.copy(d)
    
    # Compute the chi-square value of the solution.
    chi_tmp = weights*(mag - m[ind1] - z[ind2] - a[ind3]*sny - b[ind3]*csy - c[ind3]*snx - d[ind3]*csx)**2
    chisq = np.sum(chi_tmp)/(npoints-npars)
    
    if use_weights:
        return m, z, a, b, c, d, sigma, niter, chisq, npoints, npars
    else:
        return m, z, a, b, c, d, niter, chisq, npoints, npars
