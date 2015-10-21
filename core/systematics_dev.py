#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from coarse_decor import find_sigma

def coarse_positions(ind1, ind2, mag, emag, ind3=None, x=None, y=None, use_weights=True, maxiter=100, eps=1e-3, verbose=False):
    """ 
    Given data this code will fit a model of the form:
    mag = m[ind1] + z[ind2] + a[ind3]*np.sin(2*np.pi*y) + b[ind3]*np.cos(2*np.pi*y)
    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    ipx_x = False
    ipx_y = False
    
    if ind3 is not None:
        npars3 = 2*len(np.unique(ind3))
        
        if x is not None:
            ipx_x = True
            npars = npars + npars3
            
        if y is not None:
            ipx_y = True
            npars = npars + npars3
        
        if not ipx_x and not ipx_y:
            print 'Warning: intrapixel indices set but no positions given.'
            exit()
    
    # Create the necessary arrays.
    weights = 1./emag**2
    z = np.zeros(np.amax(ind2)+1)
    
    if ipx_y:
        sny = np.sin(2*np.pi*y)
        csy = np.cos(2*np.pi*y)
        b = np.zeros(np.amax(ind3)+1)
    sol2 = 0.
    
    if ipx_x:
        snx = np.sin(2*np.pi*x)
        csx = np.cos(2*np.pi*x)
        d = np.zeros(np.amax(ind3)+1)
    sol3 = 0.
    
    # Start the iterative solution.
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Computation of parameters.
        m = np.bincount(ind1, weights*(mag-z[ind2]-sol2-sol3))/np.bincount(ind1, weights)
        z = np.bincount(ind2, weights*(mag-m[ind1]-sol2-sol3))/np.bincount(ind2, weights)
        
        sol1 = m[ind1] + z[ind2]
        
        if ipx_y:
            a = np.bincount(ind3, weights*(mag-sol1-b[ind3]*csy-sol3)*sny)/np.bincount(ind3, weights*sny**2)
            b = np.bincount(ind3, weights*(mag-sol1-a[ind3]*sny-sol3)*csy)/np.bincount(ind3, weights*csy**2)
            sol2 = a[ind3]*sny + b[ind3]*csy
        
        if ipx_x:
            c = np.bincount(ind3, weights*(mag-sol1-sol2-d[ind3]*csx)*snx)/np.bincount(ind3, weights*snx**2)
            d = np.bincount(ind3, weights*(mag-sol1-sol2-c[ind3]*snx)*csx)/np.bincount(ind3, weights*csx**2)
            sol3 = c[ind3]*snx + d[ind3]*csx
        
        if use_weights:
            res = mag - sol1 - sol2 - sol3
            sigma = find_sigma(ind1, res, emag)
            weights = 1./(emag**2+(sigma**2)[ind1])
        
        if (niter > 0):
            
            crit1 = np.nanmax(np.abs(sol1-sol1_old))
            crit2 = np.nanmax(np.abs(sol2-sol2_old))
            crit3 = np.nanmax(np.abs(sol3-sol3_old))
            
            if verbose:
                print ' crit1 = %g, crit2 = %g, crit3 = %g'%(crit1, crit2, crit3)
            
            if (crit1 < eps) & (crit2 < eps) & (crit3 < eps):
                break
        
        sol1_old = np.copy(sol1)
        sol2_old = np.copy(sol2)
        sol3_old = np.copy(sol3)
    
    # Compute the chi-square value of the solution.
    chi_tmp = weights*(mag - sol1 - sol2 - sol3)**2
    chisq = np.sum(chi_tmp)/(npoints-npars)
    
    return m, z, niter, chisq, npoints, npars
    
