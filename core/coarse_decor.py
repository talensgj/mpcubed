#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
    
    
def find_sigma(ind1, res, err, maxiter=50, eps=1e-3):
    
    err1 = np.zeros(np.amax(ind1)+1)
    err2 = 2*np.ones(np.amax(ind1)+1)
    
    wgt = 1./(err*err + (err1[ind1])**2)
    term2 = np.bincount(ind1, res*res*wgt*wgt)
    term1 = np.bincount(ind1, wgt)
    diff1 = term2 - term1
    
    args1, = np.where(diff1 < 1e-10)
    
    wgt = 1./(err*err + (err2[ind1])**2)
    term2 = np.bincount(ind1, res*res*wgt*wgt)
    term1 = np.bincount(ind1, wgt)
    diff2 = term2 - term1
    
    args2, = np.where(diff2 > 1e-10)
    
    for niter in range(maxiter):
        err3 = (err2 + err1)/2.
        
        wgt = 1./(err*err + (err3[ind1])**2)
        term2 = np.bincount(ind1, res*res*wgt*wgt)
        term1 = np.bincount(ind1, wgt)
        diff3 = term2 - term1
        
        here = (diff3 > 1e-10 )
        err1[here] = err3[here]
        diff1[here] = diff3[here]
        err2[~here] = err3[~here]
        diff2[~here] = diff3[~here]
        
        if np.all((err2-err1) < eps):
            break
    
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.
            
    return err3
    
    
def coarse_decorrelation(ind1, ind2, values, errors, maxiter=100, eps=1e-3, verbose=False, use_weights=True):
    """ 
    Given data this code will fit a model of the form:
    values = m[ind1] + z[ind2]
    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    # Create the necessary arrays.
    weights = 1./(errors**2)
    z = np.zeros(np.amax(ind2)+1)
    sigma2 = np.zeros(np.amax(ind2)+1)
    
    # Start the iterative solution.
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Computation of parameters.
        m = np.bincount(ind1, weights*(values-z[ind2]))/np.bincount(ind1, weights)
        z = np.bincount(ind2, weights*(values-m[ind1]))/np.bincount(ind2, weights)
        
        if use_weights:
            res = values - m[ind1] - z[ind2]
            sigma1 = find_sigma(ind1, res, np.sqrt(errors**2+(sigma2**2)[ind2]))
            sigma2 = find_sigma(ind2, res, np.sqrt(errors**2+(sigma1**2)[ind1]))
            weights = 1./(errors**2+(sigma1**2)[ind1]+(sigma2**2)[ind2])
        
        if (niter > 0):
            
            # Check if the solution has converged.
            critm = np.nanmax(np.abs(m-m_old))
            critz = np.nanmax(np.abs(z-z_old))
            
            if verbose:
                print ' critm = %g, critz = %g'%(critm, critz)
            
            if (critm < eps) & (critz < eps):
                break

        m_old = np.copy(m)
        z_old = np.copy(z)
    
    # Compute the chi-square value of the solution.
    chi_tmp = weights*(values - m[ind1] - z[ind2])**2
    chisq = np.sum(chi_tmp)/(npoints-npars)
    
    if use_weights:
        return m, z, sigma1, sigma2, niter, chisq, npoints, npars
    else:
        return m, z, niter, chisq, npoints, npars


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
    a = np.zeros(np.amax(ind3)+1)
    b = np.zeros(np.amax(ind3)+1)
    c = np.zeros(np.amax(ind3)+1)
    d = np.zeros(np.amax(ind3)+1)
    
    # Start the iterative solution.
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Computation of parameters.
        m = np.bincount(ind1, weights*(mag-z[ind2]-a[ind3]*sny-b[ind3]*csy-c[ind3]*snx-d[ind3]*csx))/np.bincount(ind1, weights)
        z = np.bincount(ind2, weights*(mag-m[ind1]-a[ind3]*sny-b[ind3]*csy-c[ind3]*snx-d[ind3]*csx))/np.bincount(ind2, weights)
        
        a = np.bincount(ind3, weights*(mag-m[ind1]-z[ind2]-b[ind3]*csy-c[ind3]*snx-d[ind3]*csx)*sny)/np.bincount(ind3, weights*sny**2)
        b = np.bincount(ind3, weights*(mag-m[ind1]-z[ind2]-a[ind3]*sny-c[ind3]*snx-d[ind3]*csx)*csy)/np.bincount(ind3, weights*csy**2)
        
        c = np.bincount(ind3, weights*(mag-m[ind1]-z[ind2]-a[ind3]*sny-b[ind3]*csy-d[ind3]*csx)*snx)/np.bincount(ind3, weights*snx**2)
        d = np.bincount(ind3, weights*(mag-m[ind1]-z[ind2]-a[ind3]*sny-b[ind3]*csy-c[ind3]*snx)*csx)/np.bincount(ind3, weights*csx**2)
        
        if use_weights:
            res = mag - m[ind1] - z[ind2] - a[ind3]*sny - b[ind3]*csy - c[ind3]*snx - d[ind3]*csx
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
