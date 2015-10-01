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
    
    
def coarse_decorrelation(ind1, ind2, values, errors, z=None, maxiter=100, eps=1e-3, verbose=False, use_weights=True):
    """ 
    Given data this code will fit a model of the form:
    values = m[ind1] + z[ind2]
    """
    
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    length1 = np.amax(ind1)+1
    length2 = np.amax(ind2)+1
    
    weights = 1./(errors**2)
    
    if z is None:
        z = np.zeros(length2)
        
    sigma1 = np.zeros(length1)
    sigma2 = np.zeros(length2)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
            
        m = np.bincount(ind1, weights*(values-z[ind2]))/np.bincount(ind1, weights)
        z = np.bincount(ind2, weights*(values-m[ind1]))/np.bincount(ind2, weights)
        
        if use_weights:
            res = values - m[ind1] - z[ind2]
            sigma1 = find_sigma(ind1, res, np.sqrt(errors**2+sigma2[ind2]**2))
            sigma2 = find_sigma(ind2, res, np.sqrt(errors**2+sigma1[ind1]**2))
            weights = 1./(errors**2+sigma1[ind1]**2+sigma2[ind2]**2)
            
        sol = m[ind1] + z[ind2]
        
        if (niter > 0):
            
            critsol = np.nanmax(np.abs(sol-sol_old))
            crits1 = np.nanmax(np.abs(sigma1-sigma1_old))
            crits2 = np.nanmax(np.abs(sigma2-sigma2_old))
            
            if verbose:
                print ' critsol = %g, crits1 = %g, crits2 = %g'%(critsol, crits1, crits2)
            
            if (critsol < eps):
                break

        sol_old = np.copy(sol)
        sigma1_old = np.copy(sigma1)
        sigma2_old = np.copy(sigma2)
            
    chi_tmp = weights*(values - m[ind1] - z[ind2])**2
    chisq = np.sum(chi_tmp)/(npoints-npars)
    
    if use_weights:
        return m, z, sigma1, sigma2, niter, chisq, npoints, npars
    else:
        return m, z, niter, chisq, npoints, npars


def coarse_positions(ind1, ind2, ind3, y, mag, emag, maxiter=100, eps=1e-3, verbose=False, use_weights=True):
    """ 
    Given data this code will fit a model of the form:
    mag = m[ind1] + z[ind2] + a[ind3]*np.sin(2*np.pi*y) + b[ind3]*np.cos(2*np.pi*y)
    """
    
    npoints = len(mag)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars3 = 2*len(np.unique(ind3))
    npars = npars1 + npars2 + npars3
    
    length1 = np.amax(ind1)+1
    length2 = np.amax(ind2)+1
    length3 = np.amax(ind3)+1
    
    sn = np.sin(2*np.pi*y)
    cs = np.cos(2*np.pi*y)
    
    weights = 1./emag**2
    
    z = np.zeros(length2)
    a = np.zeros(length3)
    b = np.zeros(length3)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        m = np.bincount(ind1, weights*(mag-z[ind2]-a[ind3]*sn-b[ind3]*cs))/np.bincount(ind1, weights)
        z = np.bincount(ind2, weights*(mag-m[ind1]-a[ind3]*sn-b[ind3]*cs))/np.bincount(ind2, weights)
        
        a = np.bincount(ind3, weights*(mag-m[ind1]-z[ind2]-b[ind3]*cs)*sn)/np.bincount(ind3, weights*sn**2)
        b = np.bincount(ind3, weights*(mag-m[ind1]-z[ind2]-a[ind3]*sn)*cs)/np.bincount(ind3, weights*cs**2)
        
        if use_weights:
            res = mag - m[ind1] - z[ind2] - a[ind3]*sn - b[ind3]*cs
            sigma = find_sigma(ind1, res, emag)
            weights = 1./(emag**2+sigma[ind1]**2)
        
        if (niter > 0):

            critm = np.nanmax(np.abs(m-m_old))
            critz = np.nanmax(np.abs(z-z_old))
            crita = np.nanmax(np.abs(a-a_old))
            critb = np.nanmax(np.abs(b-b_old))
            
            if verbose:
                print ' critm = %g, critz = %g, crita = %g, critb = %g'%(critm, critz, crita, critb)
            
            if (critm < eps) & (critz < eps) & (crita < eps) & (critb < eps):
                break
        
        m_old = np.copy(m)
        z_old = np.copy(z)
        a_old = np.copy(a)
        b_old = np.copy(b)
    
    chi_tmp = weights*(mag - m[ind1] - z[ind2] - a[ind3]*sn - b[ind3]*cs)**2
    chisq = np.sum(chi_tmp)/(npoints-npars)
    
    if use_weights:
        return m, z, a, b, sigma, niter, chisq, npoints, npars
    else:
        return m, z, a, b, niter, chisq, npoints, npars
