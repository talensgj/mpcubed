#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def coarse_decor(idx1, idx2, value, error, maxiter=100, dtol=1e-3, verbose=True):
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars1 = np.amax(idx1) + 1
    npars2 = np.amax(idx2) + 1
    npars = npars1 + npars2
    
    # Create arrays.
    weights = 1./error**2
    par2 = np.zeros(npars2)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Compute the parameters.
        par1 = np.bincount(idx1, weights*(value - par2[idx2]))/np.bincount(idx1, weights)
        par2 = np.bincount(idx2, weights*(value - par1[idx1]))/np.bincount(idx2, weights)
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit1 = np.nanmax(np.abs(par1 - par1_old))
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            
            if (dcrit1 < dtol) & (dcrit2 < dtol):
                break
        
        par1_old = np.copy(par1)
        par2_old = np.copy(par2)
    
    # Compute the chi-square of the fit.
    chisq = weights*(value - par1[idx1] - par2[idx2])**2        
    chisq = np.sum(chisq)
    
    return par1, par2, niter, chisq, npoints, npars

def coarse_decor_intrapix(idx1, idx2, idx3, value, error, x, y, maxiter=100, dtol=1e-3, verbose=True):
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars1 = np.amax(idx1) + 1
    npars2 = np.amax(idx2) + 1
    npars3 = 4*(np.amax(idx3) + 1)
    npars = npars1 + npars2 + npars3
    
    # Create arrays.
    weights = 1./error**2
    par2 = np.zeros(npars2)
    
    snx = np.sin(2*np.pi*x)
    csx = np.cos(2*np.pi*x)
    sny = np.sin(2*np.pi*y)
    csy = np.cos(2*np.pi*y)
    
    sol2 = 0
    sol3 = 0
    b = np.zeros(npars3/4)
    d = np.zeros(npars3/4)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Compute the parameters.
        par1 = np.bincount(idx1, weights*(value - par2[idx2] - sol2 - sol3))/np.bincount(idx1, weights)
        par2 = np.bincount(idx2, weights*(value - par1[idx1] - sol2 - sol3))/np.bincount(idx2, weights)
        
        sol1 = par1[idx1] + par2[idx2]
        
        a = np.bincount(idx3, weights*(value - sol1 - b[idx3]*csx - sol3)*snx)/np.bincount(idx3, weights*snx**2)
        b = np.bincount(idx3, weights*(value - sol1 - a[idx3]*snx - sol3)*csx)/np.bincount(idx3, weights*csx**2)
        
        sol2 = a[idx3]*snx + b[idx3]*csx
        
        c = np.bincount(idx3, weights*(value - sol1 - sol2 - d[idx3]*csy)*sny)/np.bincount(idx3, weights*sny**2)
        d = np.bincount(idx3, weights*(value - sol1 - sol2 - c[idx3]*sny)*csy)/np.bincount(idx3, weights*csy**2)
        
        sol3 = c[idx3]*sny + d[idx3]*csy
        
        par3 = np.vstack([a, b, c, d]).T
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit1 = np.nanmax(np.abs(par1 - par1_old))
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            dcrit3 = np.nanmax(np.abs(par3 - par3_old))
            
            if (dcrit1 < dtol) & (dcrit2 < dtol) & (dcrit3 < dtol):
                break
        
        par1_old = np.copy(par1)
        par2_old = np.copy(par2)
        par3_old = np.copy(par3)
    
    # Compute the chi-square of the fit.
    chisq = weights*(value - sol1 - sol2 - sol3)**2        
    chisq = np.sum(chisq)
    
    return par1, par2, par3, niter, chisq, npoints, npars

# UNTESTED
def sigma_function(idx, residuals, weights):
    
    term2 = np.bincount(idx, (residuals**2)*(weights**2))
    term1 = np.bincount(idx, weights)
    
    return term2 - term1

# UNTESTED
def find_sigma(idx, residuals, error, maxiter=10, eps=1e-3):
    
    # Search for a solution between 0 and 2.
    err1 = np.zeros(np.amax(idx) + 1)
    err2 = 2*np.ones(np.amax(idx) + 1)
    
    # Compute the value of the fucntion at the beginning the interval.
    weights = 1/(error**2 + (err1**2)[idx])
    diff1 = sigma_function(idx, residuals, weights)
    args1, = np.where(diff1 < 1e-10)
    
    # Compute the value of the fucntion at the end the interval.
    weights = 1/(error**2 + (err2**2)[idx])
    diff2 = sigma_function(idx, residuals, weights)
    args2, = np.where(diff2 > 1e-10)
    
    # Find the solution.
    for niter in range(maxiter):
        
        err3 = (err2 + err1)/2.
        
        weights = 1/(error**2 + (err3**2)[idx])
        diff3 = sigma_function(idx, residuals, weights)
        
        here = (diff3 > 1e-10 )
        err1[here] = err3[here]
        diff1[here] = diff3[here]
        err2[~here] = err3[~here]
        diff2[~here] = diff3[~here]
        
        if np.all((err2 - err1) < eps):
            break
    
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.
            
    return err3

# UNTESTED
def coarse_decor_sigmas(idx1, idx2, value, error, maxiter=100, dtol=1e-3, verbose=True):
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars1 = np.amax(idx1) + 1
    npars2 = np.amax(idx2) + 1
    npars = npars1 + npars2
    
    # Create arrays.
    weights = 1./error**2
    par2 = np.zeros(npars2)
    sigma2 = np.zeros(npars2)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Compute the parameters.
        par1 = np.bincount(idx1, weights*(value - par2[idx2]))/np.bincount(idx1, weights)
        par2 = np.bincount(idx2, weights*(value - par1[idx1]))/np.bincount(idx2, weights)
        
        sigma1 = find_sigma(idx1, value - par1[idx1] - par2[idx2], np.sqrt(error**2 + (sigma2**2)[idx2]))
        sigma2 = find_sigma(idx2, value - par1[idx1] - par2[idx2], np.sqrt(error**2 + (sigma1**2)[idx1]))
        weights = 1/(error**2 + (sigma1**2)[idx1] + (sigma2**2)[idx2])
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit1 = np.nanmax(np.abs(par1 - par1_old))
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            
            if (dcrit1 < dtol) & (dcrit2 < dtol):
                break
        
        par1_old = np.copy(par1)
        par2_old = np.copy(par2)
    
    # Compute the chi-square of the fit.
    chisq = weights*(value - par1[idx1] - par2[idx2])**2        
    chisq = np.sum(chisq)
    
    return par1, par2, sigma1, sigma2, niter, chisq, npoints, npars
