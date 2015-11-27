#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from coarse_decor import find_sigma

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
    
    return par1, par2, niter, chisq, npoints, npars

@profile
def spatial_decor(ind1, ind2, ind3, mag, emag, x, y, maxiter=100, eps=1e-3, verbose=True):
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars3 = 4*len(np.unique(ind3))
    npars = npars1 + npars2 + npars3
    
    # Create the necessary arrays.
    weights = 1/emag**2
    z = np.zeros(np.amax(ind2) + 1)
    
    snx = np.sin(2*np.pi*x)
    csx = np.cos(2*np.pi*x)
    sny = np.sin(2*np.pi*y)
    csy = np.cos(2*np.pi*y)
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
    
    return m, z, A, niter, chisq, npoints, npars


@profile
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


