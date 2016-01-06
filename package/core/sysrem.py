#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def sysrem(ind1, ind2, values, errors, maxiter=100, eps=1e-3, verbose=False):
    """ 
    Given data this code will fit a model of the form:
    values = a1[ind1]*a2[ind2]
    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    # Create the necessary arrays.
    weights = 1./errors**2
    a2 = np.ones(np.amax(ind2)+1)
    
    # Start the iterative solution.
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Computation of parameters.
        a1 = np.bincount(ind1, weights*values*a2[ind2])/np.bincount(ind1, weights*(a2**2)[ind2])
        a2 = np.bincount(ind2, weights*values*a1[ind1])/np.bincount(ind2, weights*(a1**2)[ind1])
        
        if (niter > 0):
            
            # Check if the solution has converged.
            crit1 = np.nanmax(np.abs((a1_old-a1)/a1_old))
            crit2 = np.nanmax(np.abs((a2_old-a2)/a2_old))
            
            if verbose:
                print ' crit1 = %g, crit2 = %g'%(crit1, crit2)
            
            if (crit1 < eps) & (crit2 < eps):
                break

        a1_old = np.copy(a1)
        a2_old = np.copy(a2)
    
    # Compute the chi-square value of the solution.
    chi_tmp = weights*(values-a1[ind1]*a2[ind2])**2
    chisq = np.sum(chi_tmp)/(npoints-npars)
    
    return a1, a2, niter, chisq, npoints, npars


def intrarem(ind1, ind2, ind3, y, flux, eflux, maxiter=100, eps=1e-3, verbose=False):
    """
    Given data this code will fit a model of the form:
    flux = F[ind1]*T[ind2]*(a[ind3]*np.sin(2*np.pi*y)+b[ind3]*np.cos(2*np.pi*y)+1)
    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(flux)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars3 = 2*len(np.unique(ind3))
    npars = npars1 + npars2 + npars3
    
    # Create the necessary arrays.
    weights = 1./eflux**2
    cs = np.cos(2.*np.pi*y)
    sn = np.sin(2.*np.pi*y)
    T = np.ones(np.amax(ind2)+1)
    sol2 = np.ones(len(flux))
    b = np.zeros(np.amax(ind3)+1)
 
    # Start the iterative solution.
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = %i'%niter
        
        # Computation of parameters.
        F = np.bincount(ind1, weights*flux*sol2*T[ind2])/np.bincount(ind1, weights*(T**2)[ind2]*(sol2)**2)
        T = np.bincount(ind2, weights*flux*sol2*F[ind1])/np.bincount(ind2, weights*(F**2)[ind1]*(sol2)**2)
        
        sol1 = F[ind1]*T[ind2]
        
        a = np.bincount(ind3, weights*(flux-sol1*(b[ind3]*cs+1))*sol1*sn)/np.bincount(ind3, weights*(sol1*sn)**2)
        b = np.bincount(ind3, weights*(flux-sol1*(a[ind3]*sn+1))*sol1*cs)/np.bincount(ind3, weights*(sol1*cs)**2)
    
        sol2 = a[ind3]*sn + b[ind3]*cs + 1
    
        if (niter > 0):

            # Check if the solution has converged.
            critF = np.nanmax(np.abs((F-F_old)/F_old))
            critT = np.nanmax(np.abs((T-T_old)/T_old))
            crita = np.nanmax(np.abs((a-a_old)/a_old))
            critb = np.nanmax(np.abs((b-b_old)/b_old))
            
            if verbose:
                print ' critF = %g, critT = %g, crita = %g, critb = %g'%(critF, critT, crita, critb)
            
            if (critF < eps) & (critT < eps) & (crita < eps) & (critb < eps):
                break
            
        F_old = np.copy(F)
        T_old = np.copy(T)
        a_old = np.copy(a)
        b_old = np.copy(b)
    
    # Compute the chi-square value of the solution.
    chi_tmp = weights*(flux - F[ind1]*T[ind2]*(a[ind3]*sn+b[ind3]*cs+1))**2
    chisq = np.sum(chi_tmp)/(npoints-npars)
            
    return F, T, a, b, niter, chisq, npoints, npars