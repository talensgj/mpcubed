#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def sysrem(ind1, ind2, values, errors, a2=None, maxiter=50, eps=1e-3):
    
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    if a2 is None:
        a2 = np.ones(np.amax(ind2)+1)
    
    s = values/errors**2
    r = 1./errors**2
    
    niter = 0
    end_crit = True
    
    while (niter < maxiter) & (end_crit):
        
        print niter
        a1 = np.bincount(ind1, s*a2[ind2])/np.bincount(ind1, r*(a2**2)[ind2])
        a2 = np.bincount(ind2, s*a1[ind1])/np.bincount(ind2, r*(a1**2)[ind1])
        
        #sol = a1[ind1]*a2[ind2]
        
        if niter == 0:
            end_crit = True
        else:
            
            crit1 = np.nanmax(np.abs((a1o-a1)/a1o))
            crit2 = np.nanmax(np.abs((a2o-a2)/a2o))
            #crit3 = np.nanmax(np.abs((solo-sol)/solo))
            print '%.4f %.4f'%(crit1, crit2)
            
            end_crit = (crit1 > eps) | (crit2 > eps)
        
        #solo = np.copy(sol)
        a1o = np.copy(a1)
        a2o = np.copy(a2)
            
        niter += 1
    
    chi_tmp = r*(values-a1[ind1]*a2[ind2])**2
    
    chisq = np.sum(chi_tmp)/(npoints-npars)
    chisq_pbin1 = np.bincount(ind1, chi_tmp)
    chisq_pbin2 = np.bincount(ind2, chi_tmp)
    
    return a1, a2, niter, chisq, chisq_pbin1, chisq_pbin2, npoints, npars

def smooth(ind1, ind2, values, errors, a2=None, maxiter=50, eps=1e-3):
    
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    if a2 is None:
        a2 = np.ones(np.amax(ind2)+1)
    
    s = values/errors**2
    r = 1./errors**2
    
    niter = 0
    end_crit = True
    
    while (niter < maxiter) & (end_crit):
        
        print niter
        a1 = np.bincount(ind1, s*a2[ind2])/np.bincount(ind1, r*(a2**2)[ind2])
        a2 = np.bincount(ind2, s*a1[ind1])/np.bincount(ind2, r*(a1**2)[ind1])
        
        if niter == 0:
            end_crit = True
        else:
            
            crit1 = np.nanmax(np.abs((a1o-a1)/a1o))
            crit2 = np.nanmax(np.abs((a2o-a2)/a2o))
            print crit1, crit2
            end_crit = (crit1 > eps) | (crit2 > eps)
        
        a1o = np.copy(a1)
        a2o = np.copy(a2)
            
        niter += 1
    
    chi_tmp = r*(values-a1[ind1]*a2[ind2])**2
    
    chisq = np.sum(chi_tmp)/(npoints-npars)
    chisq_pbin1 = np.bincount(ind1, chi_tmp)
    chisq_pbin2 = np.bincount(ind2, chi_tmp)
    
    return a1, a2, niter, chisq, chisq_pbin1, chisq_pbin2, npoints, npars
