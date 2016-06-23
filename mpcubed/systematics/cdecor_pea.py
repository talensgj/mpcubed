#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import linalg

import sigmas

from collections import namedtuple

Quality = namedtuple('Quality', 'niter chisq npoints npars') 

def cdecor(idx1, idx2, idx3, idx4, mag, error, x, y, maxiter=100, dtol=1e-3, verbose=False):
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(mag)
    npars1 = np.amax(idx1) + 1
    npars2 = np.amax(idx2) + 1
    npars3 = np.amax(idx3) + 1
    npars4 = 4*(np.amax(idx4) + 1)
    npars = npars2 + npars3 + npars4
    
    # Create arrays.
    weights = 1./error**2
    
    strides = np.cumsum(np.bincount(idx4))
    strides = np.append(0, strides)
    mat = np.vstack([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T
    ipx = np.zeros(len(mag))

    par2 = np.bincount(idx2, weights*mag)/np.bincount(idx2, weights)
    par3 = np.zeros(npars3)
    par4 = np.zeros((npars4/4, 4))
    
    err1 = np.zeros(npars1)
    err3 = np.zeros(npars3)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        err1 = sigmas.find_sigma(idx1, mag - par2[idx2] - par3[idx3] - ipx, error**2 + (err3**2)[idx3])
        par3, err3 = sigmas.find_par_sigma(idx3, mag - par2[idx2] - ipx, error**2 + (err1**2)[idx1])
        
        weights = 1/(error**2 + (err1**2)[idx1] + (err3**2)[idx3])
        
        par2 = np.bincount(idx2, weights*(mag - par3[idx3] - ipx))/np.bincount(idx2, weights)
        par3 = np.bincount(idx3, weights*(mag - par2[idx2] - ipx))/np.bincount(idx3, weights)
        
        res = mag - par2[idx2] - par3[idx3]
        wsqrt = np.sqrt(weights)
        for i in range(npars4/4):
            par4[i] = linalg.lstsq(mat[strides[i]:strides[i+1],:]*wsqrt[strides[i]:strides[i+1]:,None], res[strides[i]:strides[i+1]]*wsqrt[strides[i]:strides[i+1]])[0]
            ipx[strides[i]:strides[i+1]] = np.dot(mat[strides[i]:strides[i+1],:], par4[i])
        
        # Check if the solution has converged.
        if (niter > 0):
            
            tmp2 = np.abs(par2 - par2_old)
            tmp3 = np.abs(par3 - par3_old)
            tmp4 = np.abs(par4 - par4_old) 
            
            dcrit2 = np.nanmax(tmp2)
            dcrit3 = np.nanmax(tmp3)
            dcrit4 = np.nanmax(tmp4)
            
            if verbose:
                print '{}/{}, {:.3f}'.format(np.sum(tmp2<dtol), npars2, dcrit2)
                print '{}/{}, {:.3f}'.format(np.sum(tmp3<dtol), npars3, dcrit3)
                print '{}/{}, {:.3f}'.format(np.sum(tmp4<dtol), npars4, dcrit4)
            
            if (dcrit2 < dtol) & (dcrit3 < dtol) & (dcrit4 < dtol):
                print 'Solution has converged, ending the iterations.'
                break
                
        if (niter > 1):
            
            tmp2 = np.abs(par2 - par2_older)
            tmp3 = np.abs(par3 - par3_older)
            tmp4 = np.abs(par4 - par4_older) 
            
            dcrit2 = np.nanmax(tmp2)
            dcrit3 = np.nanmax(tmp3)
            dcrit4 = np.nanmax(tmp4)
            
            if (dcrit2 < dtol) & (dcrit3 < dtol) & (dcrit4 < dtol):
                print 'Solution is oscillating, ending the iterations.'
                break
        
        if (niter > 0):
            par2_older = np.copy(par2_old)
            par3_older = np.copy(par3_old)
            par4_older = np.copy(par4_old)
        
        par2_old = np.copy(par2)
        par3_old = np.copy(par3)
        par4_old = np.copy(par4)
    
    # Compute the chi-square of the fit.
    chisq = weights*(mag - par2[idx2] - par3[idx3] - ipx)**2        
    chisq = np.sum(chisq)
    
    return par2, par3, par4, err1, err3, Quality(niter, chisq, npoints, npars)
