#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gc

import numpy as np
from scipy import linalg

import sigmas

from collections import namedtuple

Quality = namedtuple('Quality', 'niter chisq npoints npars') 

def cdecor(idx2, value, error, maxiter=100, dtol=1e-3, verbose=True):
    """ Perform a coarse decorrelation.

    Args:
        idx1 (int): Indices along which to calculate the first parameter.
        idx2 (int): Indices along which to calculate the second parameter.
        value (float): Values to fit.
        error (float): Measurement errors corresponding to the values.
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.

    Returns:
        par1 (float): The parameters corresponding to idx1.
        par2 (float): The parameters corresponding to idx2.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.

    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars2 = np.amax(idx2) + 1
    npars = npars2
    
    # Create arrays.
    weights = 1./error**2
    par2 = np.zeros(npars2)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        par2 = np.bincount(idx2, weights*value)/np.bincount(idx2, weights)
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            
            if (dcrit2 < dtol):
                break
        
        par2_old = np.copy(par2)
    
    # Compute the chi-square of the fit.
    chisq = weights*(value - par2[idx2])**2        
    chisq = np.sum(chisq)
    
    return par2, Quality(niter, chisq, npoints, npars)

def cdecor_intrapix(idx2, idx3, value, error, x, y, maxiter=100, dtol=1e-3, verbose=True):
    """ Perform a coarse decorrelation with intrapixel variations.
    
    Args:
        idx1 (int): Indices along which to calculate the first parameter.
        idx2 (int): Indices along which to calculate the second parameter.
        idx3 (int): Indices along which to calculate the intrapixel variations.
        value (float): Values to fit.
        error (float): Measurement errors corresponding to the values.
        x (float): The x position corresponding to the values.
        y (float): The y position corresponding to the values. 
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.
        
    Returns:
        par1 (float): The parameters corresponding to idx1.
        par2 (float): The parameters corresponding to idx2.
        par3 (float): The amplitudes of the intrapixel variations corresponding
            to idx3.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.
    
    """
    
    sort = np.argsort(idx3)
    idx2 = idx2[sort]
    idx3 = idx3[sort]
    value = value[sort]
    error = error[sort]
    x = x[sort]
    y = y[sort]
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars2 = np.amax(idx2) + 1
    npars3 = 4*(np.amax(idx3) + 1)
    npars = npars2 + npars3
    
    # Create arrays.
    weights = 1./error**2
    par2 = np.zeros(npars2)
    par3 = np.zeros((npars3/4, 4))
    
    strides = np.cumsum(np.bincount(idx3))
    strides = np.append(0, strides)
    mat = np.vstack([np.sin(2*np.pi*x), np.cos(2*np.pi*x), np.sin(2*np.pi*y), np.cos(2*np.pi*y)]).T
    ipx = np.zeros(len(value))
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
        
        # Compute the parameters.
        par2 = np.bincount(idx2, weights*(value - ipx))/np.bincount(idx2, weights)
        
        res = value - par2[idx2]
        wsqrt = np.sqrt(weights)
        for i in range(npars3/4):
            par3[i] = np.linalg.lstsq(mat[strides[i]:strides[i+1],:]*wsqrt[strides[i]:strides[i+1]:,None], res[strides[i]:strides[i+1]]*wsqrt[strides[i]:strides[i+1]])[0]
            #ipx[strides[i]:strides[i+1]] = np.dot(mat[strides[i]:strides[i+1],:], par3[i])
            ipx[strides[i]:strides[i+1]] = np.sum(mat[strides[i]:strides[i+1],:]*par3[i], axis=1)
        
        gc.collect()        
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            dcrit3 = np.nanmax(np.abs(par3 - par3_old))
            
            if (dcrit2 < dtol) & (dcrit3 < dtol):
                break
        
        par2_old = np.copy(par2)
        par3_old = np.copy(par3)
    
    # Compute the chi-square of the fit.
    chisq = weights*(value - par2[idx2] - ipx)**2        
    chisq = np.sum(chisq)
    
    return par2, par3, Quality(niter, chisq, npoints, npars)
    
def cdecor_sigmas(idx1, idx2, value, error, sigma1, sigma2, maxiter=100, dtol=1e-3, verbose=True):
    """ Perform a coarse decorrelation with extra error terms.
    
    Args:
        idx1 (int): Indices along which to calculate the first parameter.
        idx2 (int): Indices along which to calculate the second parameter.
        error (float): Measurement errors corresponding to the values.
        sigma1 (float): Initial value for the extra error corresponding to
            idx1.
        sigma2 (float): Initial value for the extra error corresponding to
            idx2.
        maxiter (int): The maximum number of iterations to perform. Defualt
            is 100.
        dtol (float): Maximum allowed change in the parameters, iteration
            terminates if the cahnge falls below this value. Default is 1e-3.
        verbose (bool): Output the current iteration. Default is True.
        
    Returns:
        par1 (float): The parameters corresponding to idx1.
        par2 (float): The parameters corresponding to idx2.
        sigma1 (float): The extra error corresponding to idx1.
        sigma2 (float): The extra error corresponding to idx2.
        quality: A named tuple with fields niter, chisq, npoints and npars
            describing the number of iterations, the chi-square value, the
            number of datapoints and the number of parameters of the fit.
    
    """
    
    # Determine the number of datapoints and parameters to fit.
    npoints = len(value)
    npars1 = np.amax(idx1) + 1
    npars2 = np.amax(idx2) + 1
    npars = npars2
    
    # Create arrays.
    par1 = np.zeros(npars1)
    par2 = np.zeros(npars2)
    
    for niter in range(maxiter):
        
        if verbose:
            print 'niter = {}'.format(niter)
            
        # Compute the parameters.
        sigma1 = sigmas.find_sigma(idx1, value - par2[idx2], error**2 + (sigma2**2)[idx2])
        par2, sigma2 = sigmas.find_par_sigma(idx2, value, error**2 + (sigma1**2)[idx1])
        
        # Check if the solution has converged.
        if (niter > 0):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_old))
            
            if (dcrit2 < dtol):
                break
        
        # Check if the solution is oscillating?
        if (niter > 1):
            
            dcrit2 = np.nanmax(np.abs(par2 - par2_older))
            
            if (dcrit2 < dtol):
                break
        
        if (niter > 0):
            par2_older = np.copy(par2_old)
        
        par2_old = np.copy(par2)
        
    # Compute the chi-square of the fit.
    chisq = (value - par2[idx2])**2/(error**2 + (sigma1**2)[idx1] + (sigma2**2)[idx2])     
    chisq = np.sum(chisq)
    
    return par2, sigma1, sigma2, Quality(niter, chisq, npoints, npars)
