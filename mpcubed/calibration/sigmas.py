#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numba import jit    


@jit(nopython=True)
def _par_sigma_function(idx, res, errsq, err):
    
    weights = 1/(errsq + (err**2)[idx])
    par = np.bincount(idx, weights*res)/np.bincount(idx, weights)
    diff = np.bincount(idx, weights**2*(res - par[idx])**2 - weights)
    
    return par, diff


@jit(nopython=True)
def find_par_sigma(idx, res, errsq, maxiter=10):
    
    # Search for a solution between 0 and 2.
    nvals = np.amax(idx) + 1
    err1 = np.zeros(nvals)
    err2 = 2*np.ones(nvals)
    
    # Compute the value of the function at the beginning the interval.
    par, diff1 = _par_sigma_function(idx, res, errsq, err1)
    args1, = np.where(diff1 < 1e-10)
    
    # Compute the value of the function at the end the interval.
    par, diff2 = _par_sigma_function(idx, res, errsq, err2)
    args2, = np.where(diff2 > 1e-10)
    
    # Find the solution.
    for niter in range(maxiter):
        
        err3 = (err1 + err2)/2.
        par, diff3 = _par_sigma_function(idx, res, errsq, err3)
    
        err1 = np.where(diff3 > 1e-10, err3, err1)
        err2 = np.where(diff3 > 1e-10, err2, err3)
        
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.
    
    par, _ = _par_sigma_function(idx, res, errsq, err3)
    
    return par, err3


@jit(nopython=True)
def _sigma_function(idx, ressq, errsq, err):
    
    weights = 1/(errsq + (err**2)[idx])
    term = ressq*weights**2 - weights
    term = np.bincount(idx, term)
    
    return term


@jit(nopython=True)
def find_sigma(idx, residuals, errsq, maxiter=10):
    
    # Search for a solution between 0 and 2.
    nvals = np.amax(idx) + 1
    err1 = np.zeros(nvals)
    err2 = 2*np.ones(nvals)
    
    ressq = residuals*residuals
    
    # Compute the value of the function at the beginning the interval.
    diff1 = _sigma_function(idx, ressq, errsq, err1)
    args1, = np.where(diff1 < 1e-10)
    
    # Compute the value of the function at the end the interval.
    diff2 = _sigma_function(idx, ressq, errsq, err2)
    args2, = np.where(diff2 > 1e-10)

    # Find the solution.
    for niter in range(maxiter):
        
        err3 = (err2 + err1)/2.
        diff3 = _sigma_function(idx, ressq, errsq, err3)
        
        err1 = np.where(diff3 > 1e-10, err3, err1)
        err2 = np.where(diff3 > 1e-10, err2, err3)
    
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.

    return err3
