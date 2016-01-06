#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import numexpr as ne
from core import coarse_decor

import matplotlib.pyplot as plt

def sigma_function(idx, ressq, errsq, err):
    
    weights = 1/(errsq + (err**2)[idx])
    term = ressq*weights**2 - weights
    term = np.bincount(idx, term)
    
    return term

def find_sigma(idx, residuals, error, maxiter=10):
    
    # Search for a solution between 0 and 2.
    N = np.amax(idx) + 1
    err1 = np.zeros(N)
    err2 = np.full(N, 2)
    
    ressq = residuals*residuals
    errsq = error*error
    
    # Compute the value of the fucntion at the beginning the interval.
    diff1 = sigma_function(idx, ressq, errsq, err1)
    args1, = np.where(diff1 < 1e-10)
    
    # Compute the value of the fucntion at the end the interval.
    diff2 = sigma_function(idx, ressq, errsq, err2)
    args2, = np.where(diff2 > 1e-10)
    
    # Find the solution.
    for niter in range(maxiter):
        
        err3 = (err2 + err1)/2.
        
        diff3 = sigma_function(idx, ressq, errsq, err3)
        
        err1 = np.where(diff3 > 1e-10, err3, err1)
        err2 = np.where(diff3 > 1e-10, err2, err3)
    
    err3 = (err2 + err1)/2.
    err3[args1] = 0.
    err3[args2] = 2.
            
    return err3
    
def sigma_function2(ressq, errsq, err):
    
    weights = 1/(errsq + err**2)
    term = ressq*weights**2 - weights
    term = np.sum(term)
    
    return term
    
def find_sigma2(idx, residuals, error, maxiter=10, eps=1e-3):
    
    # Search for a solution between 0 and 2.
    N = np.amax(idx) + 1
    err1 = np.zeros(N)
    err2 = 2.*np.ones(N)
    err3 = np.zeros(N)
    
    ressq = residuals*residuals
    errsq = error*error
    
    for i in range(N):
    
        here = (idx == i)
    
        # Compute the value of the fucntion at the beginning the interval.
        diff1 = sigma_function2(ressq[here], errsq[here], err1[i])
        if (diff1 < 1e-10):
            err3[i] = err1
            continue
        
        # Compute the value of the fucntion at the end the interval.
        diff2 = sigma_function2(ressq[here], errsq[here], err2[i])
        if (diff2 > 1e-10):
            err3[i] = err2
            continue
        
        tmpa, tmpb = ressq[here], errsq[here]
        # Find the solution.
        for niter in range(maxiter):
            
            err3[i] = (err2[i] + err1[i])/2.
            
            diff3 = sigma_function2(tmpa, tmpb, err3[i])
            
            if (diff3 > 1e-10):
                err1[i] = err3[i]
            else:
                err2[i] = err3[i]

    return err3
    
def call_fs():
    
    idx = np.repeat(np.arange(100), 50)
    residuals = np.random.randn(5000)
    error = .2*np.ones(5000)
    
    for i in range(1000):
        sigma1 = find_sigma(idx, residuals, error)
        sigma3 = find_sigma2(idx, residuals, error)
        sigma2 = coarse_decor.find_sigma(idx, residuals, error)
        
        print np.allclose(sigma1, sigma2)
        
    return
