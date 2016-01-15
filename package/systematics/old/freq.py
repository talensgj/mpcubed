#!/usr/bin/env python
# -*- coding: utf-8 -*-

def remove_freq(x, y, freq, n, weights = None):
    
    nrows = len(x)
    ncols = 2*n + 1
    
    A = np.zeros((nrows, ncols))
    
    A[:,0] = 1
    for i in range(n):
        A[:,2*i + 1] = np.sin(2*np.pi*(i + 1)*freq*x)
        A[:,2*i + 2] = np.cos(2*np.pi*(i + 1)*freq*x)

    if weights is None:
        pars = np.linalg.lstsq(A, y)[0]
    else:
        Aw = A*np.sqrt(weights[:,np.newaxis])
        yw = y*np.sqrt(weights)
        pars = np.linalg.lstsq(Aw, yw)[0]
        
    fit = np.dot(A, pars)
    
    return pars, fit

def lstcor(idx, x, y, freq, n, weights=None):
    
    nrows = len(x)
    ncols = np.amax(idx) + 1
    
    baseline = np.zeros((nrows, ncols))
    baseline[range(nrows),idx] = 1. 

    sine_terms = np.zeros((nrows, n))
    for i in range(n):
        sine_terms[i] = np.sin(2*np.pi*(i + 1)*freq*x)
        
    cosine_terms = np.zeros((nrows, n))
    for i in range(n):
        cosine_terms[i] = np.cos(2*np.pi*(i + 1)*freq*x)
        
    matrix = np.hstack([baseline, sine_terms, cosine_terms])
    
    pars = np.linalg.lstsq(matrix, y)[0]
    
    fit = np.dot(matrix, pars)
    
    return pars, fit
