#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def fourier_mat(x, freq, nf):
    
    nx = len(x)
    
    sin_mat = np.zeros((nx, nf))
    cos_mat = np.zeros((nx, nf))
    for i in range(nf):
        sin_mat[:,i] = np.sin(2*np.pi*(i+1)*freq*x)
        cos_mat[:,i] = np.cos(2*np.pi*(i+1)*freq*x)

    return np.hstack([sin_mat, cos_mat])

def fourier_fit(x, y, freq, nf, weights):
    
    fmat = fourier_mat(x, freq, nf)
    fmat = np.hstack([np.ones((len(x),1)), fmat])
    
    pars = np.linalg.lstsq(fmat*np.sqrt(weights[:,None]), y*np.sqrt(weights))[0]
    fit = np.dot(fmat, pars)
    
    return pars, fit
    
def lst_trend(idx, x, y, freq, nf, weights, maxiter=10):
    
    fmat = fourier_mat(x, freq, nf)
    fit2 = np.zeros(len(x))
    
    for niter in range(maxiter):
        
        pars1 = np.bincount(idx, weights*(y - fit2))/np.bincount(idx, weights)
        fit1 = pars1[idx]
        
        pars2 = np.linalg.lstsq(fmat*np.sqrt(weights[:,None]), (y - fit1)*np.sqrt(weights))[0]
        fit2 = np.dot(fmat, pars2)
        
    return pars1, pars2, fit1 + fit2

if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    x = np.linspace(0, 720, 500)
    freq = 1/360.
    nf = 5
    np.random.seed(19910909)
    pars0 = np.random.rand(2*nf)
    print pars0
    print 
    
    fmat = fourier_mat(x, freq, nf)
    y = np.dot(fmat, pars0)
    sigma = .2*np.ones(500)
    y = y + sigma*np.random.randn(500)
    weights = 1/sigma**2
    
    pars, fit = fourier_fit(x, y, freq, nf, weights)
    print pars[1:]
    print
    
    plt.errorbar(x, y, yerr=sigma, fmt='o')
    plt.plot(x, fit)
    plt.show()

    idx = np.repeat(np.arange(2), 250)
    offsets0 = 5*np.random.rand(2)
    print pars0
    print offsets0
    print
    
    y = y + offsets0[idx]
    
    offsets, pars, fit = lst_trend(idx, x, y, freq, nf, weights)
    print pars
    print offsets
    
    plt.errorbar(x, y, yerr=sigma, fmt='o')
    plt.plot(x, fit)
    plt.show()
    
    
    
