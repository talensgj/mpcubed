#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

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

if __name__ == '__main__':
    
    x = np.linspace(0, 1, 100)
    y = 7.32 + .23*np.cos(2*np.pi*x) + .45*np.sin(2*np.pi*x)

    print remove_freq(x, y, 1., 1)
