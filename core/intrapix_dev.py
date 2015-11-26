#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from scipy.linalg import solve, lstsq

import matplotlib.pyplot as plt

from time import time

@profile
def fit_ipx_exact(ind, x, y, values, weights):
    
    snx = np.sin(x)
    csx = np.cos(x)
    
    sny = np.sin(y)
    csy = np.cos(y)
    
    b = np.zeros((4, len(x)))
    b[0] = snx
    b[1] = csx
    b[2] = sny
    b[3] = csy
    
    A = np.einsum('ik,jk->ijk', b, b)
    
    Adev = np.zeros((4, 4, np.amax(ind)+1))
    for i in range(4):
        for j in range(4):
            Adev[i,j] = np.bincount(ind, A[i,j])
    
    for i in range(500):
        
        bdev = np.zeros((4, np.amax(ind)+1))
        for j in range(4):
            bdev[j] = np.bincount(ind, weights*values*b[j])
            
        for j in range(Adev.shape[2]):
            np.linalg.lstsq(Adev[:,:,j], bdev[:,j])
    
    return
    
@profile
def fit_ipx_iterative(ind, x, y, values, weights):
    
    snx = np.sin(x)
    csx = np.cos(x)
    
    sny = np.sin(y)
    csy = np.cos(y)
    
    b = np.zeros(np.amax(ind)+1)
    c = np.zeros(np.amax(ind)+1)
    d = np.zeros(np.amax(ind)+1)
    
    for niter in range(500):
        
        a = np.bincount(ind, weights*(values - b[ind]*csx - c[ind]*sny - d[ind]*csy)*snx)/np.bincount(ind, weights*snx**2)
        b = np.bincount(ind, weights*(values - a[ind]*snx - c[ind]*sny - d[ind]*csy)*csx)/np.bincount(ind, weights*csx**2)
        c = np.bincount(ind, weights*(values - a[ind]*snx - b[ind]*csx - d[ind]*csy)*sny)/np.bincount(ind, weights*sny**2)
        d = np.bincount(ind, weights*(values - a[ind]*snx - b[ind]*csx - c[ind]*sny)*csy)/np.bincount(ind, weights*csy**2)
    
        sol = np.array([a, b, c, d])
    
        #if niter > 0:
            #crit = np.abs(sol - sol_old)
            
            #if np.all(crit < 1e-3):
                #print niter
                #break
    
        sol_old = np.copy(sol)
        
    return
    
x = np.linspace(0, 400, 5000)
y = np.linspace(0, 300, 5000)
values = .1*np.sin(x) + .2*np.cos(x) + .3*np.sin(y) + .4*np.cos(y)
weights = np.random.rand(5000)
ind = np.repeat(np.arange(100), 50)

print ind

fit_ipx_exact(ind, x, y, values, weights)
fit_ipx_iterative(ind, x, y, values, weights)



