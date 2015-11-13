#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import matplotlib.pyplot as plt

from time import time

def fit_sine_iterative(x, y, weights):
    
    sn = np.sin(x)
    cs = np.cos(x)
    
    b = 0.
    for niter in range(50):
        
        a = np.sum(weights*(y - b*cs)*sn)/np.sum(weights*sn**2)
        b = np.sum(weights*(y - a*sn)*cs)/np.sum(weights*cs**2)
        
        if niter > 0:
            crit1 = np.abs(a - a_old)
            crit2 = np.abs(b - b_old)
            
            if (crit1 < 1e-3) & (crit2 < 1e-3):
                print niter
                break
        
        a_old = np.copy(a)
        b_old = np.copy(b)
    
    return a, b
    
    
def fit_sine_exact(x, y, weights):
    
    sn = np.sin(x)
    cs = np.cos(x)
    
    b = np.zeros((2, len(x)))
    b[0] = sn
    b[1] = cs

    b = np.sum(weights*y*b, axis=1)

    A = np.zeros((2, 2, len(x)))
    A[0,0] = sn**2
    A[1,1] = cs**2
    A[1,0] = sn*cs
    A[0,1] = sn*cs
    
    A = np.sum(weights*A, axis=2)
    
    a, b = np.linalg.solve(A, b)
    
    return a, b

def fit_ipx_exact(ind, x, y, mag, emag):
    
    snx = np.sin(2*np.pi*x)
    csx = np.cos(2*np.pi*x)
    sny = np.sin(2*np.pi*y)
    csy = np.cos(2*np.pi*y)
    weights = 1/emag**2

    b = np.zeros((len(mag), 4))
    b[:,0] = snx
    b[:,1] = csx
    b[:,2] = sny
    b[:,3] = csy
    
    bnew = np.zeros((np.amax(ind)+1, 4))
    for i in range(4):
        bnew[:,i] = np.bincount(ind, weights*mag*b[:,i])
    
    #b = np.sum(weights*mag*b, axis=1)
    
    A = np.zeros((len(mag), 4, 4))
    A[:,0,0] = snx**2
    A[:,1,1] = csx**2
    A[:,2,2] = sny**2
    A[:,3,3] = csy**2
    
    A[:,0,1] = A[:,1,0] = snx*csx
    A[:,0,2] = A[:,2,0] = snx*sny
    A[:,0,3] = A[:,3,0] = snx*csy
    A[:,1,2] = A[:,2,1] = csx*sny
    A[:,1,3] = A[:,3,1] = csx*csy
    A[:,2,3] = A[:,3,2] = sny*csy
    
    Anew = np.zeros((np.amax(ind)+1, 4, 4))
    for i in range(4):
        for j in range(4):
            Anew[:,i,j] = np.bincount(ind, weights*A[:,i,j])
    
    args1, = np.where(np.linalg.det(Anew) < 1e-3)
    args2, = np.where(np.linalg.det(Anew) >= 1e-3)

    #A = np.sum(weights*A, axis=2)
    
    sol = np.zeros((np.amax(ind)+1, 4))
    if len(args2) > 0:
        sol[args2] = np.linalg.solve(Anew[args2], bnew[args2])
    if len(args1) > 0:
        for i in args1:
            sol[i] = np.linalg.lstsq(np.squeeze(Anew[i]), np.squeeze(bnew[i]))[0]
    
    return sol[:,0], sol[:,1], sol[:,2], sol[:,3]

#x = np.linspace(-50*np.pi, 48*np.pi, 50)
#y = .3*np.sin(x + 2.34)
#weights = np.random.rand(50)

#plt.plot(x, y)
#plt.show()

#start = time()
#a, b = fit_sine_exact(x, y, weights)
#print time()-start

#print a, b
#print np.sqrt(a**2 + b**2), np.arctan2(b, a)

#start = time()
#a, b = fit_sine_iterative(x, y, weights)
#print time()-start

#print a, b
#print np.sqrt(a**2 + b**2), np.arctan2(b, a)

ind = np.repeat(np.arange(50), 40)

x = np.linspace(0, 500, 2000)
y = np.linspace(0, 300, 2000)

mag = .82*np.sin(2*np.pi*x) + .45*np.cos(2*np.pi*x) + .39*np.sin(2*np.pi*y) + .67*np.cos(2*np.pi*y)
mag = mag + .2*np.random.randn(2000)
emag = .2*np.ones(2000)

a, b, c, d = fit_ipx_exact(ind, x, y, mag, emag)
print a, b, c, d
plt.errorbar(x, mag, yerr=emag, fmt='.')
plt.plot(x, a[ind]*np.sin(2*np.pi*x) + b[ind]*np.cos(2*np.pi*x) + c[ind]*np.sin(2*np.pi*y) + d[ind]*np.cos(2*np.pi*y))
plt.show()


