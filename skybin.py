#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from fLCfile import fLCfile
from core.coarse_decor_dev import new_coarse_decor_sigmas, coarse_decor_sigmas
from core.coordinate_grids import HealpixGrid
from usefull_functions_dev import flux2mag

import matplotlib.pyplot as plt
from viridis import viridis

f = fLCfile('/data2/talens/2015Q2/LPE/fLC_201504ALPE.hdf5')

ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])
nobs = nobs.astype('int')

hg = HealpixGrid(8)

skyidx = hg.find_gridpoint(ra, dec)

skyidx, idx = np.unique(skyidx, return_inverse=True)

here = (idx == 3)
staridx, = np.where(here)
staridx = np.repeat(staridx, nobs[here])
flux0, eflux0, lstseq, sky, flags = f.read_data(['flux0', 'eflux0', 'lstseq', 'sky', 'flag'], ascc[here], nobs[here])
lstseq = lstseq.astype('int')

here = (flux0 > 0) & (eflux0 > 0) & (sky > 0) & (flags < 1)
flux0 = flux0[here]
eflux0 = eflux0[here]
lstseq = lstseq[here]
staridx = staridx[here]

mag0, emag0 = flux2mag(flux0, eflux0)

staridx, idx1 = np.unique(staridx, return_inverse=True)
lstseq, idx2 = np.unique(lstseq, return_inverse=True)

sigma1 = np.zeros(np.amax(idx1) + 1)
sigma2 = np.zeros(np.amax(idx2) + 1)

m, s0, sigma1, sigma2a, niter, chisq, npoints, npars = coarse_decor_sigmas(idx1, idx2, mag0, emag0, sigma1, sigma2, verbose = True,)
m, s1, sigma1, sigma2b, niter, chisq, npoints, npars = new_coarse_decor_sigmas(idx1, idx2, mag0, emag0, sigma1, sigma2, verbose = True)

ax = plt.subplot(211)
plt.plot(s0, '.')
plt.plot(s1, '.')
plt.subplot(212, sharex=ax)
plt.plot(s0 - s1, '.')
plt.show()

ax = plt.subplot(211)
plt.plot(sigma2a, '.')
plt.plot(sigma2b, '.')
plt.subplot(212, sharex=ax)
plt.plot(sigma2a - sigma2b, '.')
plt.show()

exit()
x0 = s[7961]
y0 = sigma2[7961]
print x0, y0
plt.plot(s, '.')
plt.show()

nx = 50
ny = 51

xgrid = np.linspace(-.055, -.035, nx + 1)
ygrid = np.linspace(0, .05, ny + 1)

x = (xgrid[:-1] + xgrid[1:])/2.
y = (ygrid[:-1] + ygrid[1:])/2.

chisq = np.zeros((nx, ny))
sigmaf = np.zeros((nx, ny))
sigmaf2 = np.zeros((nx, ny))
cloudf = np.zeros((nx, ny))
logL = np.zeros((nx ,ny))

for i in range(nx):
    s[7961] = x[i]
    ressq = (mag0 - m[idx1] - s[idx2])**2
    for j in range(ny):
        
        sigma2[7961] = y[j]
        
        chisq[i,j] = np.sum((mag0 - m[idx1] - s[idx2])**2/(emag0**2 + sigma1[idx1]**2 + sigma2[idx2]**2))
        
        logL[i,j] = chisq[i,j] + np.sum(np.log(emag0**2 + sigma1[idx1]**2 + sigma2[idx2]**2))
        
        weights = 1/(emag0**2 + (sigma1[idx1])**2 + (y[j])**2)
        
        term = np.bincount(idx2, ressq*weights**2 - weights)
        
        sigmaf[i,j] = term[7961]
        
        term = np.bincount(idx2, (mag0 - m[idx1] - s[idx2])/(emag0**2 + sigma1[idx1]**2 + sigma2[idx2]**2))
        
        cloudf[i,j] = term[7961]
   
plt.figure(figsize=(10,10))
        
plt.subplot(221)
plt.title(r'$\chi^2$')
plt.pcolormesh(xgrid, ygrid, chisq.T)
plt.xlim(-.055, -.035)
plt.ylim(0, .05)
plt.xlabel(r'$s$')
plt.ylabel(r'$\sigma$')

plt.subplot(222)
plt.title(r'$-\ln(L)$')
plt.pcolormesh(xgrid, ygrid, logL.T)
plt.contour(x, y, logL.T, [np.amin(logL)+1], colors='k')
plt.scatter(x0, y0, marker='+', s=40, c='w')
plt.xlim(-.055, -.035)
plt.ylim(0, .05)
plt.xlabel(r'$s$')
plt.ylabel(r'$\sigma$')

plt.subplot(223)
plt.title(r'$\partial_{\sigma}\ln(L)$')
plt.pcolormesh(xgrid, ygrid, sigmaf.T)
plt.contour(x, y, sigmaf.T, [0.], colors='k')
plt.scatter(x0, y0, marker='+', s=40, c='w')
plt.xlim(-.055, -.035)
plt.ylim(0, .05)
plt.xlabel(r'$s$')
plt.ylabel(r'$\sigma$')

plt.subplot(224)
plt.title(r'$\partial_{s}\ln(L)$')
plt.pcolormesh(xgrid, ygrid, cloudf.T)
plt.contour(x, y, cloudf.T, [0.], colors='k')
plt.scatter(x0, y0, marker='+', s=40, c='w')
plt.xlim(-.055, -.035)
plt.ylim(0, .05)
plt.xlabel(r'$s$')
plt.ylabel(r'$\sigma$')

plt.tight_layout()
plt.show()

exit()
