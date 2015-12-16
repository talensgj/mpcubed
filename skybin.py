#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from fLCfile import fLCfile
from core.coarse_decor_dev import coarse_decor, coarse_decor_sigmas
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

m, s, sigma1, sigma2, niter, chisq, npoints, npars = coarse_decor_sigmas(idx1, idx2, mag0, emag0, sigma1, sigma2, verbose = True)
print s[7961], sigma2[7961]
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
        
plt.subplot(221)
plt.pcolormesh(xgrid, ygrid, chisq.T)
plt.xlim(-.055, -.035)
plt.ylim(0, .05)

plt.subplot(222)
plt.pcolormesh(xgrid, ygrid, logL.T)
plt.contour(x, y, logL.T, [np.amin(logL)+1], colors='k')

plt.subplot(223)
plt.pcolormesh(xgrid, ygrid, sigmaf.T)
plt.contour(x, y, sigmaf.T, [0.], colors='k')

plt.subplot(224)
plt.pcolormesh(xgrid, ygrid, cloudf.T)
plt.contour(x, y, cloudf.T, [0.], colors='k')

plt.show()

exit()

# 12, 7961

nx = 50
ny = 51

xgrid = np.linspace(-.1, 0, nx + 1)
ygrid = np.linspace(0, .05, ny + 1)

x = (xgrid[:-1] + xgrid[1:])/2.
y = (ygrid[:-1] + ygrid[1:])/2.



for i in range(nx):
    s[7961] = x[i]
    
    for j in range(ny):
        sigma2[7961] = y[j]
        
        prefactor[i,j] = np.sum(np.log(emag0**2 + sigma1[idx1]**2 + sigma2[idx2]**2))
        chisq[i,j] = np.sum((mag0 - m[idx1] - s[idx2])**2/(emag0**2 + sigma1[idx1]**2 + sigma2[idx2]**2))

logL = prefactor + chisq
logL = logL - np.amin(logL)

print np.amin(logL), np.amax(logL)

plt.subplot(111)
plt.pcolormesh(xgrid, ygrid, logL.T, vmin=0, vmax=10)
plt.contour(x, y, logL.T, np.linspace(0, 1, 10), colors = 'k')
plt.scatter([-0.0437, -0.0454], [0.0, .0322])
plt.show()

exit()

plt.subplot(121)
plt.imshow(prefactor)
plt.subplot(122)
plt.imshow(chisq)
plt.show()

plt.pcolormesh(ygrid, xgrid, -(chisq + prefactor))
plt.colorbar()
plt.scatter([0.0, .0322], [-0.0454, -0.0437])
plt.show()
