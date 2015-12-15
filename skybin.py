#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from fLCfile import fLCfile
from core.coarse_decor import coarse_decor, coarse_decor_sigmas
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

plt.plot(s, '.')
plt.show()
exit()

# 12, 7961

xgrid = np.linspace(-.1, 0, 101)
ygrid = np.linspace(0, .05, 51)

x = (xgrid[:-1] + xgrid[1:])/2.
y = (ygrid[:-1] + ygrid[1:])/2.

chisq = np.zeros((100,50))
prefactor = np.zeros((100,50))
for i in range(100):
    print i
    s[7961] = x[i]
    for j in range(50):
        sigma2[7961] = y[j]
        prefactor[i,j] = np.sum(np.log(emag0**2 + sigma1[idx1]**2 + sigma2[idx2]**2))
        chisq[i,j] = np.sum((mag0 - m[idx1] - s[idx2])**2/(emag0**2 + sigma1[idx1]**2 + sigma2[idx2]**2))

plt.subplot(121)
plt.imshow(prefactor)
plt.subplot(122)
plt.imshow(chisq)
plt.show()

plt.pcolormesh(ygrid, xgrid, -(chisq + prefactor))
plt.colorbar()
plt.scatter([0.0, .0322], [-0.0454, -0.0437])
plt.show()
