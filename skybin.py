#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from fLCfile import fLCfile
from core.coarse_decor import coarse_decor, coarse_decor_sigmas
from core.coordinate_grids import HealpixGrid
from usefull_functions_dev import flux2mag

import matplotlib.pyplot as plt

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

print len(idx1), len(idx2), len(mag0), len(emag0)
m, s, sigma1, sigma2, niter, chisq, npoints, npars = coarse_decor_sigmas(idx1, idx2, mag0, emag0, verbose = True)
#m, s, niter, chisq, npoints, npars = coarse_decor(idx1, idx2, mag0, emag0, verbose = True)
plt.plot(s, '.')
plt.show()
