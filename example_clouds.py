#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from package import IO
from package.coordinates import grids

import matplotlib.pyplot as plt

f = IO.fLCfile('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5')

ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])
nobs = nobs.astype('int')

hg = grids.HealpixGrid(8)
idx = hg.radec2idx(ra, dec)

arg, = np.where(ascc == '807144')
print idx[arg] 

select = (idx == 266)
print select
flux0, lstseq = f.read_data(['flux0', 'lstseq'], ascc[select], nobs[select])
lstseq = lstseq.astype('int')

lstseq, time_idx = np.unique(lstseq, return_inverse=True)

nstars = len(ascc[select])
ntimes = len(lstseq)

star_idx = np.repeat(np.arange(nstars), nobs[select])

array = np.full((nstars, ntimes), fill_value=np.nan)
array[star_idx, time_idx] = flux0

array = array/np.nanmean(array, axis=1, keepdims=True)

plt.imshow(array, origin='lower', interpolation='None', aspect='auto', vmin=0, vmax=2)
plt.show()

