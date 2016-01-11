#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from package import IO
from package.coordinates import grids

import matplotlib.pyplot as plt

f = IO.fLCfile('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5')

ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])
nobs = nobs.astype('int')

#sys = IO.SysFile('/data2/talens/2015Q2/LPE/sys0_201506ALPE.hdf5')
#pgcam, trans, nobs_trans = sys.read_trans(ravel = True)
#pgipx, a, b, c, d, nobs_ipx = sys.read_intrapix(ravel = True)
#hg, clouds, nobs_clouds, lstmin, lstmax = sys.read_clouds()

select = np.arange(1000,1100)

flux0, lst, lstseq = f.read_data(['flux0', 'lst', 'lstseq'], ascc[select], nobs[select])
lstseq = lstseq.astype('int')

#ra = np.repeat(ra[select], nobs[select])
#dec = np.repeat(dec[select], nobs[select])
#ha = np.mod(lst*15. - ra, 360.)

#camidx, decidx = pgcam.radec2idx(ha, dec)
#ipxidx, decidx = pgipx.radec2idx(ha, dec)
#skyidx = hg.radec2idx(ra, dec)

#trans = trans[camidx, decidx]
#clouds = clouds[skyidx, lstseq - lstmin]

nstars = len(ascc[select])
idx_star = np.repeat(np.arange(nstars), nobs[select])
lstseq, idx_time = np.unique(lstseq, return_inverse=True) 
ntimes = len(lstseq)

array = np.full((nstars, ntimes), fill_value=np.nan)
array[idx_star, idx_time] = flux0

array = array/np.nanmean(array, axis=1, keepdims=True)
plt.imshow(array, origin='lower', aspect='auto', interpolation='None', vmin=0, vmax=2)
plt.show()
