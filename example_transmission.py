#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from package import IO
from package.coordinates import grids
from package.core import cdecor
from package import misc

import matplotlib.pyplot as plt

f = IO.fLCfile('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5')

ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])
nobs = nobs.astype('int')

pg = grids.PolarGrid(13500, 720)
idx1, idx2 = pg.radec2idx(ra, dec)

arg, = np.where(ascc == '807144')
print idx2[arg] 

select = (idx2 == 451)

flux, eflux, lstseq, lst, flag, x, y = f.read_data(['flux0', 'eflux0', 'lstseq', 'lst', 'flag', 'x', 'y'], ascc[select], nobs[select])
lstseq = lstseq.astype('int')

mag, emag = misc.flux2mag(flux, eflux)

nstars = len(ascc[select])
staridx = np.repeat(np.arange(nstars), nobs[select])

lstseq, timeidx = np.unique(lstseq, return_inverse=True)
ntimes = len(lstseq)

ra = np.repeat(ra[select], nobs[select])
dec = np.repeat(dec[select], nobs[select])
ha = np.mod(lst*15. - ra, 360.)

camidx, decidx = pg.radec2idx(ha, dec)

here = (flux >0) & (eflux > 0) & (flag < 1)
x = x[here]
y = y[here]
mag = mag[here]
emag = emag[here]
lst = lst[here]
staridx = staridx[here]
camidx = camidx[here]
timeidx = timeidx[here]

staridx, idx1 = np.unique(staridx, return_inverse=True)
camidx, idx2 = np.unique(camidx, return_inverse=True)

m, z1, quality = cdecor.cdecor(idx1, idx2, mag, emag)

m, z2, A, quality = cdecor.cdecor_intrapix(idx1, idx2, idx2//50, mag, emag, x, y)

print quality

array = np.full((nstars, ntimes), fill_value = np.nan)
array[staridx[idx1], timeidx] = mag - m[idx1] - z2[idx2]

plt.imshow(array, origin='lower', interpolation='None', aspect='auto', vmin=0, vmax=2)
plt.show()

#lstseq, time_idx = np.unique(lstseq, return_inverse=True)

#nstars = len(ascc[select])
#ntimes = len(lstseq)

#star_idx = np.repeat(np.arange(nstars), nobs[select])

#array = np.full((nstars, ntimes), fill_value=np.nan)
#array[star_idx, time_idx] = flux0

#array = array/np.nanmean(array, axis=1, keepdims=True)

#plt.imshow(array, origin='lower', interpolation='None', aspect='auto', vmin=0, vmax=2)
#plt.show()

#ra = np.repeat(ra[select], nobs[select])
#dec = np.repeat(dec[select], nobs[select])
#ha = np.mod(lst*15. - ra, 360.)

#nstars = len(ascc[select])
#star_idx = np.repeat(np.arange(nstars), nobs[select])
#dayidx = lstseq // 13500
#dayidx, day_idx = np.unique(dayidx, return_inverse=True)

#xidx = np.ravel_multi_index([day_idx, star_idx], (16, nstars))

#idx1, idx2 = pg.radec2idx(ha, dec)

#array = np.full((len(np.unique(xidx)), 13500), fill_value=np.nan)
#_, xidx = np.unique(xidx, return_inverse=True)
#array[xidx, idx1] = flux0
#array = array/np.nanmean(array, axis=1, keepdims=True)

#plt.imshow(array, origin='lower', interpolation='None', aspect='auto', vmin=0, vmax=2)
#plt.show()


