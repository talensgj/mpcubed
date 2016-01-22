#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from package import IO
from package.systematics import cdecor
from package.coordinates import grids
from package import misc

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

f = IO.fLCfile('/data2/talens/inj_signals/reference/fLC_201506ALPE.hdf5')
fields = ['ascc', 'ra', 'dec', 'vmag', 'nobs']
ascc, ra, dec, vmag, nobs = f.read_header(fields)
nobs = nobs.astype('int')

staridx = np.arange(len(ascc))

pgcam = grids.PolarGrid(13500, 720)
raidx, decidx = pgcam.radec2idx(ra, dec)

arg, = np.where(ascc == '807144')
select = (decidx == decidx[arg])

ascc = ascc[select]
ra = ra[select]
dec = dec[select]
nobs = nobs[select]
staridx = staridx[select]

fields = ['flux0', 'eflux0', 'lst', 'jdmid', 'sky', 'flag']
flux0, eflux0, lst, jdmid, sky, flag = f.read_data(fields, ascc, nobs)

staridx = np.repeat(staridx, nobs)
ra = np.repeat(ra, nobs)
dec = np.repeat(dec, nobs)
ha = np.mod(lst*15. -ra, 360.)

haidx, decidx = pgcam.radec2idx(ha, dec)

dayidx = np.floor(jdmid).astype('int')

staridx, idx1 = np.unique(staridx, return_inverse=True)
dayidx, idx3 = np.unique(dayidx, return_inverse=True)

nstars = len(staridx)
nha = 13500
ndays = len(dayidx)

here = (flux0 > 0) & (eflux0 > 0) & (sky > 0) & (flag < 1)
mag0, emag0 = misc.flux2mag(flux0[here], eflux0[here])

m, z, quality = cdecor.cdecor(idx1[here], haidx[here], mag0, emag0)

image = np.full((nstars, ndays, nha), fill_value=np.nan)
image[idx1, idx3, haidx] = misc.flux2mag(flux0) - m[idx1]

image = image.reshape((nstars*ndays, nha))

res = np.full((nstars, ndays, nha), fill_value=np.nan)
res[idx1, idx3, haidx] = misc.flux2mag(flux0) - m[idx1] - z[haidx]

res = res.reshape((nstars*ndays, nha))

fig = plt.figure(figsize=(13,9))

gs = gridspec.GridSpec(4, 2, width_ratios = [15, .5], height_ratios = [2,10,10,10])

ax = plt.subplot(gs[1,0])
im = plt.imshow(image, aspect='auto', vmin=-.5, vmax=.5)

cax = plt.subplot(gs[1,1])
cb = plt.colorbar(im, cax = cax)

plt.subplot(gs[2,0], sharex=ax)
plt.plot(z, '.')
plt.ylim(-.5, .5)

plt.subplot(gs[3,0], sharex=ax)
im = plt.imshow(res, aspect='auto', vmin=-.5, vmax=.5)

cax = plt.subplot(gs[3,1])
cb = plt.colorbar(im, cax = cax)

plt.show()


