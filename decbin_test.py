#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid
from core.index_functions import index_statistics

from core import intrapix_dev
from core import coarse_decor

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

f = fLCfile('/data2/talens/3mEast/fLC_20150601LPE.hdf5')
ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])

pg = PolarGrid(13500, 720)
decidx = pg.find_decidx(dec)
print decidx[ascc=='807144']

here = (decidx == 451)
ascc = ascc[here]
ra = ra[here]
dec = dec[here]
nobs = nobs[here]

lst, x, y, flux, eflux = f.read_data(['lst', 'x', 'y', 'flux0', 'eflux0'], ascc, nobs)

nobs = nobs.astype('int')

ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
haidx = pg.find_raidx(ha)

staridx = np.repeat(np.arange(len(ascc)), nobs)

here = (flux > 0)&(eflux > 0)
lst = lst[here]
staridx = staridx[here]
haidx = haidx[here]
x = x[here]
y = y[here]
flux = flux[here]
eflux = eflux[here]

mag = -2.5*np.log10(flux)
emag = 2.5/np.log(10)*eflux/flux

m, z, niter, chisq, npoints, npars = intrapix_dev.coarse_positions(staridx, haidx, mag, emag, verbose=True, use_weights=True)
print niter, chisq, npoints, npars
    
plt.plot(z, '.')
plt.show()
 
