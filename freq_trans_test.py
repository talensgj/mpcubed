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

ra = np.repeat(ra[select], nobs[select])
dec = np.repeat(dec[select], nobs[select])
ha = np.mod(lst*15. - ra, 360.)

pgipx = grids.PolarGrid(270, 720)

haidx, decidx = pgipx.radec2idx(ha, dec)
haidx, idx = np.unique(haidx, return_inverse=True)

# Transmission?
freq = 1/360.

sin_mat = np.zeros((len(mag), 5))
sin_mat[:,0] = np.sin(2*np.pi*freq*ha)
sin_mat[:,1] = np.sin(4*np.pi*freq*ha)
sin_mat[:,2] = np.sin(6*np.pi*freq*ha)
sin_mat[:,3] = np.sin(8*np.pi*freq*ha)
sin_mat[:,4] = np.sin(10*np.pi*freq*ha)

cos_mat = np.zeros((len(mag), 5))
cos_mat[:,0] = np.cos(2*np.pi*freq*ha)
cos_mat[:,1] = np.cos(4*np.pi*freq*ha)
cos_mat[:,2] = np.cos(6*np.pi*freq*ha)
cos_mat[:,3] = np.cos(8*np.pi*freq*ha)
cos_mat[:,4] = np.cos(10*np.pi*freq*ha)

# Intrapixel variations.
snx_mat = np.zeros((len(mag), np.amax(idx)+1))
snx_mat[np.arange(len(mag)), idx] = np.sin(2*np.pi*x)

csx_mat = np.zeros((len(mag), np.amax(idx)+1))
csx_mat[np.arange(len(mag)), idx] = np.cos(2*np.pi*x)

sny_mat = np.zeros((len(mag), np.amax(idx)+1))
sny_mat[np.arange(len(mag)), idx] = np.sin(2*np.pi*y)

csy_mat = np.zeros((len(mag), np.amax(idx)+1))
csy_mat[np.arange(len(mag)), idx] = np.cos(2*np.pi*y)

n = 5
tra_mat = np.hstack([sin_mat[:,:n], cos_mat[:,:n]])

here = (flux > 0) & (eflux > 0) & (flag < 1)
tra_mat = tra_mat[here]
mag = mag[here]
staridx = staridx[here]

fit2 = np.zeros(len(mag))
for i in range(5):

    pars1 = np.bincount(staridx, mag - fit2)/np.bincount(staridx)
    fit1 = pars1[staridx]
    
    pars2 = np.linalg.lstsq(tra_mat, mag - fit1)[0]
    fit2 = np.dot(tra_mat, pars2)
    
print np.sum((mag - fit1 - fit2)**2)
    
for i in range(nstars):
    here = staridx == i

    plt.subplot(211)
    plt.title(i)
    plt.plot(mag[here], '.')
    plt.plot(fit[here], '.')
    plt.subplot(212)
    plt.plot((mag - fit)[here], '.')
    plt.show()


