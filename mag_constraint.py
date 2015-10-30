#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

f = fLCfile('/data2/talens/3mEast/fLC_20150611LPE.hdf5')
pg = PolarGrid(13500, 720)

Mascc, Mra, Mdec, Mnobs, Mv = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
Mnobs = Mnobs.astype('int')

Mstaridx = np.arange(len(Mascc))
decidx, decuni = pg.find_decidx(Mdec, compact=True)
nbins = len(decidx)      

m = np.zeros(len(Mascc))
z = np.full((13502*722,), fill_value=np.nan)

niter = np.zeros(nbins)
chisq = np.zeros(nbins)
npoints = np.zeros(nbins)
npars = np.zeros(nbins)

for ind in range(nbins):
        
    here = (decuni == ind)
    ascc = Mascc[here]
    ra = Mra[here]
    dec = Mdec[here]
    nobs = Mnobs[here]
    staridx = Mstaridx[here]

    lst, flux, eflux, sky, flag = f.read_data(['lst', 'flux0', 'eflux0', 'sky', 'flag'], ascc, nobs)

    # Build indices.    
    staridx = np.repeat(staridx, nobs)
    
    ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
    camtransidx = pg.find_gridpoint(ha, np.repeat(dec, nobs))
    
    # Flag bad data.
    here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    flux = flux[here]
    eflux = eflux[here]
    
    staridx = staridx[here]
    camtransidx = camtransidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    # Get unique indices.
    staridx, staruni = np.unique(staridx, return_inverse=True)
    camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)

    sj = np.zeros(len(camtransidx))
    for i in range(50):
        #mi = (np.bincount(staruni, (mag - sj[camtransuni])/emag**2) + np.bincount(staruni, Mv[staridx][staruni]/emag**2))/(np.bincount(staruni, 1/emag**2) + np.bincount(staruni, 1/emag**2))
        mi = np.bincount(staruni, (mag - sj[camtransuni])/emag**2)/np.bincount(staruni, 1/emag**2)
        sj = np.bincount(camtransuni, (mag - mi[staruni])/emag**2)/np.bincount(camtransuni, 1/emag**2)
        
    offset = np.median(mi-Mv[staridx])
    sj = sj + offset
    mi = mi - offset
            
    m[staridx] = mi
    z[camtransidx] = sj
    
z = z.reshape((13502, 722))
vmin = np.nanpercentile(z, 1)
vmax = np.nanpercentile(z, 99)

plt.plot(Mv, m, '.')
plt.show()

plt.plot(Mdec, Mv - m, '.')
plt.show()

mask = np.all(np.isnan(z), axis=1)
z = z[~mask]
mask = np.all(np.isnan(z), axis=0)
z = z[:, ~mask]

plt.imshow(z.T, aspect='auto', cmap=viridis, vmin=vmin, vmax=vmax)
plt.colorbar()
plt.show()
    
    
    
    
