#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid, CartesianGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

import multiprocessing as mp
import time

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

def decbin(decidx, staridx, hauni, mag, emag):
    
    return (decidx, systematics_dev.trans(staridx, hauni, mag, emag, verbose=False, use_weights=False))

def single(decidx, staridx, hauni, mag, emag):
    decuni = np.unique(decidx)
    
    for ind in decuni:
        decbin(ind, staridx[decidx==ind], hauni[decidx==ind], mag[decidx==ind], emag[decidx==ind])
    
    return
    

def multi(processes, decidx, staridx, hauni, mag, emag):
    
    decuni = np.unique(decidx)
    
    pool = mp.Pool(processes=processes)
    results = [pool.apply_async(decbin, args=(ind, staridx[decidx==ind], hauni[decidx==ind], mag[decidx==ind], emag[decidx==ind])) for ind in decuni]
    results = [p.get() for p in results]
    results.sort() # to sort the results by input window width
    
    return

f = fLCfile('/data2/talens/3mEast/fLC_20150611LPE.hdf5')
ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])
lst, flux, eflux, sky, flag = f.read_data(['lst', 'flux0', 'eflux0', 'sky', 'flag'], ascc, nobs)

nobs = nobs.astype('int')

staridx = np.repeat(np.arange(len(ascc)), nobs)

ha = np.mod(lst*15 - np.repeat(ra, nobs), 360.)
pg = PolarGrid(13500, 720)
haidx, hauni = pg.find_gridpoint(ha, np.repeat(dec, nobs), compact=True)
decidx = pg.find_decidx(dec)
decidx = np.repeat(decidx, nobs)

here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
lst = lst[here]
staridx = staridx[here]
hauni = hauni[here]
decidx = decidx[here]
flux = flux[here]
eflux = eflux[here]

mag = -2.5*np.log10(flux)
emag = 2.5/np.log(10)*eflux/flux

start = time.time()
single(decidx, staridx, hauni, mag, emag)
print time.time() - start

start = time.time()
multi(4, decidx, staridx, hauni, mag, emag)
print time.time() - start
