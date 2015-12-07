#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

import multiprocessing as mp

from core.coordinate_grids import PolarGrid, HealpixGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

from time import time

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

def wrap(part, staruni, camtransuni, mag, emag):

    return part, systematics_dev.trans(staruni, camtransuni, mag, emag, verbose=False)

def serial(parts, staruni, camtransuni, mag, emag):
    return [wrap(i, staruni[i], camtransuni[i], mag[i], emag[i]) for i in range(parts)]
    
def multiprocess(parts, staruni, camtransuni, mag, emag, processes=1):
    
    pool = mp.Pool(processes=processes)
    results = [pool.apply_async(wrap, args=(i, staruni[i], camtransuni[i], mag[i], emag[i])) for i in range(parts)]
    results = [p.get() for p in results]
    results.sort() #
    
    return results

f = fLCfile('/data2/talens/3mEast/fLC_20150716LPE.hdf5')
ascc, ra, dec, nobs = f.read_header(['ascc', 'ra', 'dec', 'nobs'])
nobs = nobs.astype('int')

pg = PolarGrid(13500, 720)
decidx = pg.find_decidx(dec)

amin = np.amin(decidx)
amax = np.amax(decidx)
#parts = amax-amin
parts = 256
print parts
array1 = []
array2 = []
array3 = []
array4 = []

start = time()
for i in range(parts):

    this_part = np.arange(amin + i, amax, parts)
    here = np.in1d(decidx, this_part)
    this_ascc = ascc[here]
    this_ra = ra[here]
    this_dec = dec[here]
    this_nobs = nobs[here]
    
    lstidx, lst, x, y, flux, eflux, sky, flag = f.read_data(['lstidx', 'lst', 'x', 'y', 'flux0', 'eflux0', 'sky', 'flag'], this_ascc, this_nobs)
    lstidx = lstidx.astype('int')

    # Build indices
    staridx = np.repeat(np.arange(len(this_ascc)), this_nobs)

    ha = np.mod(lst*15 - np.repeat(this_ra, this_nobs), 360.)

    camtransidx = pg.find_gridpoint(ha, np.repeat(this_dec, this_nobs))

    #pg2 = PolarGrid(270, 720)
    #intrapixidx = pg2.find_gridpoint(ha, np.repeat(this_dec, this_nobs))

    # Flag bad data.
    here = (flux > 0) & (eflux > 0) & (sky > 0) & (flag < 1)
    x = x[here]
    y = y[here]
    flux = flux[here]
    eflux = eflux[here]

    staridx = staridx[here]
    camtransidx = camtransidx[here]
    #intrapixidx = intrapixidx[here]

    # Convert flux to magnitudes
    mag = 25 - 2.5*np.log10(flux)
    emag = 2.5/np.log(10)*eflux/flux

    # Get unique indices.
    staridx, staruni = np.unique(staridx, return_inverse=True)
    camtransidx, camtransuni = np.unique(camtransidx, return_inverse=True)
    #intrapixidx, intrapixuni = np.unique(intrapixidx, return_inverse=True)
    
    print len(mag)
    
    array1.append(staruni)
    array2.append(camtransuni)
    array3.append(mag)
    array4.append(emag)
print time() - start

start = time()
results = serial(parts, array1, array2, array3, array4)
print results
print time()-start

start = time()
results = multiprocess(parts, array1, array2, array3, array4, processes=8)
print time()-start

        
        
        
