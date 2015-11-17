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

# Initialize reader and coordinate grids.
f = fLCfile('/data2/talens/Orientation/fLC_20150618LPE.hdf5')

# Read header data.
ascc, ra, dec, nobs, vmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
nobs = nobs.astype('int')

tmp1 = 0
tmp2 = 0
for i in range(len(ascc)):
    
    flux, eflux = f.read_data(['flux0', 'eflux0'], [ascc[i]], [nobs[i]])
    
    here = (flux > 0) & (eflux > 0)
    flux = flux[here]
    eflux = eflux[here]
    
    mag = -2.5*np.log10(flux)
    emag = 2.5/np.log(10.)*eflux/flux
    
    weights = 1/emag**2
    
    tmp1 = tmp1 + np.sum(weights*(mag - vmag[i]))
    tmp2 = tmp2 + np.sum(weights)
    
    print -tmp1/tmp2
