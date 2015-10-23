#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid, HealpixGrid, CartesianGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile
from usefull_functions_dev import mad

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis
import matplotlib.gridspec as gridspec

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

f = fLCfile('/data2/talens/3mEast/fLC_20150611LPE.hdf5')
ascc, ra, dec, nobs, vmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])

with h5py.File('/data2/talens/3mEast/variations/cam_20150611LPE+flags+ipxcartesian+weights.hdf5', 'r') as f:
    
    staridx = f['data/staridx'].value
    sigma = f['data/sigma'].value

here = sigma>0
    
plt.semilogy(vmag[staridx[here]], sigma[here], '.', alpha=.2)
plt.ylim(1e-3, 1)
plt.xlim(8.4, 2)
plt.xlabel('V')
plt.ylabel('MAD')
plt.show()
    
    
