#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib import rcParams

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPS.hdf5') as f:
    lc = f['data/1413051']
    
    plt.subplot(211)
    plt.plot(lc['jdmid'], lc['flux0'], '.')
    plt.plot(lc['jdmid'], lc['flux1'], '.')
    plt.subplot(212)
    plt.plot(lc['jdmid'], lc['flux1']/lc['flux0'], '.')
    plt.show()
