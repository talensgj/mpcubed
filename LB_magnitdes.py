#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid
from core.index_functions import index_statistics

from core import systematics_dev

from fLCfile import fLCfile

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 
import matplotlib.gridspec as gridspec

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/3mEast/LBtests/June1.hdf5', 'r') as f:
    vmag = f['header_table/vmag'].value
    dec = f['header_table/dec'].value
    ra = f['header_table/ra'].value

with h5py.File('/data2/talens/3mEast/LBtests/camip_June1_iter4.hdf5', 'r') as f:
    mz = f['data/magnitudes/m'].value

with h5py.File('/data2/talens/3mEast/LBtests/skyip_June1_iter4.hdf5', 'r') as f:
    ms = f['data/magnitudes/m'].value

plt.subplot(211)
plt.plot(vmag, mz-vmag, '.')
plt.ylim(-.5, .5)

plt.subplot(212)
plt.plot(dec, ms-vmag, '.')
plt.ylim(-.5, .5)

plt.show()

plt.subplot(211)
plt.plot(ra, mz-vmag, '.')
plt.ylim(-.5, .5)

plt.subplot(212)
plt.plot(ra, ms-vmag, '.')
plt.ylim(-.5, .5)

plt.show()

plt.subplot(121)
plt.scatter(ra, dec, c=mz-vmag, vmin=-.1, vmax=.1)

plt.subplot(122)
plt.scatter(ra, dec, c=ms-vmag, vmin=-.1, vmax=.1)

plt.show()

plt.plot(mz, mz-ms, '.')
plt.ylim(-.1, .1)
plt.show()
