#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np


import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from viridis import viridis
from sysfile import SysFile

f = SysFile('/data2/talens/2015Q2/LPE/nosigmas/sys0_201506LPE.hdf5')
pg, trans0, nobs = f.read_trans()

f = SysFile('/data2/talens/2015Q2/LPE/nosigmas/sys0_201506BLPE.hdf5')
pg, trans1, nobs = f.read_trans()

trans0 = trans0[1:-1,1:-1]
trans0 = np.ma.masked_invalid(trans0)

trans1 = trans1[1:-1,1:-1]
trans1 = np.ma.masked_invalid(trans1)

image = (trans0 - trans1).T

vmin = np.nanpercentile(image, 1)
vmax = np.nanpercentile(image, 99)
print vmin, vmax

plt.pcolormesh(pg.bins1, pg.bins2, image, vmin = vmin, vmax = vmax, cmap = viridis)
plt.colorbar()
plt.xlim(270, 350)
plt.ylim(-20, 65)
plt.show()
