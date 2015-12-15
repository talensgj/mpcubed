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

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201504ALPE_test.hdf5')
pg, trans0, nobs = f.read_trans()

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201504ALPE_test2.hdf5')
pg, trans1, nobs = f.read_trans()

trans0 = trans0[1:-1,1:-1]
trans0 = np.ma.masked_invalid(trans0)

trans1 = trans1[1:-1,1:-1]
trans1 = np.ma.masked_invalid(trans1)

image = (trans0 - trans1).T
image = image - np.nanmean(image, axis=1, keepdims=True)[:, None]

#for i in range(720):
    #print image[i].mask
    #if np.all(image[i].mask): continue

    #ax = plt.subplot(211)
    #plt.plot(trans0[:,i], '.')
    #plt.plot(trans1[:,i], '.')
    #plt.subplot(212, sharex=ax)
    #plt.plot(image[i], '.')
    #plt.ylim(-5e-3, 5e-3)
    #plt.show()
    #plt.close()

#here = np.abs(image) < 1e-3
#image[:,:] = 1
#image[here] = 0

plt.pcolormesh(pg.bins1, pg.bins2, image, vmin = -1e-3, vmax = 1e-3, cmap = viridis)
plt.colorbar()
plt.xlim(270, 350)
plt.ylim(-20, 65)
plt.show()

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201504ALPE.hdf5')
pg, a0, b0, c0, d0, nobs = f.read_intrapix()

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201504ALPE_test.hdf5')
pg, a1, b1, c1, d1, nobs = f.read_intrapix()

b0 = b0[1:-1,1:-1]
b0 = np.ma.masked_invalid(b0)

b1 = b1[1:-1,1:-1]
b1 = np.ma.masked_invalid(b1)

image = (b0 - b1).T

plt.pcolormesh(pg.bins1, pg.bins2, image, cmap = viridis)
plt.colorbar()
plt.xlim(270, 350)
plt.ylim(-20, 65)
plt.show()

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201504ALPE.hdf5')
pg, clouds0, nobs, lstmin, lstmax = f.read_clouds()

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201504ALPE_test.hdf5')
pg, clouds1, nobs, lstmin, lstmax = f.read_clouds()

for i in range(100):
    plt.plot(clouds0[i] - clouds1[i], '.')
    plt.show()
    plt.close()

