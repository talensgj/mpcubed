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

from sysfile import SysFile

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201506LPE.hdf5')
pg, trans, nobs = f.read_trans()

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201506ALPE.hdf5')
pg, transA, nobs = f.read_trans()

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201506BLPE.hdf5')
pg, transB, nobs = f.read_trans()

plt.subplot(211)
plt.plot(trans[:,360], '.')
plt.plot(transA[:,360], '.')
plt.plot(transB[:,360], '.')
plt.ylabel('Transmission')
plt.subplot(212)
plt.plot(transA[:,360] - trans[:,360], '.')
plt.plot(transB[:,360] - trans[:,360], '.')
plt.xlabel('Hour Angle')
plt.ylabel('Transmission')

plt.tight_layout()
plt.show()
plt.close()

plt.subplot(211)
plt.plot(trans[12000], '.')
plt.plot(transA[12000], '.')
plt.plot(transB[12000], '.')
plt.ylabel('Transmission')
plt.subplot(212)
plt.plot(transA[12000] - trans[12000], '.')
plt.plot(transB[12000] - trans[12000], '.')
plt.xlabel('Declination')
plt.ylabel('Transmission')

plt.tight_layout()
plt.show()
plt.close()


f = SysFile('/data2/talens/2015Q2/LPE/sys0_201506LPE.hdf5')
hg, clouds, nobs, lstmin, lstmax = f.read_clouds()

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201506ALPE.hdf5')
hg, cloudsA, nobsA, lstminA, lstmaxA = f.read_clouds()

f = SysFile('/data2/talens/2015Q2/LPE/sys0_201506BLPE.hdf5')
hg, cloudsB, nobsB, lstminB, lstmaxB = f.read_clouds()

plt.subplot(211)
plt.plot(np.arange(lstmin, lstmax + 1), clouds[39], '.')
plt.plot(np.arange(lstminA, lstmaxA + 1), cloudsA[39], '.')
plt.plot(np.arange(lstminB, lstmaxB + 1), cloudsB[39], '.')
plt.ylabel('Clouds')

plt.subplot(212)
plt.plot(np.arange(lstminA, lstmaxA + 1), cloudsA[39] - clouds[39, lstminA-lstmin: lstmaxA-lstmin+1], '.')
plt.plot(np.arange(lstminB, lstmaxB + 1), cloudsB[39] - clouds[39, lstminB-lstmin: lstmaxB-lstmin+1], '.')
plt.xlabel('Time')
plt.ylabel('Clouds')

plt.tight_layout()
plt.show()
plt.close()

plt.subplot(211)
plt.plot(clouds[:, 50], '.')
plt.plot(cloudsA[:, lstmin - lstminA + 50], '.')
plt.ylabel('Clouds')

plt.subplot(212)
plt.plot(clouds[:, 50] - cloudsA[:, lstmin - lstminA + 50], '.')
plt.xlabel('Position')
plt.ylabel('Clouds')

plt.tight_layout()
plt.show()
plt.close()
