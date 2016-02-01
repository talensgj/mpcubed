#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import gridspec

path = '/data2/talens/inj_signals/signals/'

with h5py.File(path + 'signals_index.hdf5', 'r') as f:
    P = f['P'].value
    depth = f['delta'].value

with h5py.File(path + 'HFperiods_njd5_2BL_nlst5_24h.hdf5', 'r') as f:
    P_rec = f['Prec'].value
    depth_rec = f['depth'].value
    flag = f['flag'].value
    hchisq = f['hchisq'].value
    dchisq = f['dchisq'].value
    
flag = flag + np.where(np.abs(P_rec/P - 1) > .1, 2, 0)
flag = flag + np.where(np.abs(P_rec/1. - 1) < .1, 4, 0)

print sum(flag == 0)

fig = plt.figure()

gs = gridspec.GridSpec(2, 2, width_ratios=[10,1])

plt.subplot(gs[0,0])
im = plt.scatter(P, P_rec, c=flag, vmin=0, vmax=7)
plt.plot([0, 15], [0,15], c='k')
plt.plot([0, 15], [0, 7.5], c='k')
plt.plot([0, 15], [0, 30], c='k')
plt.xlim(0, 15)
plt.ylim(0, 15)

plt.subplot(gs[1,0])
plt.scatter(depth, depth_rec, c=flag, vmin=0, vmax=7)
plt.plot([.005, .03], [.005, .03], c='k')
plt.xlim(.005, .03)
plt.ylim(-.05, .05)

ax = plt.subplot(gs[:,1])
plt.colorbar(im, cax=ax)

plt.show()



plt.hist(depth, bins=np.linspace(.005, .03, 10))
plt.hist(depth[flag == 0], bins=np.linspace(.005, .03, 10))
plt.show()

plt.hist(P, bins=np.linspace(1, 15, 10))
plt.hist(P[flag == 0], bins=np.linspace(1, 15, 10))
plt.show()
