#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from package import IO
from package import plotting
from package.systematics import cdecor
from package.coordinates import grids
from package.red_apply import CorrectLC
from package import misc

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

#f = plotting.SysPlot('/data2/talens/inj_signals/reference/sys0_201506ALPE.hdf5')
#f.plot_trans()
#f.plot_clouds()

f = IO.fLCfile('/data2/talens/inj_signals/reference/fLC_201506ALPE.hdf5')
lc = f.read_star('807144')

g = CorrectLC('/data2/talens/inj_signals/reference/fLC_201506ALPE.hdf5', 0)
trans, ipx, clouds, flags = g.get_correction('807144')

h = IO.SysFile('/data2/talens/inj_signals/reference/sys0_201506ALPE.hdf5')
ascc, vmag, mag, sigma, nobs = h.read_magnitudes()
arg, = np.where(ascc == '807144')
vmag = vmag[arg]
mag = mag[arg]

days = np.floor(lc['jdmid'])
dayidx, = np.where(np.diff(days) > 0)
intidx, = np.where(np.diff(lc['lstseq']) > 1)

plt.figure(figsize=(16,8))
plt.suptitle('Systematics 2015-06 East', size='xx-large')

gs = gridspec.GridSpec(4, 2, width_ratios = [15,.5], height_ratios = [1,10,10,10])

mag0, emag0 = misc.flux2mag(lc['flux0'], lc['eflux0'])

ax1 = plt.subplot(gs[1,0])
ax1.invert_yaxis()
plt.title('ASCC 807144, $V = {:.1f}$'.format(vmag[0]))
plt.plot(mag0, '.')
plt.plot(mag + trans + ipx, '.')
for i in dayidx:
    plt.axvline(i + .5, c='k')
for i in intidx:
    plt.axvline(i + .5, c='k', ls='--')
plt.ylabel(r'Magnitude')

ax2 = plt.subplot(gs[2,0], sharex=ax1)
ax2.invert_yaxis()
plt.plot(mag0 - trans - ipx, '.')
plt.plot(mag + clouds, '.')
for i in dayidx:
    plt.axvline(i + .5, c='k')
for i in intidx:
    plt.axvline(i + .5, c='k', ls='--')
plt.ylabel('Magnitude')

ax3 = plt.subplot(gs[3,0], sharex=ax1)
ax3.invert_yaxis()
plt.plot(mag0 - trans - ipx - clouds, '.')
plt.xlim(-.5, len(mag0) - .5)
for i in dayidx:
    plt.axvline(i + .5, c='k')
for i in intidx:
    plt.axvline(i + .5, c='k', ls='--')
plt.xlabel('Time')
plt.ylabel('Magnitude')

plt.tight_layout()
plt.show()
