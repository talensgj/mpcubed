#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import rcParams
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package import IO

f = IO.SysFile('/data2/talens/inj_signals/reference/sys0_201504ALPE.hdf5')

#ascc, vmag, mag, sigma, nobs = f.read_magnitudes()

#ax1 = plt.subplot(221)
#plt.plot(vmag, vmag - mag, '.', alpha=.1)
#plt.xlim(2, 8.4)

#plt.subplot(222, sharey=ax1)
#plt.hist(vmag - mag, bins=np.linspace(-2,2,251), log=True, orientation='horizontal', histtype='step')
#plt.ylim(-2, 2)
#plt.xlim(1,)

#ax2 = plt.subplot(223)
#plt.plot(vmag, sigma, '.', alpha=.1)
#plt.xlim(2, 8.4)

#plt.subplot(224, sharey=ax2)
#plt.hist(sigma, bins=np.linspace(0,.5,251), log=True, orientation='horizontal', histtype='step')
#plt.ylim(0, .5)
#plt.xlim(1,)

#plt.tight_layout()
#plt.show()

#plt.hist2d(vmag - mag, sigma, bins=[np.linspace(-2,2,251), np.linspace(0,.5,251)], norm=colors.LogNorm())
#cb = plt.colorbar()

#plt.xlim(-2, 2)
#plt.ylim(0, .5)

#plt.tight_layout()
#plt.show()

hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()

select = np.where(~np.isnan(clouds))

plt.subplot(121)
plt.hist(clouds[select], bins=np.linspace(-5,5,1001), log=True, histtype='step')
plt.xlim(-5,5)
plt.ylim(1,)

plt.subplot(122)
plt.hist(sigma[select], bins=np.linspace(0,2,1001), log=True, histtype='step')
plt.xlim(0,2)
plt.ylim(1,)

plt.tight_layout()
plt.show()

plt.hist2d(clouds[select], sigma[select], bins=[np.linspace(-5,5,1001), np.linspace(0,2,1001)], norm=colors.LogNorm())
cb = plt.colorbar()

plt.xlim(-5,5)
plt.ylim(0,2)

cb.set_label('Count')
plt.xlabel('clouds')
plt.ylabel('sigma')

plt.tight_layout()
plt.show()
