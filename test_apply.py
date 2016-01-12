#!/usr/bin/env python
# -*- coding: utf-8 -*-

from package.red_apply_dev import SysCorr
from package import IO
from package import misc

import matplotlib.pyplot as plt

f = SysCorr('/data2/talens/2015Q2/LPW/fLC_201506ALPW.hdf5', 0, '/data2/talens/2015Q2/LPW/test.hdf5')
trans, intrapix, clouds, flags = f.get_correction('807144')

f = IO.fLCfile('/data2/talens/2015Q2/LPW/fLC_201506ALPW.hdf5')
lc = f.read_star('807144')

mag, emag = misc.flux2mag(lc['flux0'], lc['eflux0'])

plt.subplot(311)
plt.plot(mag, '.')
plt.subplot(312)
plt.plot(clouds + trans + intrapix, '.')
plt.subplot(313)
plt.plot(mag - clouds - trans - intrapix, '.')
plt.show()
