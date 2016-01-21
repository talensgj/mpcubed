#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from package import IO

f = IO.SysFile('/data2/talens/inj_signals/reference/sys0_201504ALPE.hdf5')
hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()

sigma = sigma[np.isfinite(sigma)]
print sigma.size, sum(sigma == 0), sum(sigma == 2)
