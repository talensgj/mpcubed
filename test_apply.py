#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from package.red_apply_dev import SysCorr
from package import IO
from package import misc
from package import plotting

import matplotlib.pyplot as plt

#f = plotting.SysPlot('/data2/talens/2015Q2/LPE/test_06A_iter1.hdf5')
#f.plot_trans()
#f.plot_intrapix()

#f = plotting.SysPlot('/data2/talens/2015Q2/LPE/test_06A.hdf5')
#f.plot_trans()
#f.plot_intrapix()

f = SysCorr('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5', 0, '/data2/talens/2015Q2/LPE/test_06A.hdf5')
trans, intrapix, clouds, flags = f.get_correction('807144')
data = f.write_corrected_lightcurves('/data2/talens/2015Q2/LPE/test_tmp_06A.hdf5')


