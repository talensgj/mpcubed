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
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

fig = plotting.SysPlot('/data2/talens/inj_signals/reference/sys0_201506ALPE.hdf5')
#fig.plot_trans()
fig.plot_clouds()
