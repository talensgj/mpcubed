#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from package import plotting

systematics = plotting.SysPlot('/data2/talens/2015Q2/LPE/sys0_201504ALPE.hdf5')

systematics.plot_magnitudes()
systematics.plot_trans()
systematics.plot_intrapix()
