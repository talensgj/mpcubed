#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from package import plotting

systematics = plotting.SysPlot('/data2/talens/2015Q2_vmag/LPC/sys0_vmag_201505BLPC.hdf5')

#systematics.plot_magnitudes(savefig=True)
systematics.plot_trans(savefig=True)
systematics.plot_intrapix(savefig=True)
