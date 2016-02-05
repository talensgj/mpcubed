#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from package import plotting

systematics = plotting.SysPlot('/data2/talens/2015Q2/LPC/sys0_201506ALPC.hdf5')

systematics.plot_magnitudes(savefig=True)
systematics.plot_trans(savefig=True)
systematics.plot_intrapix(savefig=True)
