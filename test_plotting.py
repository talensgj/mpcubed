#!/usr/bin/env python
# -*- coding: utf-8 -*-

from package import plotting

#f = plotting.fLCplot('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5')
#f.plot_coverage()
#f.plot_star('807144')
#f.plot_as_array('807144', 'peak')

f = plotting.SysPlot('/data2/talens/2015Q2/LPW/test.hdf5')
#f.plot_magnitudes()
#f.plot_trans()
#f.plot_intrapix()
f.plot_clouds()
