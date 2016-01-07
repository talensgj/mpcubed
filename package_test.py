#!/usr/bin/env python
# -*- coding: utf-8 -*-

from package import IO
from package import plotting

f = IO.fLCfile('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5')
f = plotting.fLCplot('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5')
