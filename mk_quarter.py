#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

from package import IO

filelist = glob.glob('/data2/talens/2015Q2/LPN/red0_*.hdf5')
print filelist
IO.make_quarter(filelist, '/data2/talens/2015Q2/LPN/red0_2015Q2LPN.hdf5')
