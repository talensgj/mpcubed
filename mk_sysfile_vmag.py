#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing as mp

from package.red_decor_vmag import CoarseDecorrelation
from package.red_apply_vmag import CorrectLC

def systematics(filename):
    
    CoarseDecorrelation(filename, 0)
    
    return

def correct(LBfile):
    
    # Create reduced lightcurves.
    f = CorrectLC(LBfile, 0)
    redfile = f.make_redfile()

    return

data = ['/data2/talens/2015Q2/LPN/fLC_201504ALPN.hdf5',
        '/data2/talens/2015Q2/LPN/fLC_201504BLPN.hdf5',
        '/data2/talens/2015Q2/LPN/fLC_201505ALPN.hdf5',
        '/data2/talens/2015Q2/LPN/fLC_201505BLPN.hdf5',
        '/data2/talens/2015Q2/LPN/fLC_201506ALPN.hdf5',
        '/data2/talens/2015Q2/LPN/fLC_201506BLPN.hdf5',
        '/data2/talens/2015Q2/LPS/fLC_201504ALPS.hdf5',
        '/data2/talens/2015Q2/LPS/fLC_201504BLPS.hdf5',
        '/data2/talens/2015Q2/LPS/fLC_201505ALPS.hdf5',
        '/data2/talens/2015Q2/LPS/fLC_201505BLPS.hdf5',
        '/data2/talens/2015Q2/LPS/fLC_201506ALPS.hdf5',
        '/data2/talens/2015Q2/LPS/fLC_201506BLPS.hdf5']

pool = mp.Pool(processes = 6)
for i in range(12):
    pool.apply_async(correct, args = (data[i],))
pool.close()
pool.join()
