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

data = ['/data2/talens/2015Q2/LPE/fLC_201504ALPE.hdf5',
        '/data2/talens/2015Q2/LPE/fLC_201504BLPE.hdf5',
        '/data2/talens/2015Q2/LPC/fLC_201504ALPC.hdf5',
        '/data2/talens/2015Q2/LPC/fLC_201504BLPC.hdf5',
        '/data2/talens/2015Q2/LPW/fLC_201504ALPW.hdf5',
        '/data2/talens/2015Q2/LPW/fLC_201504BLPW.hdf5',
        '/data2/talens/2015Q2/LPE/fLC_201505ALPE.hdf5',
        '/data2/talens/2015Q2/LPE/fLC_201505BLPE.hdf5',
        '/data2/talens/2015Q2/LPC/fLC_201505ALPC.hdf5',
        '/data2/talens/2015Q2/LPC/fLC_201505BLPC.hdf5',
        '/data2/talens/2015Q2/LPW/fLC_201505ALPW.hdf5',
        '/data2/talens/2015Q2/LPW/fLC_201505BLPW.hdf5']

pool = mp.Pool(processes = 6)
for i in range(12):
    pool.apply_async(correct, args = (data[i],))
pool.close()
pool.join()
