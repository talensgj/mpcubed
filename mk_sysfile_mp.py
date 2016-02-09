#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp

from package.red_decor import CoarseDecorrelation
from package.red_apply import CorrectLC

def systematics(LBfile):
    
    # Perform coarse decorrelation.
    CoarseDecorrelation(LBfile, 0)
    
    return

def correct(LBfile):
    
    # Create reduced lightcurves.
    f = CorrectLC(LBfile, 0)
    redfile = f.make_redfile()

    return
    
def reduction(LBfile):
    
    # Perform coarse decorrelation.
    CoarseDecorrelation(LBfile, 0)
    
    # Create reduced lightcurves.
    f = CorrectLC(LBfile, 0)
    redfile = f.make_redfile()

    return

data = ['/data2/talens/2015Q2/LPC/fLC_201504ALPC.hdf5',
        '/data2/talens/2015Q2/LPC/fLC_201504BLPC.hdf5',
        '/data2/talens/2015Q2/LPC/fLC_201505ALPC.hdf5',
        '/data2/talens/2015Q2/LPC/fLC_201505BLPC.hdf5',
        '/data2/talens/2015Q2/LPC/fLC_201506ALPC.hdf5',
        '/data2/talens/2015Q2/LPC/fLC_201506BLPC.hdf5']

pool = mp.Pool(processes = 6)
for i in range(6):
    pool.apply_async(correct, args = (data[i],))
pool.close()
pool.join()
