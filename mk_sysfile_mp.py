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

inj = ['/data2/talens/2015Q2/LPW/fLC_201504ALPW.hdf5',
       '/data2/talens/2015Q2/LPW/fLC_201504BLPW.hdf5',
       '/data2/talens/2015Q2/LPW/fLC_201505ALPW.hdf5',
       '/data2/talens/2015Q2/LPW/fLC_201505BLPW.hdf5',
       '/data2/talens/2015Q2/LPW/fLC_201506ALPW.hdf5',
       '/data2/talens/2015Q2/LPW/fLC_201506BLPW.hdf5']

pool = mp.Pool(processes = 6)
for i in range(6):
    pool.apply_async(systematics, args = (inj[i],))
pool.close()
pool.join()
