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

data = ['/data3/talens/2015Q3/LPE/fLC_201507ALPE.hdf5',
        '/data3/talens/2015Q3/LPE/fLC_201507BLPE.hdf5',
        '/data3/talens/2015Q3/LPE/fLC_201508ALPE.hdf5',
        '/data3/talens/2015Q3/LPE/fLC_201508BLPE.hdf5',
        '/data3/talens/2015Q3/LPE/fLC_201509ALPE.hdf5',
        '/data3/talens/2015Q3/LPS/fLC_201507ALPS.hdf5',
        '/data3/talens/2015Q3/LPS/fLC_201507BLPS.hdf5',
        '/data3/talens/2015Q3/LPS/fLC_201508ALPS.hdf5',
        '/data3/talens/2015Q3/LPS/fLC_201508BLPS.hdf5',
        '/data3/talens/2015Q3/LPS/fLC_201509ALPS.hdf5',
        '/data3/talens/2015Q3/LPW/fLC_201507ALPW.hdf5',
        '/data3/talens/2015Q3/LPW/fLC_201507BLPW.hdf5',
        '/data3/talens/2015Q3/LPW/fLC_201508ALPW.hdf5',
        '/data3/talens/2015Q3/LPW/fLC_201508BLPW.hdf5',
        '/data3/talens/2015Q3/LPW/fLC_201509ALPW.hdf5']

pool = mp.Pool(processes = 15)
for i in range(15):
    pool.apply_async(systematics, args = (data[i],))
pool.close()
pool.join()
