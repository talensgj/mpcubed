#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp

from time import time

from package.red_decor import CoarseDecorrelation
from package.red_apply import CorrectLC

def reduction(LBfile):
    
    # Perform coarse decorrelation.
    CoarseDecorrelation(LBfile, 0)
    
    # Create reduced lightcurves.
    f = CorrectLC(LBfile, 0)
    redfile = f.make_redfile()

    return

ref = ['/data2/talens/inj_signals/reference/fLC_201504ALPE.hdf5',
           '/data2/talens/inj_signals/reference/fLC_201504BLPE.hdf5',
           '/data2/talens/inj_signals/reference/fLC_201505ALPE.hdf5',
           '/data2/talens/inj_signals/reference/fLC_201505BLPE.hdf5',
           '/data2/talens/inj_signals/reference/fLC_201506ALPE.hdf5',
           '/data2/talens/inj_signals/reference/fLC_201506BLPE.hdf5']
           
start = time()
pool = mp.Pool(processes = 6)
for i in range(6):
    pool.apply_async(reduction, args = (ref[i],))
pool.close()
pool.join()
print time() - start

inj = ['/data2/talens/inj_signals/signals/fLC_201504ALPE.hdf5',
           '/data2/talens/inj_signals/signals/fLC_201504BLPE.hdf5',
           '/data2/talens/inj_signals/signals/fLC_201505ALPE.hdf5',
           '/data2/talens/inj_signals/signals/fLC_201505BLPE.hdf5',
           '/data2/talens/inj_signals/signals/fLC_201506ALPE.hdf5',
           '/data2/talens/inj_signals/signals/fLC_201506BLPE.hdf5']

start = time()
procs = [mp.Process(target=reduction, args = (inj[i],)) for i in range(6)]
for p in procs:
    p.start()
for p in procs:
    p.join()
print time() - start
