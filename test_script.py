#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp

from package import IO
from package.red_decor import CoarseDecorrelation
from package.red_apply import CorrectLC

def ensure_dir(path):
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    return

def twoweek_baseline(date, camera, mode, filepath, outpath):
    
    if (mode == 0):
        part = 'A'
        dates = [date + '%.2i'%i + camera for i in range(1, 16)]
    elif (mode == 1):
        part = 'B'
        dates = [date + '%.2i'%i + camera for i in range(16, 32)]
    elif (mode == 2):
        part = ''
        dates = [date + '%.2i'%i + camera for i in range(1, 32)]
    else:
        print 'Unknown value for mode.'
        print 'exiting...'
        exit()
        
    outfile = os.path.join(outpath, 'fLC_%s%s%s.hdf5'%(date, part, camera))
    filelist = [os.path.join(filepath, '%s/fLC/fLC_%s.hdf5'%(date, date)) for date in dates]
    
    IO.make_baseline(filelist, outfile)
    
    return outfile

LBfile = '/data2/talens/2015Q2/LPE/fLC_201504ALPE.hdf5'
    
# Perform coarse decorrelation.
CoarseDecorrelation(LBfile, 0)

# Create reduced lightcurves.
f = CorrectLC(LBfile, 0)
redfile = f.make_redfile()

print redfile
