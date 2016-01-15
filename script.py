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

# Create baseline.
filepath = '/data2/mascara/LaPalma'
outpath = '/data2/talens/2015Q2/LPE'

ensure_dir(outpath)

LBfile = twoweek_baseline('201504', 'LPE', 0, filepath, outpath)

#pool = mp.Pool(5)
#dates = ['201504', '201505', '201505', '201506', '201506']
#mode = [1, 0, 1, 0, 1]
#results = [pool.apply_async(twoweek_baseline, args=(dates[i], 'LPE', mode[i], filepath, outpath)) for i in range(5)] 
#output = [p.get() for p in results]
#print output
    
# Perform coarse decorrelation.
CoarseDecorrelation(LBfile, 0)

# Create reduced lightcurves.
f = CorrectLC(LBfile, aperture)
redfile = f.make_redfile()
