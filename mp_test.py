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

filepath = '/data2/mascara/LaPalma'
outpath = '/data2/talens/2015Q2/LPE'

ensure_dir(outpath)

twoweek_baseline('201505', 'LPE', 0, filepath, outpath)

#date = ['201504', '201504', '201505', '201505', '201506', '201506']
#mode = [0, 1, 0, 1, 0, 1]

#pool = mp.Pool(processes = 6)
#for i in range(6):
    #pool.apply_async(twoweek_baseline, args = (date[i], 'LPC', mode[i], filepath, outpath))
#pool.close()
#pool.join()
    
