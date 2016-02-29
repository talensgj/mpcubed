#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp

from package import IO

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
outpath = '/data3/talens/2015Q3/LPW'

ensure_dir(outpath)

twoweek_baseline('201509', 'LPC', 1, filepath, outpath)

#date = ['201507', '201507', '201508', '201508', '201509', '201509']
#mode = [0, 1, 0, 1, 0, 1]

#pool = mp.Pool(processes = 6)
#for i in range(6):
    #pool.apply_async(twoweek_baseline, args = (date[i], 'LPW', mode[i], filepath, outpath))
#pool.close()
#pool.join()
    
