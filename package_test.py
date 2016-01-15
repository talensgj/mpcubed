#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import multiprocessing as mp

from package import detrend

def ensure_dir(path):
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    return


filepath = '/data2/mascara/LaPalma'
outpath = '/data2/talens/2015Q2/LPE'

ensure_dir(outpath)

dates = ['201504', '201504']
mode = [0, 1]

pool = mp.Pool(2)
[pool.apply_async(detrend.create_baseline, args=(dates[i], 'LPE', mode[i], filepath, outpath)) for i in range(2)]
