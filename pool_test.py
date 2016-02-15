#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import multiprocessing as mp

def function():
    
    print 'Waiting...'
    
    time.sleep(5)
    
    return 

pool = mp.Pool(processes = 6)
for i in range(14):
    pool.apply_async(function)
pool.close()
pool.join()
