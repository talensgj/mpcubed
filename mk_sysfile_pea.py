#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp

from red_decor_pea import CoarseDecorrelation
from red_apply_pea import CorrectLC

def systematics(filename, aper):
    
    # Perform coarse decorrelation.
    CoarseDecorrelation(filename, aper)
    
    return

def correct(filename, aper):
    
    # Create reduced lightcurves.
    f = CorrectLC(filename, aper)
    redfile = f.make_redfile()

    return

#def systematics_mp(jobs, nprocs):
    
    #pool = mp.Pool(nprocs)
    #for i in range(len(jobs)):
        #pool.apply_async(systematics, args=jobs[i])
    #pool.close()
    #pool.join()
    
    #return
    
#def correct_mp(jobs, nprocs):
    
    #pool = mp.Pool(nprocs)
    #for i in range(len(jobs)):
        #pool.apply_async(correct, args=jobs[i])
    #pool.close()
    #pool.join()
    
    #return

#def main(filelist, aper, mode, nprocs):

    #jobs = []
    #for filename in filelist:
        #jobs.append((filename, aper))

    #if (len(jobs) == 1):
        #if (mode == 0):
            #systematics(*jobs[0])
        #else:
            #correct(*jobs[0])
    #elif (nprocs == 1):
        #for i in range(len(jobs)):
            #if (mode == 0):
                #systematics(*jobs[i])
            #else:
                #correct(*jobs[i])
    #else:
        #if (mode == 0):
            #systematics_mp(jobs, nprocs)
        #else:
            #correct_mp(jobs, nprocs)
            
    #return
    
def systematics_mp(queue):
    
    while True:
        
        item = queue.get()
        
        if (item == 'DONE'):
            break
        else:
            systematics(*item)
    
    return
    
def correct_mp(queue):
    
    while True:
        
        item = queue.get()
        
        if (item == 'DONE'):
            break
        else:
            correct(*item)
    
    return
    
def main(filelist, aper, mode, nprocs):
    
    the_queue = mp.Queue(nprocs)
    
    if (mode == 0):
        the_pool = mp.Pool(nprocs, systematics_mp, (the_queue,))
    elif (mode == 1): 
        the_pool = mp.Pool(nprocs, correct_mp, (the_queue,))
    else:
        print 'Received an invalid value for mode.'
        return
    
    for filename in filelist:
        the_queue.put((filename, aper))
    
    for i in range(nprocs):
        the_queue.put('DONE')
        
    the_pool.close()
    the_pool.join()
        
    return

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Perform the coarse decorrelation on a baseline.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the file(s) to perform the coarse decorrelation on')
    parser.add_argument('-a', '--aper', type=int, choices=[0,1], default=0,
                        help ='the aperture to perform the coarse decorrelation on', dest='aper')
    parser.add_argument('-m', '--mode', type=int, choices=[0,1], default=0,
                        help ='compute or apply the coarse decorrelation', dest='mode')
    parser.add_argument('-n', '--nprocs', type=int, default=6,
                        help='the number of processes to use', dest='nprocs')
    args = parser.parse_args()
    
    main(args.files, args.aper, args.mode, args.nprocs)
