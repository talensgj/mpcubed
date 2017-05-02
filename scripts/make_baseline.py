#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp

from mpcubed import IO, misc

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
     
    outpath = os.path.join(outpath, camera)
    misc.ensure_dir(outpath)
        
    outfile = os.path.join(outpath, 'fLC_%s%s%s.hdf5'%(date, part, camera))
    filelist = [os.path.join(filepath, '%s/fLC/fLC_%s.hdf5'%(date, date)) for date in dates]
    
    IO.make_baseline(filelist, outfile)
    
    return outfile
    
def twoweek_baseline_mp(queue):
    
    while True:
        
        item = queue.get()
        
        if (item == 'DONE'):
            break
        else:
            twoweek_baseline(*item)
            
    return

def main(date, cameras, mode, filepath, outpath, nprocs=6):

    the_queue = mp.Queue(nprocs)
    the_pool = mp.Pool(nprocs, twoweek_baseline_mp, (the_queue,))

    for cam in cameras:
        the_queue.put((date, cam, mode, filepath, outpath))

    for i in range(nprocs):
        the_queue.put('DONE')
    
    the_pool.close()
    the_pool.join()
        
    return

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Combine fLC files to create a basline.')
    parser.add_argument('date', type=str,
                        help='a date of the form YYYYMM')
    parser.add_argument('mode', type=int, choices=[0,1,2],
                        help='integer indicating which baseline to create')
    parser.add_argument('outpath', type=str,
                        help='path where the resulting file(s) will be written')
    parser.add_argument('-d', '--data', type=str, default='/data2/mascara/LaPalma',
                        help='path to the data', dest='filepath')
    parser.add_argument('-c', '--cam', type=str, nargs='+', default=['LPN', 'LPE', 'LPS', 'LPW', 'LPC'],
                        help ='the camera(s) to perform the combination for', dest='cameras')
    parser.add_argument('-n', '--nprocs', type=int, default=6,
                        help='the number of processes to use', dest='nprocs')
    args = parser.parse_args()
    
    main(args.date, args.cameras, args.mode, args.filepath, args.outpath, args.nprocs)
