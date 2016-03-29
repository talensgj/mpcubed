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
     
    outpath = os.path.join(outpath, camera)
    ensure_dir(outpath)
        
    outfile = os.path.join(outpath, 'fLC_%s%s%s.hdf5'%(date, part, camera))
    filelist = [os.path.join(filepath, '%s/fLC/fLC_%s.hdf5'%(date, date)) for date in dates]
    
    IO.make_baseline(filelist, outfile)
    
    return outfile

def twoweek_baseline_mp(jobs, nprocs):
    
    pool = mp.Pool(nprocs)
    for i in range(len(jobs)):
        pool.apply_async(twoweek_baseline, args = jobs[i])
    pool.close()
    pool.join()
    
    return

def main(date, cameras, mode, filepath, outpath, nprocs=6):

    jobs = []

    for cam in cameras:
        jobs.append((date, cam, mode, filepath, outpath))

    if (len(jobs) == 1):
        twoweek_baseline(*jobs[0])
    elif (nprocs == 1):
        for i in range(len(jobs)):
            twoweek_baseline(*jobs[i])
    else:
        twoweek_baseline_mp(jobs, nprocs)
        
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
