#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp

from mpcubed import io, misc

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
    
    io.make_baseline(outfile, filelist)
    
    return outfile

def main(date, cameras, mode, filepath, outpath):

    for cam in cameras:
        twoweek_baseline(date, cam, mode, filepath, outpath)
        
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
    args = parser.parse_args()
    
    main(args.date, args.cameras, args.mode, args.filepath, args.outpath)
