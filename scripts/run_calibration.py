# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 16:24:07 2017

@author: talens
"""

from mpcubed import cdecor_vmag

def main(filelist, aper, nprocs):
    
    for filename in filelist:
        
        cdecor_vmag.CoarseDecorrelation(filename, aper, nprocs=nprocs)
    
    return
    
if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Perform the coarse decorrelation on a baseline.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the file(s) to perform the coarse decorrelation on')
    parser.add_argument('-a', '--aper', type=int, choices=[0,1], default=0,
                        help ='the aperture to perform the coarse decorrelation on', dest='aper')
    parser.add_argument('-n', '--nprocs', type=int, default=4,
                        help='the number of processes to use', dest='nprocs')
    args = parser.parse_args()

    main(args.files, args.aper, args.nprocs)
