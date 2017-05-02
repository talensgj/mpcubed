#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

from mpcubed import IO

def combine_quarter(filelist, outfile):
    
    filepath = os.path.dirname(filelist[0])
    outfile = os.path.join(filepath, outfile)
    
    IO.make_quarter(filelist, outfile)
    
    return 
    
def main(filelist, outfile):
    
    combine_quarter(filelist, outfile)
    
    return
    
if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Combine temporary lightcurve files.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the files to combine')
    parser.add_argument('outfile', type=str,
                        help='the name of the resulting file')
    args = parser.parse_args()
    
    main(args.files, args.outfile)

