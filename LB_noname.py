#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os 
import glob

import h5py
import numpy as np
from numpy.lib.recfunctions import stack_arrays

def make_quarterfile(filelist, redfile):
    
    nfiles = len(filelist)
    
    # Create a list of unique stars in the files.
    ascc = np.array([])
    for i in range(nfiles):
        with h5py.File(filelist[i], 'r') as f:
            ascc = np.append(ascc, f.keys())
            
    ascc = np.unique(ascc)
    
    # Read the data.
    for sID in ascc:
        first = True
        for i in range(nfiles):
        
            with h5py.File(filelist[i], 'r') as f:
                
                # Try to read the star.
                try:
                    tmp = f[sID].value
                except:
                    continue
                
                # Add the data to the lightcurve.
                if first:
                    lc = tmp
                    first = False
                else:
                    lc = stack_arrays((lc, tmp), asrecarray=True)
        
        # Write the data to the redfile.
        with h5py.File(redfile) as f:
            for key in lc.dtype.names:
                f.create_dataset('data/' + sID + '/' + key, data = lc[key])

    return
    
if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Apply the systematics corrections to a given fLC file.')
    parser.add_argument('path', help = 'Path to the temporary files.')
    parser.add_argument('aperture', type = int, help = 'The aperture to be merged.')
    parser.add_argument('-o', '--outfile', help = 'The output file.', default = None)
    args = parser.parse_args()
    
    if not os.path.exists(args.path):
        print args.path, 'does not exist.'
        print 'exiting...'
        exit()
    
    expression = os.path.join(args.path, 'tmp%i_*'%args.aperture)
    filelist = glob.glob(expression)
    filelist = np.sort(filelist)
    
    
    
    make_quarterfile(filelist)
