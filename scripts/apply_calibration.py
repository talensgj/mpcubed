# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 16:24:07 2017

@author: talens
"""

from mpcubed import apply_calibration

def main(filelist, aper):
    
    for filename in filelist:
        
        f = apply_calibration.CorrectLC(filename, aper)
        f.make_redfile()
    
    return
    
if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Perform the coarse decorrelation on a baseline.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the file(s) to perform the coarse decorrelation on')
    parser.add_argument('-a', '--aper', type=int, choices=[0,1], default=0,
                        help ='the aperture to perform the coarse decorrelation on', dest='aper')
    args = parser.parse_args()

    main(args.files, args.aper)