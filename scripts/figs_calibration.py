# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 16:24:07 2017

@author: talens
"""

from mpcubed import summarize

def main(filelist):
    
    for filename in filelist:
        
        summarize.calibration_summary(filename)
    
    return
    
if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Make figures of calibration terms.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the file(s) to create the figures for')
    args = parser.parse_args()

    main(args.files)
