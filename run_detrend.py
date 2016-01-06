#!/usr/bin/env python
# -*- coding: utf-8 -*-

from package.fLCmerge import fLCmerge
from package import detrend

def main(args):
    
    filelist = []
    for :
        fLCmerge('201506', 'LPE', 0, filepath, outpath)
        
    for filename in filelist:
        detrend.CoarseDecorrelation(filename, 0)
    
    for filename in filelist:
        detrend.SysCorr(filename, 0)
        
    detrend.make_quarterfile(filelist, redfile)
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
