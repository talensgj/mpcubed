# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 16:24:07 2017

@author: talens
"""

from mpcubed import summarize

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Make figures of the candidates for the database.')
    parser.add_argument('blsdir', type=str,
                        help='the BLS run to create the figures for')
    parser.add_argument('aper', type=int,
                        help='the aperture to use when making the figures')
    parser.add_argument('method', type=str,
                        help='the dtrending method to use when making the figures')
    args = parser.parse_args()

    summarize.boxlstsq_summary(args.blsdir, args.aper, args.method)
