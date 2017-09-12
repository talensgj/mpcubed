#!/usr/bin/env python
# -*- coding: utf-8 -*-

from mpcubed import transit_search 

def main(args):
    
    transit_search.transit_search(args.files, args.name, args.patches, args.aper, args.method, args.outdir, args.nprocs)

    return

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description='Run the box least-squares on a collection of reduced data.')
    parser.add_argument('files', type=str, nargs='+',
                        help='the file(s) to process')
    parser.add_argument('name', type=str,
                        help ='the name of this box least-squares run')
    parser.add_argument('-p', '--patches', type=int, nargs='+', default=None,
                        help ='the sky patch(es) to process', dest='patches')
    parser.add_argument('-a', '--aper', type=int, default=0,
                        help ='the aperture to use', dest='aper')
    parser.add_argument('-m', '--method', type=str, default='legendre',
                        help ='detrending method', dest='method')
    parser.add_argument('-o', '--outdir', type=str, default='/data3/talens/boxlstsq',
                        help ='the location to write the results', dest='outdir')
    parser.add_argument('-n', '--nprocs', type=int, default=6,
                        help='the number of processes to use', dest='nprocs')
    args = parser.parse_args()
    
    main(args)
