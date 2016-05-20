#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

from transit_search import read_header, read_data

def main():

    data = glob.glob('/data2/talens/red2015/*LPC.hdf5')
    data = np.sort(data)

    ascc, ra, dec, vmag, sptype, jdmin, jdmax = read_header(data)

    select = (vmag < 4) & (vmag > 3.5) & (dec > 45) & (dec < 55)
    ascc = ascc[select]

    print ascc

    lstseq, jdmid, lst, mag, emag, mask, trend, cam = read_data(data, ascc)

    sel = np.any(~mask, axis=1)
    print sel

    with h5py.File('/home/talens/similar_bpic_LPC.hdf5') as f:
        
        grp = f.create_group('header')
        grp.create_dataset('ascc', data=ascc[sel])
        grp.create_dataset('ra', data=ra[sel])
        grp.create_dataset('dec', data=dec[sel])
        grp.create_dataset('vmag', data=vmag[sel])
        grp.create_dataset('sptype', data=sptype[sel])

        for i in range(len(ascc)):

            if np.all(mask[i]): continue

            plt.plot(jdmid[~mask[i]], mag[i, ~mask[i]], '.', c='k')
            plt.plot(jdmid[~mask[i]], trend[i, ~mask[i]], '.', c='r')
            plt.show()
            plt.close()
                
            grp = f.create_group('data/'+ascc[i])
            grp.create_dataset('jdmid', data=jdmid[~mask[i]])
            grp.create_dataset('lst', data=lst[~mask[i]])
            grp.create_dataset('mag', data=mag[i, ~mask[i]])
            grp.create_dataset('emag', data=emag[i, ~mask[i]])
            grp.create_dataset('trend', data=trend[i, ~mask[i]])
            
    return

if __name__ == '__main__':
    main()
