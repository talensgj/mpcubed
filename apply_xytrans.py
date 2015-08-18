#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import os
import glob

import matplotlib.pyplot as plt
from coordinate_grids import PolarGrid

with h5py.File('/data2/talens/Jul2015/Trans0716LPC_pg2700x720.hdf5') as f:
    bins = f['Data/binnum'].value
    trans = f['Data/trans'].value

pg = PolarGrid(2700, 720)
trans = pg.put_values_on_grid(trans, ind=bins, fill_value=np.nan)
trans = np.ravel(trans)

with h5py.File('/data2/talens/Jul2015/fLC_20150714LPC.hdf5') as f:
    with h5py.File('/data2/talens/Jul2015/red_2015050714LPC.hdf5') as g:
    
        ascc = f['table_header/ascc'].value
        
        for sid in ascc:

            si = f['header/'+sid]
            lc = f['data/'+sid]

            ha = np.mod(lc['lst']*15.-np.repeat(si['ra'], si['nobs'].astype('int')), 360.)
            dec = np.repeat(si['dec'], si['nobs'].astype('int'))

            binnum = pg.find_gridpoint(ha, dec)
        
            trans0 = trans[binnum]
            cflux0 = lc['flux0']/trans0
            
            g.create_dataset('data/'+sid+'/trans0', data=trans0)
            g.create_dataset('data/'+sid+'/cflux0', data=cflux0)
