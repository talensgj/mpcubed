#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import os
import glob

from index_functions import index_statistics
from coordinate_grids import PolarGrid

def apply_xytrans(fLC, cam):
    
    with h5py.File(cam, 'r') as f:
        
        gridname = f['header'].attrs['grid']
        
        if gridname == 'polar':
            grid = PolarGrid(f['header'].attrs['nx'], f['header'].attrs['ny'])
        else:
            print 'Grid %s not implemented yet.'%gridname
            exit()
        
        binnum = f['data/binnum'].value
        trans = f['data/trans'].value
        count = f['data/count'].value
        
    trans = grid.put_values_on_grid(trans, ind=binnum, fill_value=np.nan)
    trans = np.ravel(trans)

    count = grid.put_values_on_grid(count, ind=binnum, fill_value=0)
    count = np.ravel(count)

    head, tail = os.path.split(fLC)
    tail = 'red_'+tail.rsplit('_')[-1]
    red = os.path.join(head, tail)
    
    with h5py.File(fLC, 'r') as f, h5py.File(red, 'w-') as g:

        ascc = f['table_header/ascc'].value
        
        grp = g.create_group('header')
        grp.attrs['cam'] = cam
        grp = g.create_group('data')
        
        for sid in ascc:

            si = f['header/'+sid]
            lc = f['data/'+sid]

            ha = np.mod(lc['lst']*15.-np.repeat(si['ra'], si['nobs'].astype('int')), 360.)
            dec = np.repeat(si['dec'], si['nobs'].astype('int'))

            binnum = grid.find_gridpoint(ha, dec)
        
            trans0 = trans[binnum]
            cflux0 = lc['flux0']/trans0
            ecflux0 = lc['eflux0']/trans0
            
            sflux0 = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std', keeplength=True)
            scflux0 = sflux0/trans0
            
            flags = np.where(np.isnan(trans0), 1, 0) # Flag data where no transmission coefficient exists.
            flags = flags + np.where(count[binnum]<=5, 2, 0) # Flag data where less than 5 datapoints were used to compute the transmission.
            
            record = np.rec.fromarrays([trans0, cflux0, ecflux0, sflux0, scflux0, flags], names=['trans0', 'cflux0', 'ecflux0', 'sflux0', 'scflux0', 'flags'])
            
            grp.create_dataset(sid, data=record)
    
    return 0
    
apply_xytrans('/data2/talens/Jul2015/fLC_20150714LPC.hdf5', '/data2/talens/Jul2015/cam_20150716LPC_pg2700x720.hdf5')
