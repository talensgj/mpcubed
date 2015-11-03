#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import os 
import glob
 
import h5py
import numpy as np


#with h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter5_weights.hdf5', 'r') as f, h5py.File('/data2/talens/3mEast/LBtests/camip_15day_iter5_weights_ref.hdf5') as g:
    
    #f.copy('header', g)
    
    #ascc = f['data/ascc'].value
    #m = f['data/m'].value
    #z = f['data/z'].value
    #a = f['data/a'].value
    #b = f['data/b'].value
    #c = f['data/c'].value
    #d = f['data/d'].value
    
    #grp = g.create_group('data')
    
    #grp.create_dataset('magnitudes/ascc', data = ascc)
    #grp.create_dataset('magnitudes/m', data = m)
    
    #idx, = np.where(~np.isnan(z))
    #grp.create_dataset('camtrans/idx', data=idx)
    #grp.create_dataset('camtrans/z', data=z[idx])
    
    #grp['camtrans'].attrs['grid'] = 'polar'
    #grp['camtrans'].attrs['nx'] = 13500
    #grp['camtrans'].attrs['ny'] = 720   
    
    #idx, = np.where(~np.isnan(a))
    #grp.create_dataset('intrapix/idx', data=idx)
    #grp.create_dataset('intrapix/a', data=a[idx])
    #grp.create_dataset('intrapix/b', data=b[idx])
    #grp.create_dataset('intrapix/c', data=c[idx])
    #grp.create_dataset('intrapix/d', data=d[idx])
    
    #grp['intrapix'].attrs['grid'] = 'polar'
    #grp['intrapix'].attrs['nx'] = 270
    #grp['intrapix'].attrs['ny'] = 720
    
with h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter5_weights.hdf5', 'r') as f, h5py.File('/data2/talens/3mEast/LBtests/skyip_15day_iter5_weights_ref.hdf5') as g:
    
    f.copy('header', g)
    
    ascc = f['data/ascc'].value
    m = f['data/m'].value
    s = f['data/s'].value
    #sigma1 = np.full(m.shape, fill_value=np.nan)
    #sigma2 = np.full(s.shape, fill_value=np.nan)
    sigma1 = f['data/sigma1'].value
    sigma2 = f['data/sigma2'].value

    grp = g.create_group('data')
    
    grp.create_dataset('magnitudes/ascc', data=ascc)
    grp.create_dataset('magnitudes/m', data=m)
    grp.create_dataset('magnitudes/sigma', data=sigma1)
    
    idx, lstseq = np.where(~np.isnan(s))
    grp.create_dataset('skytrans/idx', data=idx)
    grp.create_dataset('skytrans/lstseq', data=lstseq)
    grp.create_dataset('skytrans/s', data=s[idx, lstseq])
    grp.create_dataset('skytrans/sigma', data=sigma2[idx, lstseq])
    
    grp['skytrans'].attrs['grid'] = 'healpix'
    grp['skytrans'].attrs['nx'] = 8
    grp['skytrans'].attrs['lstmin'] = 0
    grp['skytrans'].attrs['lstlen'] = 15*13500
