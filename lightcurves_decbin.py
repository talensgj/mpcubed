#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from core.coordinate_grids import PolarGrid

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

pg1 = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)

with h5py.File('/data2/talens/3mEast/fLC_20150612LPE.hdf5', 'r') as f:
    
    ascc = f['header_table/ascc'].value
    ra = f['header_table/ra'].value
    dec = f['header_table/dec'].value
    nobs = f['header_table/nobs'].value.astype('int')
    vmag = f['header_table/vmag'].value
    
    decidx = pg1.find_decidx(dec)
    print decidx[ascc=='807144']
    
    here = (decidx == 451)
    ascc = ascc[here]
    ra = ra[here]
    dec = dec[here]
    nobs = nobs[here]
    vmag = vmag[here]
    
    for i in range(len(ascc)):
        lc = f['data/'+ascc[i]].value
        
        ha = np.mod(lc['lst']*15. - ra[i], 360.)
        mag = -2.5*np.log10(lc['flux0'])
        
        plt.plot(ha, mag, '.')
        
plt.show()
plt.close()
        
        
