#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib import rcParams

from coordinate_grids import PolarGrid

from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

with h5py.File('/data2/talens/Jul2015/transipE.hdf5') as f:
    
    decidx = f['header/decidx'].value
    chisq = f['header/chisq'].value
    niter = f['header/niter'].value
    npoints = f['header/npoints'].value
    npars = f['header/npars'].value
    
    pg = PolarGrid(13500, 720)
    pg1 = PolarGrid(270, 720)
    
    dec = pg.find_dec(decidx)
    
    A = np.full((722, 272), fill_value=np.nan)
    
    for ind in range(len(decidx)):
        
        try:
            data = f['data/%i'%decidx[ind]]
        except:
            pass
        
        haidx_cam = data['haidx_cam'].value
        camtrans = data['camtrans'].value
        
        haidx_ip = data['haidx_ip'].value
        a = data['a'].value
        b = data['b'].value
        
        A[decidx[ind], haidx_ip] = np.sqrt(a**2+b**2)
        
        
        #ha_cam = pg.find_ra(haidx_cam)
        ##ha_cam = np.mod(ha_cam-180, 360)
        
        #ha_ip = pg1.find_ra(haidx_ip)
        ##ha_ip = np.mod(ha_ip-180, 360)
        
        #plt.figure(figsize=(16,8))
        #ax = plt.subplot(311)
        
        #plt.title(r'Dec = %.2f, niter = %i, $\chi^2$ = %.2f'%(dec[ind], niter[ind], chisq[ind]/(npoints[ind]-npars[ind])))
        #plt.plot(ha_cam, camtrans, '.')
        #plt.ylabel('camtrans')
        #plt.subplot(312, sharex=ax)
        #plt.plot(ha_ip, a, '.')
        #plt.ylabel('a')
        #plt.subplot(313, sharex=ax)
        #plt.plot(ha_ip, b, '.')
        #plt.ylabel('b')
        #plt.xlabel('HA [deg]')
        ##plt.xlim(180-50, 180+50)
        #plt.tight_layout()
        #plt.show()
        #plt.close()

A = A[1:-1,1:-1]
#A = np.roll(A, 135, axis=1)

xlim, ylim = np.where(np.isfinite(A))

plt.imshow(A, aspect='auto', vmin=0, vmax=.1, cmap=viridis)
cb = plt.colorbar()
cb.set_label('Amplitude')
plt.xlabel('HA')
plt.ylabel('Dec')
plt.ylim(np.amin(xlim)-.5, np.amax(xlim)+.5)
plt.xlim(np.amin(ylim)-.5, np.amax(ylim)+.5)
plt.show()
plt.close()
