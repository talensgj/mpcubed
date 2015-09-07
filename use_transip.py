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

with h5py.File('/data2/talens/Jul2015/transipS.hdf5') as f:
    
    decidx = f['header/decidx'].value
    chisq = f['header/chisq'].value
    niter = f['header/niter'].value
    npoints = f['header/npoints'].value
    npars = f['header/npars'].value
    
    camtrans = np.full((722, 13502), fill_value=np.nan)
    a = np.full((722, 272), fill_value=0)
    b = np.full((722, 272), fill_value=0)
    
    for ind in decidx:
        
        data = f['data/%i'%ind]
        
        camtrans[ind, data['haidx_cam'].value] = data['camtrans'].value
        a[ind, data['haidx_ip'].value] = data['a'].value
        b[ind, data['haidx_ip'].value] = data['b'].value
    
pg = PolarGrid(13500, 720)
pg1 = PolarGrid(270, 720)
    
with h5py.File('/data2/talens/Jul2015/fLC_20150715LPS.hdf5') as f:
    
    ascc = f['table_header/ascc'].value
    ra = f['table_header/ra'].value
    dec = f['table_header/dec'].value
    
    print np.amin(ra), np.amax(ra)
    print np.amin(dec), np.amax(dec)
    
    decidx = pg.find_decidx(dec)
    
    #here = decidx == 495
    #ascc = ascc[here]
    #ra = ra[here]
    #dec = dec[here]
    
    here = ascc == '1413051'
    ascc = ascc[here]
    ra = ra[here]
    dec = dec[here]
    
    for i in range(0,len(ascc),100):
        
        lc = f['data/'+ascc[i]]
        
        ha = np.mod(lc['lst']*15.-ra[i], 360)
        
        decidx = pg.find_decidx(dec[i])
        haidx_cam = pg.find_raidx(ha) 
        haidx_ip = pg1.find_raidx(ha)
        
        fit = camtrans[decidx, haidx_cam]*(a[decidx, haidx_ip]*np.sin(2*np.pi*lc['y'])+b[decidx, haidx_ip]*np.cos(2*np.pi*lc['y'])+1)
        fit1 = camtrans[decidx, haidx_cam]
        fit2 = (a[decidx, haidx_ip]*np.sin(2*np.pi*lc['y'])+b[decidx, haidx_ip]*np.cos(2*np.pi*lc['y'])+1)
        
        plt.figure(figsize=(16,8))
        ax = plt.subplot(311)
        plt.title('ASCC %s'%ascc[i])
        plt.plot(lc['jdmid'], lc['flux0'], '.')
        plt.plot(lc['jdmid'], fit*np.median(lc['flux0']/fit), '.')
        plt.ylabel('flux0')
        ax1 = plt.subplot(312, sharex=ax)
        plt.plot(lc['jdmid'], lc['flux0']/(fit1*np.median(lc['flux0']/fit)), '.')
        plt.plot(lc['jdmid'], fit2, '.')
        plt.ylabel('camtrans+ip')
        plt.subplot(313, sharex=ax, sharey=ax1)
        plt.plot(lc['jdmid'], lc['flux0']/(fit*np.median(lc['flux0']/fit)), '.')
        plt.ylabel('cipflux0')
        plt.xlabel('Time [JD]')
        plt.tight_layout()
        plt.show()
    
    
    
    
