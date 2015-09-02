#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py

import matplotlib.pyplot as plt
from sysrem import sysrem
from index_functions import index_statistics
from coordinate_grids import PolarGrid

from matplotlib import rcParams
from viridis import viridis 

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

fLCfile = '/data2/talens/Jul2015/fLC_20150716LPC.hdf5'

with h5py.File(fLCfile) as f:
    ra = f['table_header/ra'].value
    dec = f['table_header/dec'].value
    ascc = f['table_header/ascc'].value
    nobs = f['table_header/nobs'].value.astype('int')
    
pg = PolarGrid(13500, 720)
decidx, decuni = pg.find_decidx(dec, compact=True)
print decidx[82]
ascc = ascc[decuni==82]
nobs = nobs[decuni==82]
ra = ra[decuni==82]

lst = np.array([])
flux0 = np.array([])
eflux0 = np.array([])

with h5py.File(fLCfile) as f:
    for sid in ascc:
        lc = f['data/'+sid]
        
        lst = np.append(lst, lc['lst'])
        flux0 = np.append(flux0, lc['flux0'])
        eflux0 = np.append(eflux0, lc['eflux0'])

staridx = np.repeat(np.arange(len(ascc)), nobs)

ha = np.mod(lst*15.-np.repeat(ra,nobs), 360.)
haidx, hauni = pg.find_raidx(ha, compact=True)

values = np.copy(flux0)
flux = np.zeros((5, len(ascc)))
trans = np.zeros((5, len(haidx)))

for i in range(5):

    flux[i], trans[i], niter, chisq, chisq1, chisq2, npoints, npars = sysrem(staridx, hauni, values, eflux0, maxiter=250)
    values -= flux[i, staridx]*trans[i, hauni]

    plt.figure(figsize=(16,8))
    plt.title('Component %i'%i)
    plt.plot(trans[i], '.')
    #plt.show()
    plt.savefig('component_%i.png'%i)

np.savez('/data2/talens/Jul2015/ip_test.npz', haidx=haidx, trans=trans)
  
result = np.dot(flux.T, trans)

flux_n = np.zeros((5, len(ascc)))
with h5py.File('/data2/talens/Jul2015/red_20150716LPC.hdf5') as f:

    for i in range(len(ascc)):
        
        rc = f['data/'+ascc[i]].value

        y = flux0[staridx==i]
        yerr = eflux0[staridx==i]
        comp = trans.T[hauni[staridx==i]]
        fit1 = result[i,hauni[staridx==i]]

        pars = np.linalg.lstsq(comp/yerr[:,None], y/yerr)[0]
        flux_n[:,i] = pars
        fit2 = np.dot(comp, pars)

        fit3 = rc['camtrans0']
        pars = np.sum(y*fit3/yerr**2)/np.sum(fit3**2/yerr**2)
        fit3 = fit3*pars

        fig = plt.figure(figsize=(16,8))
        ax = plt.subplot(311)
        plt.plot(y, '.')
        plt.plot(fit1, '.', c='g')
        plt.plot(fit2, '.', c='r')
        plt.plot(fit3, '.', c='k')
        ax2 = plt.subplot(312, sharex=ax)
        plt.plot(y/fit3, '.', c='k')
        plt.plot(y/fit1, '.', c='g')
        plt.subplot(313, sharex=ax, sharey=ax2)
        plt.plot(y/fit3, '.', c='k')
        plt.plot(y/fit2, '.', c='r')
        plt.ylim(.7,1.3)
        
        plt.tight_layout()
        plt.show()
        plt.close()
