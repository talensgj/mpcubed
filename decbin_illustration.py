#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

from coordinate_grids import PolarGrid
from index_functions import index_statistics

import sysrem_dev

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis 
from matplotlib.ticker import FuncFormatter, MultipleLocator, NullLocator
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

def periodic(x, pos):
        'The two args are the value and tick position'
        x = x%24
        x1 = np.floor(x)
        x2 = (x-np.floor(x))*60
        
        return '%02.fh%02.fm' % (x1, x2)

formatter = FuncFormatter(periodic)

pg1 = PolarGrid(13500, 720)
pg2 = PolarGrid(270, 720)

with h5py.File('/data2/talens/Jul2015/fLC_20150716LPC.hdf5', 'r') as f:
    
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
    
    lst = np.array([])
    y = np.array([])
    flux = np.array([])
    eflux = np.array([])
    sflux = np.array([])
    for i in range(len(ascc)):
        lc = f['data/'+ascc[i]]
        
        lst = np.append(lst, lc['lst'])
        y = np.append(y, lc['y'])
        flux = np.append(flux, lc['flux0'])
        eflux = np.append(eflux, lc['eflux0'])
        tmp = index_statistics(lc['lstidx']//50, lc['flux0'], statistic='std', keeplength=True)
        sflux = np.append(sflux, tmp)
    
eflux = sflux

ha = np.mod(lst*15.-np.repeat(ra, nobs), 360)

haidx1 = pg1.find_raidx(ha)
haidx2 = pg2.find_raidx(ha)
staridx = np.repeat(np.arange(len(ascc)), nobs)

here = (flux > 0)&(eflux > 0)
staridx = staridx[here]
haidx1 = haidx1[here]
haidx2 = haidx2[here]
y = y[here]
flux = flux[here]
eflux = eflux[here]

data = np.full((len(ascc), 13502), fill_value=np.nan)
data[staridx, haidx1] = flux
data = data[:,1:-1]
data = np.roll(data, 13500/2, axis=1)
data = data/np.nanmedian(data, axis=1, keepdims=True)

sort = np.argsort(ra)
data = data[sort]

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_locator(NullLocator())
plt.title(r'16-07-2015 LPC, Dec = 22.50$^\circ$ - 22.75$^\circ$')
plt.imshow(data, aspect='auto', cmap=viridis, extent=(-12,12,0,1), vmin=0, vmax=2.5)
cb = plt.colorbar()
plt.xlim(-3,3)
cb.set_label('Normalized Flux')
plt.xlabel('Hour Angle')
plt.tight_layout()
#plt.show()
plt.savefig('Pipeline_part1_0_data.png')

F, T, niter, chisq, npoints, npars = sysrem_dev.sysrem(staridx, haidx1, flux, eflux, verbose=True)

x = np.arange(len(T))
x = pg1.find_ra(x)
x = np.mod(x-180, 360)-180

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
plt.title(r'$\rm{niter}=%i$, $\chi^{2}=%.2f$'%(niter, chisq))
plt.plot(x/15, T, '.')
plt.xlim(-3,3)
plt.ylim(0, 2.5)
plt.xlabel('Hour Angle')
plt.ylabel('Transmission')
#plt.show()
plt.savefig('Pipeline_part1_1_transmission1.png')

data = np.full((len(ascc), 13502), fill_value=np.nan)
data[staridx, haidx1] = flux/(F[staridx]*T[haidx1])
data = data[:,1:-1]
data = np.roll(data, 13500/2, axis=1)
data = data/np.nanmedian(data, axis=1, keepdims=True)

sort = np.argsort(ra)
data = data[sort]

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_locator(NullLocator())
plt.title(r'16-07-2015 LPC, Dec = 22.50$^\circ$ - 22.75$^\circ$')
plt.imshow(data, aspect='auto', cmap=viridis, extent=(-12,12,0,1), vmin=.8, vmax=1.2)
cb = plt.colorbar()
plt.xlim(-3,3)
cb.set_label('Reduced Flux')
plt.xlabel('Hour Angle')
#plt.show()
plt.savefig('Pipeline_part1_2_result.png')

F, T, a, b, niter, chisq, npoints, npars = sysrem_dev.intrarem(staridx, haidx1, haidx2, y, flux, eflux, verbose=True)

x = np.arange(len(T))
x = pg1.find_ra(x)
x = np.mod(x-180, 360)-180

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
plt.title(r'$\rm{niter}=%i$, $\chi^{2}=%.2f$'%(niter, chisq))
plt.plot(x/15, T, '.')
plt.xlim(-3,3)
plt.ylim(0,2.5)
plt.xlabel('Hour Angle')
plt.ylabel('Transmission')
#plt.show()
plt.savefig('Pipeline_part1_3_transmission2.png')

x = np.arange(len(a))
x = pg2.find_ra(x)
x = np.mod(x-180, 360)-180

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
plt.title(r'$\rm{niter}=%i$, $\chi^{2}=%.2f$'%(niter, chisq))
plt.plot(x/15, a, '.', label='a')
plt.plot(x/15, b, '.', label='b')
plt.xlim(-3,3)
plt.ylim(-.06, .06)
plt.xlabel('Hour Angle')
plt.ylabel('Amplitude')
#plt.show()
plt.savefig('Pipeline_part1_4_amplitudes.png')

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
plt.title(r'$\rm{niter}=%i$, $\chi^{2}=%.2f$'%(niter, chisq))
plt.plot(x/15, np.sqrt(a**2+b**2), '.', label='a')
plt.xlim(-3,3)
plt.ylim(0, .06)
plt.xlabel('Hour Angle')
plt.ylabel('Amplitude')
#plt.show()
plt.savefig('Pipeline_part1_5_amplitude.png')

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
plt.title(r'$\rm{niter}=%i$, $\chi^{2}=%.2f$'%(niter, chisq))
plt.plot(x/15, np.arctan(b/a), '.', label='a')
plt.xlim(-3,3)
plt.ylim(-np.pi/2, np.pi/2)
plt.xlabel('Hour Angle')
plt.ylabel('Phase')
#plt.show()
plt.savefig('Pipeline_part1_6_phase.png')

data = np.full((len(ascc), 13502), fill_value=np.nan)
data[staridx, haidx1] = (a[haidx2]*np.sin(2*np.pi*y)+b[haidx2]*np.cos(2*np.pi*y)+1)
data = data[:,1:-1]
data = np.roll(data, 13500/2, axis=1)
data = data/np.nanmedian(data, axis=1, keepdims=True)

sort = np.argsort(ra)
data = data[sort]

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_locator(NullLocator())
plt.title(r'$\rm{niter}=%i$, $\chi^{2}=%.2f$'%(niter, chisq))
plt.imshow(data, aspect='auto', cmap=viridis, extent=(-12,12,0,1), vmin=.8, vmax=1.2)
cb = plt.colorbar()
plt.xlim(-3,3)
cb.set_label('Intrapixel Variations')
plt.xlabel('Hour Angle')
plt.tight_layout()
#plt.show()
plt.savefig('Pipeline_part1_7_intrapixel.png')

data = np.full((len(ascc), 13502), fill_value=np.nan)
data[staridx, haidx1] = flux/(F[staridx]*T[haidx1]*(a[haidx2]*np.sin(2*np.pi*y)+b[haidx2]*np.cos(2*np.pi*y)+1))
data = data[:,1:-1]
data = np.roll(data, 13500/2, axis=1)
data = data/np.nanmedian(data, axis=1, keepdims=True)

sort = np.argsort(ra)
data = data[sort]

plt.figure(figsize=(12,6))
ax = plt.subplot(111)
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_locator(NullLocator())
plt.title(r'16-07-2015 LPC, Dec = 22.50$^\circ$ - 22.75$^\circ$')
plt.imshow(data, aspect='auto', cmap=viridis, extent=(-12,12,0,1), vmin=.8, vmax=1.2)
cb = plt.colorbar()
plt.xlim(-3,3)
cb.set_label('Reduced Flux')
plt.xlabel('Hour Angle')
plt.tight_layout()
#plt.show()
plt.savefig('Pipeline_part1_8_result2.png')
