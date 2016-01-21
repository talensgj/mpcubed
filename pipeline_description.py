#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from package import IO
from package.coordinates import grids
from package import misc

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

f = IO.fLCfile('/data2/talens/2015Q2/LPE/fLC_201506ALPE.hdf5')

ascc, ra, dec, nobs, vmag = f.read_header(['ascc', 'ra', 'dec', 'nobs', 'vmag'])
nobs = nobs.astype('int')

staridx = np.arange(len(ascc))

pg = grids.PolarGrid(13500, 36)
raidx, decidx = pg.radec2idx(ra, dec)

arg, = np.where(ascc == '807144')
select = (decidx == decidx[arg])

flux0, lst, lstseq, x, y = f.read_data(['flux0', 'lst', 'lstseq', 'x', 'y'], ascc[select], nobs[select])
lstseq = lstseq.astype('int')

mag0 = misc.flux2mag(flux0)

sys = IO.SysFile('/data2/talens/2015Q2/LPE/sys0_201506ALPE.hdf5')

sID, vmag, mag, _, _ = sys.read_magnitudes()
pgcam, trans, _ = sys.read_trans()
pgipx, sinx, cosx, siny, cosy, _ = sys.read_intrapix()
hg, clouds, _, _, lstmin, lstmax = sys.read_clouds()

ra_ = np.repeat(ra[select], nobs[select])
ha = np.mod(lst*15. - ra_, 360)
dec_ = np.repeat(dec[select], nobs[select])

staridx = np.repeat(staridx[select], nobs[select])
mag0 = mag0 - mag[staridx]

idx1, idx2 = pgcam.radec2idx(ha, dec_)
trans = trans[idx1, idx2]

idx = hg.radec2idx(ra_, dec_)
clouds = clouds[idx, lstseq - lstmin]

idx1, idx2 = pgipx.radec2idx(ha, dec_)
ipx_x = sinx[idx1, idx2]*np.sin(2*np.pi*x) + cosx[idx1, idx2]*np.cos(2*np.pi*x) 
ipx_y = siny[idx1, idx2]*np.sin(2*np.pi*y) + cosy[idx1, idx2]*np.cos(2*np.pi*y)
ipx = ipx_x + ipx_y

staridx, idx1 = np.unique(staridx, return_inverse=True)
nstars = len(staridx)

lstseq, idx2 = np.unique(lstseq, return_inverse=True)
ntimes = len(lstseq)

# Put the dat in 2d arrays.
tmp = np.full((nstars, ntimes), fill_value=np.nan)
tmp[idx1, idx2] = mag0
mag0 = tmp

tmp = np.full((nstars, ntimes), fill_value=np.nan)
tmp[idx1, idx2] = trans
trans = tmp

tmp = np.full((nstars, ntimes), fill_value=np.nan)
tmp[idx1, idx2] = clouds
clouds = tmp

tmp = np.full((nstars, ntimes), fill_value=np.nan)
tmp[idx1, idx2] = ipx
ipx = tmp

# Plot the data and the transmission.
vmin = np.nanpercentile(mag0, 1)
vmax = np.nanpercentile(mag0, 99)

ax = plt.subplot(211)

plt.imshow(mag0, aspect='auto', vmin=vmin, vmax=vmax)
plt.colorbar()

args, = np.where(np.diff(lstseq)>1)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k', ls='--')
    
args, = np.where(np.diff(lstseq//13500) > 0)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k')

plt.subplot(212, sharex=ax, sharey=ax)

plt.imshow(trans, aspect='auto', vmin=vmin, vmax=vmax)

args, = np.where(np.diff(lstseq)>1)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k', ls='--')
    
args, = np.where(np.diff(lstseq//13500) > 0)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k')
    
plt.show()

# Plot the data and the clouds.
vmin = np.nanpercentile(mag0 - trans, 1)
vmax = np.nanpercentile(mag0 - trans, 99)

ax = plt.subplot(211)

plt.imshow(mag0 - trans, aspect='auto', vmin=vmin, vmax=vmax)
plt.colorbar()

args, = np.where(np.diff(lstseq)>1)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k', ls='--')
    
args, = np.where(np.diff(lstseq//13500) > 0)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k')

plt.subplot(212, sharex=ax, sharey=ax)

plt.imshow(clouds, aspect='auto', vmin=vmin, vmax=vmax)

args, = np.where(np.diff(lstseq)>1)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k', ls='--')
    
args, = np.where(np.diff(lstseq//13500) > 0)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k')
    
plt.show()

# Plot the data and the ipx.
vmin = np.nanpercentile(mag0 - trans - clouds, 1)
vmax = np.nanpercentile(mag0 - trans - clouds, 99)

ax = plt.subplot(211)

plt.imshow(mag0 - trans - clouds, aspect='auto', vmin=vmin, vmax=vmax)
plt.colorbar()

args, = np.where(np.diff(lstseq)>1)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k', ls='--')
    
args, = np.where(np.diff(lstseq//13500) > 0)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k')
    
plt.subplot(212, sharex=ax, sharey=ax)

plt.imshow(ipx, aspect='auto', vmin=vmin, vmax=vmax)

args, = np.where(np.diff(lstseq)>1)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k', ls='--')
    
args, = np.where(np.diff(lstseq//13500) > 0)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k')
    
plt.show()

# Plot the reduced data.
vmin = np.nanpercentile(mag0 - trans - clouds - ipx, 1)
vmax = np.nanpercentile(mag0 - trans - clouds - ipx, 99)

plt.imshow(mag0 - trans - clouds, aspect='auto', vmin=vmin, vmax=vmax)
plt.colorbar()

args, = np.where(np.diff(lstseq)>1)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k', ls='--')
    
args, = np.where(np.diff(lstseq//13500) > 0)
for arg in args:
    print arg
    plt.axvline(arg+.5, c='k')
    
plt.show()
