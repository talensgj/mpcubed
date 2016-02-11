#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from package import IO
from package import plotting
from package.systematics import cdecor
from package.coordinates import grids
from package.red_apply import CorrectLC
from package import misc

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

def time2sexa(time):
    
    hours = time//1.
    minutes = (time - hours)*60.//1
    seconds = (time - hours - minutes/60.)*3600.
    
    return hours, minutes, seconds
    
def angle2sexa(angle):
    
    degrees = angle//1.
    arcmin = (angle - degrees)*60//1.
    arcsec = (angle - degrees - arcmin/60.)*3600
    
    return degrees, arcmin, arcsec
    
# Read the header.
f = IO.fLCfile('/data2/talens/fLC/fLC_20150605LPE.hdf5')
ascc, ra, dec, vmag, nobs, jdstart = f.read_header(['ascc', 'ra', 'dec', 'vmag', 'nobs', 'jdstart'])
nobs = nobs.astype('int')

# Select a subset of all stars.
hg = grids.HealpixGrid(8)
skyidx = hg.radec2idx(ra, dec)

arg, = np.where(ascc == '807144')
rabin, decbin = hg.idx2radec(skyidx[arg])

select = (skyidx == skyidx[arg])
ascc = ascc[select]
ra = ra[select]
dec = dec[select]
vmag = vmag[select]
nobs = nobs[select]
jdstart = jdstart[select]

# Read the lightcurves.
lstseq, lstidx, lst, jdmid, flux0, eflux0, x, y = f.read_data(['lstseq', 'lstidx', 'lst', 'jdmid', 'flux0', 'eflux0', 'x', 'y'], ascc, nobs)
lstseq = lstseq.astype('int')
lstidx = lstidx.astype('int')
mag0, emag0 = misc.flux2mag(flux0, eflux0)

print np.amin(lstseq), np.amax(lstseq)

# Get the systematic corrections.
sys = IO.SysFile('/data2/talens/2015Q2/LPE/sys0_201506ALPE.hdf5')

# Magnitudes.
ascc_mag, vmag, mag, sigma, _ = sys.read_magnitudes()
tmp = []
for sid in ascc:
    here = ascc_mag == sid
    tmp.append(mag[here])
mag = np.array(tmp)

# Transmission.
pg, trans, _ = sys.read_trans()
ra = np.repeat(ra, nobs)
dec = np.repeat(dec, nobs)
ha = np.mod(lst*15 - ra, 360)
haidx, decidx = pg.radec2idx(ha, dec)
trans = trans[haidx, decidx]

# Intrapixel.
pg, sinx, cosx, siny, cosy, _ = sys.read_intrapix()
haidx, decidx = pg.radec2idx(ha, dec)
ipx_x = sinx[haidx, decidx]*np.sin(2*np.pi*x) + cosx[haidx, decidx]*np.cos(2*np.pi*x)
ipx_y = siny[haidx, decidx]*np.sin(2*np.pi*y) + cosy[haidx, decidx]*np.cos(2*np.pi*y)
ipx = ipx_x + ipx_y

# Clouds.
hg, clouds, sigma, _, lstmin, lstmax = sys.read_clouds()
clouds = clouds[skyidx[arg], lstseq - lstmin]

# Plot as 2d array.
nstars = len(ascc)
staridx = np.repeat(np.arange(nstars), nobs)

image1 = np.full((nstars, 13500), fill_value = np.nan)
image1[staridx, lstidx] = mag0 - np.repeat(mag, nobs) - trans - ipx
image1 = image1[np.argsort(jdstart)]

image2 = np.full((nstars, 13500), fill_value = np.nan)
image2[staridx, lstidx] = mag0 - np.repeat(mag, nobs) - trans - ipx - clouds
image2 = image2[np.argsort(jdstart)]

fig = plt.figure(figsize=(16,6))

plt.suptitle('2015-06-05 East', size='xx-large')

gs = gridspec.GridSpec(3, 2, width_ratios = [15,.2], height_ratios = [1,10,10])

vmin = np.nanpercentile(image1, 1)
vmax = np.nanpercentile(image1, 99)
delta = vmax - vmin
vmin = vmin - .01*delta
vmax = vmax + .01*delta

plt.subplot(gs[1,0], yticks=[], xticks=[])
plt.annotate(r'$(\alpha, \delta) = ({:.1f}^\circ, {:.1f}^\circ)$'.format(rabin[0], decbin[0]), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='x-large', backgroundcolor='w')
im = plt.imshow(image1, aspect='auto', vmin=vmin, vmax=vmax, cmap=cm.Greys)
plt.xlim(np.amin(lstidx)-.5, np.amax(lstidx)-.5)
#plt.xlabel('Time')

ax = plt.subplot(gs[1,1])
cb = plt.colorbar(im, cax=ax)
cb.ax.invert_yaxis()
cb.set_label('Magnitude')

vmin = np.nanpercentile(image2, 1)
vmax = np.nanpercentile(image2, 99)
delta = vmax - vmin
vmin = vmin - .01*delta
vmax = vmax + .01*delta

plt.subplot(gs[2,0], yticks=[])
im = plt.imshow(image2, aspect='auto', vmin=vmin, vmax=vmax, cmap=cm.Greys)
plt.xlim(np.amin(lstidx)-.5, np.amax(lstidx)-.5)
plt.xlabel('Time')

ax = plt.subplot(gs[2,1])
cb = plt.colorbar(im, cax=ax)
cb.ax.invert_yaxis()
cb.set_label('Magnitude')

plt.tight_layout()
plt.show()
