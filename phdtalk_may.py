#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
rcParams['xtick.labelsize'] = 'medium'
rcParams['ytick.labelsize'] = 'medium'
rcParams['axes.labelsize'] = 'large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package import IO
from package import misc
from package import plotting
from package.red_apply_vmag import CorrectLC

import fourier
import detrend
from pea_plot import plot_Polar
from transit_candidates import RedFile, find_ns

#flcfile = '/data2/mascara/LaPalma/20150619LPC/fLC/fLC_20150619LPC.hdf5'
#sysfile = '/data3/talens/2015Q2/LPC/sys0_vmag_201506BLPC.hdf5'
#redfile = '/data3/talens/2015Q2/LPC/red0_vmag_2015Q2LPC.hdf5'

flcfile = '/data2/talens/LaPalma/fLC/fLC_20150619LPC.hdf5'
sysfile = '/data2/talens/2015Q2_vmag/LPC/sys0_vmag_201506BLPC.hdf5'
redfile = '/data2/talens/red2015/red0_vmag_2015Q2LPC.hdf5'

# Figure showing a raw lightcurve.
f = IO.fLCfile(flcfile)
ascc, ra, dec, vmag = f.read_header(['ascc', 'ra', 'dec', 'vmag'])
sel = (ascc == '807144')
ra = ra[sel]
dec = dec[sel]
vmag = vmag[sel]

lstseq, jdmid, lst, flux, eflux, x, y = f.read_data(['lstseq', 'jdmid', 'lst', 'flux0', 'eflux0', 'x', 'y'], ['807144'])
lstseq = lstseq.astype('int')

mag, emag = misc.flux2mag(flux, eflux)

jdref = np.floor(jdmid[0])

fig = plt.figure(figsize=(16,5))

gs = gridspec.GridSpec(2, 1, height_ratios = [1,20])

plt.suptitle('ASCC 807144, 2015-06-19 LPC', size='xx-large')

ax = plt.subplot(gs[1])
plt.plot(jdmid - jdref, mag, '.', c='k', alpha=.5)
ax.invert_yaxis()
plt.xlabel('Time [JD - {:.0f}]'.format(jdref))
plt.ylabel('Magnitude')

plt.tight_layout()
plt.show()
plt.close()

f = IO.SysFile(sysfile)
pg, trans, nobs = f.read_trans()

ha = np.mod(lst*15-ra, 360)
idx1, idx2 = pg.radec2idx(ha, dec)
trans = trans[idx1, idx2]

pg, a, b, c, d, nobs = f.read_intrapix()

idx1, idx2 = pg.radec2idx(ha, dec)
ipx_x = a[idx1, idx2]*np.sin(2*np.pi*x) + b[idx1, idx2]*np.cos(2*np.pi*x)
ipx_y = c[idx1, idx2]*np.sin(2*np.pi*y) + d[idx1, idx2]*np.cos(2*np.pi*y)
ipx = ipx_x + ipx_y

hg, clouds, sigma, nobs, lstmin, lstmax = f.read_clouds()

idx = hg.radec2idx(ra, dec)
clouds = clouds[idx, lstseq - lstmin]

idx = lstseq//50

m0 = np.bincount(idx)
mt = np.bincount(idx, jdmid)
m1 = np.bincount(idx, mag - trans - ipx - clouds)
m2 = np.bincount(idx, (mag - trans - ipx - clouds)**2)

jdmid_bin = mt/m0
mag_bin = m1/m0
emag_bin = np.sqrt(((m2/m0) - (m1/m0)**2))

fig = plt.figure(figsize=(16,5))

gs = gridspec.GridSpec(2, 1, height_ratios = [1,20])

plt.suptitle('ASCC 807144, 2015-06-19 LPC', size='xx-large')

ax = plt.subplot(gs[1])
plt.plot(jdmid - jdref, mag - clouds, '.', c='k', alpha=.5)
ax.invert_yaxis()
plt.xlabel('Time [JD - {:.0f}]'.format(jdref))
plt.ylabel('Magnitude')

plt.tight_layout()
plt.show()
plt.close()

fig = plt.figure(figsize=(16,5))

gs = gridspec.GridSpec(2, 1, height_ratios = [1,20])

plt.suptitle('ASCC 807144, 2015-06-19 LPC', size='xx-large')

ax = plt.subplot(gs[1])
plt.plot(jdmid - jdref, mag - clouds, '.', c='k', alpha=.5)
plt.plot(jdmid - jdref, vmag + trans + ipx, c='r', lw=2)
ax.invert_yaxis()
plt.xlabel('Time [JD - {:.0f}]'.format(jdref))
plt.ylabel('Magnitude')

plt.tight_layout()
plt.show()
plt.close()

fig = plt.figure(figsize=(16,5))

gs = gridspec.GridSpec(2, 1, height_ratios = [1,20])

plt.suptitle('ASCC 807144, 2015-06-19 LPC', size='xx-large')

ax = plt.subplot(gs[1])
plt.plot(jdmid - jdref, mag - trans - ipx, '.', c='k', alpha=.5)
ax.invert_yaxis()
plt.xlabel('Time [JD - {:.0f}]'.format(jdref))
plt.ylabel('Magnitude')

plt.tight_layout()
plt.show()
plt.close()

fig = plt.figure(figsize=(16,5))

gs = gridspec.GridSpec(2, 1, height_ratios = [1,20])

plt.suptitle('ASCC 807144, 2015-06-19 LPC', size='xx-large')

ax = plt.subplot(gs[1])
plt.plot(jdmid - jdref, mag - trans - ipx, '.', c='k', alpha=.5)
plt.plot(jdmid - jdref, vmag + clouds, c='r', lw=2)
ax.invert_yaxis()
plt.xlabel('Time [JD - {:.0f}]'.format(jdref))
plt.ylabel('Magnitude')

plt.tight_layout()
plt.show()
plt.close()

# Figure showing LST correction.
f = RedFile(redfile)
lc = f.read_lightcurve('807144', ['lstseq', 'jdmid', 'lst', 'mag0', 'emag0', 'nobs'])
lstseq = lc['lstseq']
jdmid = lc['jdmid']
lst = lc['lst']
mag = lc['mag0']
emag = lc['emag0']
nobs = lc['nobs']

sel = (nobs == 50)
lstseq = lstseq[sel]
jdmid = jdmid[sel]
lst = lst[sel]
mag = mag[sel]
emag = emag[sel]

n1 = np.ptp(lstseq)
n2, wrap = find_ns(lstseq)

if wrap:
    lst = np.mod(lst+12, 24)-12

n1 = np.maximum(n1, 2)
n2 = np.maximum(n2, 2)

weights = 1/(emag)**2
freq1, freq2, pars1, pars2, trend1, trend2, chisq = detrend.new_harmonic3(jdmid, lst, mag, weights, [n1, n2])

fig = plt.figure(figsize=(16,5))

gs = gridspec.GridSpec(2, 2, height_ratios = [1,20], width_ratios=[1,2])

plt.suptitle('ASCC 807144, 2015Q2 LPC', size='xx-large')

ax = plt.subplot(gs[1, 0])
plt.plot(lst, mag - trend1, '.', c='k', alpha=.5)
ax.invert_yaxis()
plt.xlabel('Time [LST]')
plt.ylabel('Magnitude')

ax = plt.subplot(gs[1, 1])
plt.plot(jdmid - jdref, mag - trend2, '.', c='k', alpha=.5)
ax.invert_yaxis()
plt.xlabel('Time [JD - {:.0f}]'.format(jdref))
plt.ylabel('Magnitude')

plt.tight_layout()
plt.show()
plt.close()

fig = plt.figure(figsize=(16,5))

x1 = np.linspace(np.amin(jdmid), np.amax(jdmid), 1000)
mat1 = fourier.fourier_mat(x1, freq1)
y1 = np.dot(mat1, pars1)

x2 = np.linspace(np.amin(lst), np.amax(lst), 100)
mat2 = fourier.fourier_mat(x2, freq2)
y2 = np.dot(mat2, pars2)

gs = gridspec.GridSpec(2, 2, height_ratios = [1,20], width_ratios=[1,2])

plt.suptitle('ASCC 807144, 2015Q2 LPC', size='xx-large')

ax = plt.subplot(gs[1, 0])
plt.plot(lst, mag - trend1, '.', c='k', alpha=.5)
plt.plot(x2, y2, c='r', lw=2)
ax.invert_yaxis()
plt.xlabel('Time [LST]')
plt.ylabel('Magnitude')

ax = plt.subplot(gs[1, 1])
plt.plot(jdmid - jdref, mag - trend2, '.', c='k', alpha=.5)
plt.plot(x1 - jdref, y1, c='r', lw=2)
ax.invert_yaxis()
plt.xlabel('Time [JD - {:.0f}]'.format(jdref))
plt.ylabel('Magnitude')

plt.tight_layout()
plt.show()
plt.close()

exit()
# Figure showing the transmission map.
f = IO.SysFile('/data2/talens/2015Q2_vmag/LPC/sys0_vmag_201506ALPC.hdf5')
pg, trans, nobs = f.read_trans()

vmin = np.nanpercentile(trans, 1)
vmax = np.nanpercentile(trans, 99)

majorLocator1 = MultipleLocator(668)
majorLocator2 = MultipleLocator(668)

fig = plt.figure(figsize=(15,10))

gs = gridspec.GridSpec(2, 2, height_ratios = [1,20], width_ratios = [30,1])

plt.suptitle('2015-06 "B" LPC', size='xx-large')

ax = plt.subplot(gs[1,0], aspect='equal')
im = plot_Polar(pg, trans, cmap=plotting.viridis, vmin=vmin, vmax=vmax)
plt.xlim(0, 4008)
plt.ylim(0, 2672)
ax.xaxis.set_major_locator(majorLocator1)
ax.yaxis.set_major_locator(majorLocator2)

ax = plt.subplot(gs[1,1])
cbar = plt.colorbar(im, cax=ax)
cbar.set_label('Transmission')

plt.tight_layout()
plt.show()
plt.close()

# Figure showing the intrapixel amplitudes.
f = IO.SysFile('/data2/talens/2015Q2_vmag/LPC/sys0_vmag_201506BLPC.hdf5')
pg, sinx, cosx, siny, cosy, nobs = f.read_intrapix()

vmin = -.05
vmax = .05

majorLocator1 = MultipleLocator(668)
majorLocator2 = MultipleLocator(668)

fig = plt.figure(figsize=(15,10))

gs = gridspec.GridSpec(3, 3, height_ratios = [1,10,10], width_ratios = [15,15,1])

plt.suptitle('2015-06 "B" LPC', size='xx-large')

ax = plt.subplot(gs[1,0], aspect='equal')
im = plot_Polar(pg, sinx, cmap=plotting.viridis, vmin=vmin, vmax=vmax)
plt.xlim(0, 4008)
plt.ylim(0, 2672)
ax.xaxis.set_major_locator(majorLocator1)
ax.yaxis.set_major_locator(majorLocator2)

ax = plt.subplot(gs[1,1], aspect='equal')
im = plot_Polar(pg, cosx, cmap=plotting.viridis, vmin=vmin, vmax=vmax)
plt.xlim(0, 4008)
plt.ylim(0, 2672)
ax.xaxis.set_major_locator(majorLocator1)
ax.yaxis.set_major_locator(majorLocator2)

ax = plt.subplot(gs[2,0], aspect='equal')
im = plot_Polar(pg, siny, cmap=plotting.viridis, vmin=vmin, vmax=vmax)
plt.xlim(0, 4008)
plt.ylim(0, 2672)
ax.xaxis.set_major_locator(majorLocator1)
ax.yaxis.set_major_locator(majorLocator2)

ax = plt.subplot(gs[2,1], aspect='equal')
im = plot_Polar(pg, cosy, cmap=plotting.viridis, vmin=vmin, vmax=vmax)
plt.xlim(0, 4008)
plt.ylim(0, 2672)
ax.xaxis.set_major_locator(majorLocator1)
ax.yaxis.set_major_locator(majorLocator2)

ax = plt.subplot(gs[1:3,2])
cbar = plt.colorbar(im, cax=ax)
cbar.set_label('Amplitude')

plt.tight_layout()
plt.show()
plt.close()

# Figure showing the lightcurve.
f = IO.fLCfile('/data2/talens/2015Q2_pea/LPC/fLC_201506BLPC.hdf5')
lstseq, lst, jdmid, flux, eflux = f.read_data(['lstseq', 'lst', 'jdmid', 'flux0', 'eflux0'], ['807144'])

jdmid = jdmid - 2400000.5
mag, emag = misc.flux2mag(flux, eflux)

f = CorrectLC('/data2/talens/2015Q2_pea/LPC/fLC_201506BLPC.hdf5', 0, '/data2/talens/2015Q2_vmag/LPC/sys0_vmag_201506BLPC.hdf5')
trans, intrapix, clouds, flags = f.get_correction('807144')

lstseq = lstseq.astype('int')

lstday = lstseq//13500
lstidx = lstseq%13500

lstday = lstday - np.amin(lstday)
ndays = np.amax(lstday) + 1

tmp = np.full((13500, ndays), fill_value=np.nan)
tmp[lstidx, lstday] = mag
y = tmp

tmp = np.full((13500, ndays), fill_value=np.nan)
tmp[lstidx, lstday] = mag - trans - intrapix - clouds
res = tmp

vmin = np.nanpercentile(y, 1)
vmax = np.nanpercentile(y, 99)

plt.figure(figsize=(16,8))

gs = gridspec.GridSpec(3, 2, height_ratios = [1,10,10], width_ratios = [30,1])

plt.suptitle('ASCC 807144, 2015-06 "B" LPC', size='xx-large')

ax1 = plt.subplot(gs[1,0], aspect='auto')
im1 = plt.imshow(y.T, aspect='auto', extent=(0, 24, -.5, ndays-.5), cmap=plotting.viridis, vmin=vmin, vmax=vmax)
plt.xlim(17, 23)
plt.ylabel('Sidereal Day')

cax = plt.subplot(gs[1,1])
cbar = plt.colorbar(im1, cax=cax)
cbar.set_label('$\Delta m$')

ax2 = plt.subplot(gs[2,0], aspect='auto')
im2 = plt.imshow(res.T, aspect='auto', extent=(0, 24, -.5, ndays-.5), cmap=plotting.viridis, vmin=-.5, vmax=.5)
plt.xlim(17, 23)
plt.xlabel('Time [LST]')
plt.ylabel('Sidereal Day')

cax = plt.subplot(gs[2,1])
cbar = plt.colorbar(im2, cax=cax)
cbar.set_label('$\Delta m$')

plt.tight_layout()
plt.show()

tmp = np.full((13500, ndays), fill_value=np.nan)
tmp[lstidx, lstday] = mag
y = tmp

tmp = np.full((13500, ndays), fill_value=np.nan)
tmp[lstidx, lstday] = trans
res = tmp

vmin = np.nanpercentile(y, 1)
vmax = np.nanpercentile(y, 99)

plt.figure(figsize=(16,8))

gs = gridspec.GridSpec(3, 2, height_ratios = [1,10,10], width_ratios = [30,1])

plt.suptitle('ASCC 807144, 2015-06 "B" LPC', size='xx-large')

ax1 = plt.subplot(gs[1,0], aspect='auto')
im1 = plt.imshow(y.T, aspect='auto', extent=(0, 24, -.5, ndays-.5), cmap=plotting.viridis, vmin=vmin, vmax=vmax)
plt.xlim(17, 23)
plt.ylabel('Sidereal Day')

ax2 = plt.subplot(gs[2,0], aspect='auto')
im2 = plt.imshow(res.T, aspect='auto', extent=(0, 24, -.5, ndays-.5), cmap=plotting.viridis, vmin=vmin, vmax=vmax)
plt.xlim(17, 23)
plt.xlabel('Time [LST]')
plt.ylabel('Sidereal Day')

cax = plt.subplot(gs[1:3,1])
cbar = plt.colorbar(im1, cax=cax)
cbar.set_label('$\Delta m$')

plt.tight_layout()
plt.show()

tmp = np.full((13500, ndays), fill_value=np.nan)
tmp[lstidx, lstday] = mag - trans
y = tmp

tmp = np.full((13500, ndays), fill_value=np.nan)
tmp[lstidx, lstday] = clouds
res = tmp

vmin = np.nanpercentile(y, 1)
vmax = np.nanpercentile(y, 99)

plt.figure(figsize=(16,8))

gs = gridspec.GridSpec(3, 2, height_ratios = [1,10,10], width_ratios = [30,1])

plt.suptitle('ASCC 807144, 2015-06 "B" LPC', size='xx-large')

ax1 = plt.subplot(gs[1,0], aspect='auto')
im1 = plt.imshow(y.T, aspect='auto', extent=(0, 24, -.5, ndays-.5), cmap=plotting.viridis, vmin=-.5, vmax=.5)
plt.xlim(17, 23)
plt.ylabel('Sidereal Day')

ax2 = plt.subplot(gs[2,0], aspect='auto')
im2 = plt.imshow(res.T, aspect='auto', extent=(0, 24, -.5, ndays-.5), cmap=plotting.viridis, vmin=-.5, vmax=.5)
plt.xlim(17, 23)
plt.xlabel('Time [LST]')
plt.ylabel('Sidereal Day')

cax = plt.subplot(gs[1:3,1])
cbar = plt.colorbar(im1, cax=cax)
cbar.set_label('$\Delta m$')

plt.tight_layout()
plt.show()

tmp = np.full((13500, ndays), fill_value=np.nan)
tmp[lstidx, lstday] = mag - trans - clouds
y = tmp

tmp = np.full((13500, ndays), fill_value=np.nan)
tmp[lstidx, lstday] = intrapix
res = tmp

vmin = np.nanpercentile(y, 1)
vmax = np.nanpercentile(y, 99)

plt.figure(figsize=(16,8))

gs = gridspec.GridSpec(3, 2, height_ratios = [1,10,10], width_ratios = [30,1])

plt.suptitle('ASCC 807144, 2015-06 "B" LPC', size='xx-large')

ax1 = plt.subplot(gs[1,0], aspect='auto')
im1 = plt.imshow(y.T, aspect='auto', extent=(0, 24, -.5, ndays-.5), cmap=plotting.viridis, vmin=-.5, vmax=.5)
plt.xlim(17, 23)
plt.ylabel('Sidereal Day')

ax2 = plt.subplot(gs[2,0], aspect='auto')
im2 = plt.imshow(res.T, aspect='auto', extent=(0, 24, -.5, ndays-.5), cmap=plotting.viridis, vmin=-.5, vmax=.5)
plt.xlim(17, 23)
plt.xlabel('Time [LST]')
plt.ylabel('Sidereal Day')

cax = plt.subplot(gs[1:3,1])
cbar = plt.colorbar(im1, cax=cax)
cbar.set_label('$\Delta m$')

plt.tight_layout()
plt.show()
