#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams
from viridis import viridis

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['axes.titlesize'] = 'xx-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'

filelist = glob.glob('/data2/talens/3mEast/fLC_*LPE.hdf5')
filelist = np.sort(filelist)

ndays = len(filelist)
flux = np.full((ndays, 13500), fill_value=np.nan)

#['759956' '796928' '818447' '753845' '803519' '814935' '801095' '764977'
 #'823206' '823942' '747840' '748717' '760392' '802147' '828815' '742410'
 #'741539' '758439' '743525' '739767' '748332' '826202' '803473' '802365'
 #'744559' '818742' '764851' '757199' '753753' '754304' '816042' '795152'
 #'802987' '766876' '743279' '802837' '759515' '807144' '807784' '748128'
 #'749452' '825385' '739935' '799327' '751482' '742293' '761260' '815351'
 #'822579' '800733' '744842' '797908' '801953' '742860' '802815' '748480'
 #'753744' '764750' '805488' '763502' '805979' '749810' '750059' '802172'
 #'822178' '739407' '754670' '750837' '752649' '766412' '748382' '802943'
 #'819876' '752532' '762102' '808160' '818493' '755319' '803163' '748353']

ascc = '822579'
for i in range(ndays):
    with h5py.File(filelist[i]) as f:
        try: lc = f['data/'+ascc].value
        except: continue
        ra = f['header/'+ascc]['ra']
    
    lstidx = lc['lstidx'].astype('int')
    flux[i, lstidx] = lc['flux0']

mag = -2.5*np.log10(flux)
raidx = np.floor(ra/360*13500).astype('int')
mag = mag - mag[10]
mag = mag - mag[:,11500][:,None]
#mag = mag - np.nanmean(mag, axis=1, keepdims=True)

print raidx
mag = np.roll(mag, -raidx, axis=1)

ylim, xlim = np.where(np.isfinite(flux))

plt.figure(figsize=(16,8))
plt.title('ASCC '+ascc)
plt.imshow(mag, aspect='auto', cmap=viridis, vmin=-.5, vmax=.5, extent=(0,24,0,1))
plt.colorbar().set_label(r'$\Delta m$')
plt.xlabel('Hour Angle')
plt.xlim(18.75, 23)
#plt.ylim(np.amin(ylim)-.5, np.amax(ylim)+.5)
plt.tight_layout()
plt.show()
