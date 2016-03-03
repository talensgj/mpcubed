#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

# For A4 paper.
rcParams['figure.figsize'] = (11.69,8.27)
rcParams['xtick.labelsize'] = 'medium'
rcParams['ytick.labelsize'] = 'medium'
rcParams['axes.labelsize'] = 'large'
rcParams['image.interpolation'] = 'nearest'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package.models import transit

from transit_search import read_header, read_data

def plot_candidates(filename, data, outdir):
    
    ASCC, ra, dec, vmag, sptype = read_header(data)
    
    # Read the box least-squares results.
    with h5py.File(filename, 'r') as f:
        
        grp = f['header']
        ascc = grp['ascc'].value
        chisq0 = grp['chisq0'].value
        period = grp['period'].value
        best_depth = grp['depth'].value
        best_epoch = grp['epoch'].value
        best_duration = grp['duration'].value
        best_nt = grp['nt'].value
        flag = grp['flag'].value
        
        grp = f['data']
        freq = grp['freq'].value
        dchisq = grp['dchisq'].value
        depth = grp['depth'].value
        duration = grp['duration'].value
        epoch = grp['epoch'].value
        nt = grp['nt'].value
    
    # Determine if there are any transit candidates.
    select = (flag < 1)
    if np.sum(select) < 1: 
        return
    
    # Read the data.
    jdmid, lst, mag, emag, mask, trend = read_data(data, ascc)
    jdmid = jdmid - np.amin(jdmid)
    
    for i in range(len(ascc)):
        
        if (flag[i] > 0): continue 
        
        phase = np.mod((jdmid - best_epoch[i])/period[i] + .5, 1.) - .5
        
        # The best fit box-model.
        model = transit.box_model(jdmid, period[i], best_epoch[i], -best_depth[i], best_duration[i])
        weights = np.where(mask, 0., 1./emag**2)
        mean = np.sum(weights[i]*(mag[i] - model))/np.sum(weights[i])
        model = model + mean
        
        sort = np.argsort(phase)
        mphase = phase[sort]
        model = model[sort]
        
        # Bin the phase-folded lightcurve.
        nbins = np.ceil(9*period[i]/best_duration[i])
        bins = np.linspace(-.5, .5, nbins + 1)
        
        x0, _ = np.histogram(phase, bins, weights=weights[i])
        x1, _ = np.histogram(phase, bins, weights=weights[i]*mag[i])
        
        bphase = (bins[:-1] + bins[1:])/2.
        bmag = x1/x0
        ebmag = 1/np.sqrt(x0)
        
        # Make the figure.
        fig = plt.figure()
        
        gs = gridspec.GridSpec(3, 2, height_ratios = [.5,10,10])
    
        plt.suptitle('ASCC {}, {}, $V={:.1f}$\n$\delta={:.1f}$%, $P={:.2f}$ days, $\eta={:.2f}$ days, $n_t = {:d}$'.format(ascc[i], sptype[ASCC==ascc[i]][0], vmag[ASCC==ascc[i]][0], best_depth[i]*100, period[i], best_duration[i], int(best_nt[i])), size='xx-large')
        
        ax = plt.subplot(gs[1,:])
        plt.plot(freq, dchisq[:,i], c='k')
        #plt.axvline(1/period[i], c='g', ls='--')
        for n in range(1, 5):
            plt.axvline(n/period[i], c='r', ls='--')
            plt.axvline(1/(n*period[i]), c='r', ls='--')
        plt.xlim(0, 1.8)
        plt.xlabel('Frequency [day$^{-1}$]')
        plt.ylabel(r'$\Delta\chi^2$')
        
        plt.subplot(gs[2,:])
        plt.errorbar(phase[~mask[i]], mag[i,~mask[i]], emag[i,~mask[i]], fmt='.', c='k')
        #plt.errorbar(bphase, bmag, ebmag, fmt='.', c='g')
        plt.plot(mphase, model, c='r', lw=2)
        plt.xlim(-.5, .5)
        plt.ylim(.1, -.1)
        plt.xlabel('Phase')
        plt.ylabel(r'$\Delta m$')
        
        plt.tight_layout()
        #plt.savefig('/home/talens/a4.png')
        plt.show()
        plt.close()
        
    return


def main():
    
    filename = '/data2/talens/2015Q2_vmag/boxlstsq/bls0_vmag_2015Q2_patch266.hdf5'
    data = ['/data2/talens/2015Q2_vmag/LPN/red0_vmag_2015Q2LPN.hdf5',
            '/data2/talens/2015Q2_vmag/LPE/red0_vmag_2015Q2LPE.hdf5',
            '/data2/talens/2015Q2_vmag/LPS/red0_vmag_2015Q2LPS.hdf5',
            '/data2/talens/2015Q2_vmag/LPW/red0_vmag_2015Q2LPW.hdf5',
            '/data2/talens/2015Q2_vmag/LPC/red0_vmag_2015Q2LPC.hdf5']
    outdir = 'bla'
    
    plot_candidates(filename, data, outdir)
    
    return

if __name__ == '__main__':
    main()
