#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package import plotting
from package import IO

def plot_trans(filename):
    
    systematics = plotting.SysPlot(filename)
    systematics.plot_trans()
    
    return
    
def plot_intrapix():
    
    systematics = plotting.SysPlot('/data2/talens/sys0_201506ALPW_2.hdf5')
    systematics.plot_intrapix()

    return

def plot_clouds():
    
    with h5py.File('/data2/talens/201506ALPS.hdf5', 'r') as f:
        
        grp = f['data/clouds']
        idx1 = grp['idx'].value
        lstseq = grp['lstseq'].value
        clouds = grp['clouds'].value
        
    lstseq, idx2 = np.unique(lstseq, return_inverse=True)
    N = np.amax(idx1) + 1
    
    print lstseq[6000], lstseq[7000]
    
    tmp = np.full((len(lstseq), N), fill_value=np.nan)
    tmp[idx2, idx1] = clouds

    plt.figure(figsize=(16,6))
    plt.imshow(tmp.T, aspect='auto', interpolation='nearest', vmin=-.5, vmax=.5)
    plt.colorbar()
    plt.show()
    
    return

def plot_clouds2():
    
    import matplotlib.gridspec as gridspec
    from matplotlib import cm
    from pea_grid import polar_eqarea_caps
    
    with h5py.File('/data2/talens/pea.hdf5', 'r') as f:
        
        grp = f['data/clouds']
        idx1 = grp['idx'].value
        lstseq = grp['lstseq'].value
        clouds = grp['clouds'].value
        
        lstmin = grp.attrs['lstmin']
        lstmax = grp.attrs['lstmax']
        lstlen = grp.attrs['lstlen']
    
    lstseq = lstseq - lstmin
    
    N = np.amax(idx1) + 1
    tmp = np.full((N, lstlen), fill_value=np.nan)
    tmp[idx1, lstseq] = clouds
    clouds = tmp
    
    xedges = np.linspace(0, 4008, 3*167)
    yedges = np.linspace(0, 2672, 2*167)
        
    xcenter = (xedges[:-1] + xedges[1:])/2.
    ycenter = (yedges[:-1] + yedges[1:])/2.
    
    xcenter, ycenter = np.meshgrid(xcenter, ycenter)
    systematics = plotting.SysPlot('/data2/talens/pea.hdf5')
    ha, dec = systematics._xy2hadec(xcenter, ycenter)
    
    for i in range(11955686-lstmin, 11956686-lstmin):
        
        if np.all(np.isnan(clouds[:,i])):
            continue
            
        print i + lstmin
        
        lst = ((i + lstmin) % 13500)*24./13500
        ra = np.mod(lst*15. - ha, 360.)
        
        tmp = ra.shape
        ra = ra.ravel()
        dec = dec.ravel()
        
        ring, cell, N = polar_eqarea_caps(ra, dec, 23)
        print np.sum(N)
        idx = (np.cumsum(N) - N)[ring] + cell 
        idx = idx.reshape(tmp)
        
        array = clouds[idx, i]
        array = np.ma.masked_invalid(array)
        
        fig = plt.figure(figsize=(14,9))

        #plt.suptitle('Clouds 2015-06-08 East', size='xx-large')

        gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
        
        plt.subplot(gs[1,0], aspect='equal')
        plt.annotate('seq = {}'.format(i + lstmin), (0,1), xycoords='axes fraction', ha = 'left', va='top', xytext=(5, -5), textcoords='offset points', size='large', backgroundcolor='w')
        
        im = plt.pcolormesh(xedges, yedges, array, cmap=cm.Greys, vmin=-.25, vmax=.25)
        systematics._add_hadecgrid()
        plt.xlim(0, 4008)
        plt.ylim(0, 2672)
        
        cax = plt.subplot(gs[1,1])
        cb = plt.colorbar(im, cax = cax)
        cb.ax.invert_yaxis()
        cb.set_label('Magnitude')
        
        plt.tight_layout()
        plt.savefig('/data2/talens/clouds/img_{}.png'.format(i+lstmin), dpi=160)
        #plt.show()
        plt.close()

    return
    
def compare():

    f = IO.SysFile('/data2/talens/polar_eqarea/201506ALPW.hdf5')
    pg, trans0, nosb = f.read_trans()

    f = IO.SysFile('/data2/talens/2015Q2_vmag/LPW/sys0_vmag_201506ALPW.hdf5')
    pg, trans1, nosb = f.read_trans()


    plt.imshow(trans0 - trans1, aspect='auto', interpolation='nearest', cmap=plotting.viridis, vmin=-.05, vmax=.05)
    plt.colorbar()
    plt.show()
    
    return
    
def main(args):
    
    #compare()
    
    plot_trans('/data2/talens/sys0_201506ALPW_23e13f.hdf5')
    #plot_intrapix()
    #plot_clouds()
    #plot_clouds2()
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
