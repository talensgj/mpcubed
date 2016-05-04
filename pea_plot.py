#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

import matplotlib.pyplot as plt

from package import plotting
from package import IO
from package import misc

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
  
#def plot_PolarEqArea(grid, skymap, **kwargs):
    
    #from matplotlib import patches
    
    #N = grid.N
    #strides = np.cumsum(N)
    #strides = np.append(0, strides)
    
    #for i in range(grid.nrings + 2):
        
        #if (N[i] == 1):
            #pass
            
        #else:
            #xedges = np.linspace(0, 360, N[i]+1)
            #yedges = grid.yedges[i:i+2]
            
            #data = skymap[strides[i]:strides[i+1]]
            #data = data[:,None].T
            #data = np.ma.masked_invalid(data)
            
            #xedges, yedges = np.meshgrid(xedges, yedges)
            #yedges = yedges[::-1]
            
            #plt.pcolormesh(xedges*np.pi/180-np.pi, yedges*np.pi/180, data, **kwargs)
    
    #return
    
def _hadec2xy(ha, dec):
        
    from mascara import observer
        
    # Read the pointing from the file.
    #f = IO.SysFile(self.sysfile)
    #alt0, az0, th0, x0, y0 = f.read_pointing()
    
    # Initialize the coordinate transformations.
    site = observer.Site('LaPalma')
    cam = observer.Camera('central')
    #cam = observer.Camera(altitude=alt0, azimuth=az0, orientation=th0, Xo=x0, Yo=y0, nx=4008, ny=2672)
    
    tmp = ha.shape
    ha, dec = ha.ravel(), dec.ravel()
    
    # Perfrom the coordinate transformations.
    alt, az = site.hadec2altaz(ha, dec, degree=True)
    phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
    x, y = cam.PhiThe2XY(phi, theta)
    
    x, y = x.reshape(tmp), y.reshape(tmp)
    
    return x, y
    
def plot_Polar(grid, skymap, **kwargs):
    
    xedges = grid.xedges
    yedges = grid.yedges
    
    xedges, yedges = np.meshgrid(xedges, yedges)
    xedges, yedges = _hadec2xy(xedges, yedges)
    
    skymap = skymap[1:-1,1:-1].T
    skymap = np.ma.masked_invalid(skymap)
    
    plt.pcolormesh(xedges, yedges, skymap, **kwargs)
    
    return  
    
def plot_Healpix(grid, skymap, size=400, **kwargs):
    
    xedges = np.linspace(0, 360, size)
    yedges = np.linspace(-90, 90, size)
    
    x = (xedges[:-1] + xedges[1:])/2
    y = (yedges[:-1] + yedges[1:])/2
    
    x, y = np.meshgrid(x, y)
    
    idx = grid.radec2idx(x, y)
    
    data = skymap[idx]
    data = np.ma.masked_invalid(data)
    
    plt.pcolormesh(xedges*np.pi/180-np.pi, yedges*np.pi/180, data, **kwargs)
    
    return
    
def plot_PolarEqArea(grid, skymap, size=50, **kwargs):
    
    xedges = np.linspace(0, 360, size)
    yedges = np.linspace(-90, 90, size)
    
    x = (xedges[:-1] + xedges[1:])/2
    y = (yedges[:-1] + yedges[1:])/2
    
    x, y = np.meshgrid(x, y)
    x = x.ravel()
    y = y.ravel()
    
    _, _, idx = grid.radec2idx(x, y)
    
    data = skymap[idx]
    data = data.reshape((size-1, size-1))
    data = np.ma.masked_invalid(data)
    
    plt.pcolormesh(xedges*np.pi/180-np.pi, yedges*np.pi/180, data, **kwargs)
    
    return
    
def compare():

    file0 = '/data2/talens/2015Q2_vmag/LPS/sys0_vmag_201506ALPS.hdf5'
    #file0 = '/data2/talens/2015Q2_pea/LPE/sys0_pea_201506ALPE_hacells.hdf5'
    file1 = '/data2/talens/2015Q2_pea/LPS/sys0_pea_201506ALPS_hacells.hdf5'

    from pea_grid import PolarEAGrid

    f = IO.SysFile(file0)
    pgcam, trans0, nobs = f.read_trans()
    pg, sinx0, cosx0, siny0, cosy0, nobs = f.read_intrapix()
    
    hg, clouds0, sigma, nobs, lstmin, lstmax = f.read_clouds()
    
    #with h5py.File(file0, 'r') as f:
        #grp = f['data/clouds']
        #idx = grp['idx'].value
        #lstseq = grp['lstseq'].value
        #clouds0 = grp['clouds'].value
        #lstmin = grp.attrs['lstmin']
        #lstlen = grp.attrs['lstlen']
        
    #pea = PolarEAGrid(23)
        
    #tmp = np.full((pea.npix, lstlen), fill_value=np.nan)
    #tmp[idx, lstseq-lstmin] = clouds0
    #clouds0 = tmp

    f = IO.SysFile(file1)
    pg, trans1, nobs = f.read_trans()
    pg, sinx1, cosx1, siny1, cosy1, nobs = f.read_intrapix()
    
    with h5py.File(file1, 'r') as f:
        grp = f['data/clouds']
        idx = grp['idx'].value
        lstseq = grp['lstseq'].value
        clouds1 = grp['clouds'].value
        lstmin = grp.attrs['lstmin']
        lstlen = grp.attrs['lstlen']
        
    pea = PolarEAGrid(23)
        
    tmp = np.full((pea.npix, lstlen), fill_value=np.nan)
    tmp[idx, lstseq-lstmin] = clouds1
    clouds1 = tmp
    
    #ax = plt.subplot(211)
    #plt.imshow(clouds0, aspect='auto', interpolation='nearest', cmap=plotting.viridis, vmin=-.5, vmax=.5)
    #plt.colorbar()
    #plt.subplot(212, sharex=ax)
    #plt.imshow(clouds1, aspect='auto', interpolation='nearest', cmap=plotting.viridis, vmin=-.5, vmax=.5)
    #plt.colorbar()

    #plt.show()
    
    #frame = 0
    #for i in range(40000, 46000):
        
        #if np.all(np.isnan(clouds0[:,i])):
            #continue
    
        #fig = plt.figure(figsize=(10,20))
        
        #ax = plt.subplot(211, projection='mollweide')
        #plt.title('Sky 201506ALPE, idx = {}'.format(i))
        #plot_Healpix(hg, clouds0[:,i], vmin=-.5, vmax=.5, cmap=plotting.viridis)
        #ax.grid(True)
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        
        #ax = plt.subplot(212, projection='mollweide')
        #plot_PolarEqArea(pea, clouds1[:,i], vmin=-.5, vmax=.5, cmap=plotting.viridis)
        #ax.grid(True)
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        
        #plt.savefig('/data2/talens/sky/frame_{:05d}.png'.format(frame))
        #plt.close()
        
        #frame += 1
    
    #exit()
    
    with h5py.File('/data2/talens/2015Q2/LPS/fLC_201506ALPS.hdf5', 'r') as f:
        
        grp = f['header_table']
        ascc = grp['ascc'].value
        ra = grp['ra'].value
        dec = grp['dec'].value
        vmag = grp['vmag'].value
        nobs = grp['nobs'].value
        
        sort = np.argsort(dec)
        ascc = ascc[sort]
        ra = ra[sort]
        dec = dec[sort]
        vmag = vmag[sort]
        nobs = nobs[sort]
        
        for i in range(0, len(ascc), 1000):
            
            lc = f['data/'+ascc[i]].value
            
            ha = np.mod(lc['lst']*15. - ra[i], 360.)
            dec_ = np.repeat(dec[i], nobs[i])
            
            idx1, idx2 = pgcam.radec2idx(ha, dec_)
            
            mag, emag = misc.flux2mag(lc['flux0'], lc['eflux0'])
            mag = mag - vmag[i]
            
            t0 = trans0[idx1, idx2]
            t1 = trans1[idx1, idx2]
            
            idx1, idx2 = pg.radec2idx(ha, dec_)
            
            ipx0 = sinx0[idx1, idx2]*np.sin(2*np.pi*lc['x']) + cosx0[idx1, idx2]*np.sin(2*np.pi*lc['x']) + siny0[idx1, idx2]*np.sin(2*np.pi*lc['y']) + cosy0[idx1, idx2]*np.cos(2*np.pi*lc['y'])
            ipx1 = sinx1[idx1, idx2]*np.sin(2*np.pi*lc['x']) + cosx1[idx1, idx2]*np.sin(2*np.pi*lc['x']) + siny1[idx1, idx2]*np.sin(2*np.pi*lc['y']) + cosy1[idx1, idx2]*np.cos(2*np.pi*lc['y'])
            
            idx0 = hg.radec2idx(ra[i], dec[i])
            #_, _, idx1 = pea.radec2idx(np.array([ra[i]]), np.array([dec[i]]))
            _, _, idx1 = pea.radec2idx(ha, dec_)
            
            c0 = clouds0[idx0, lc['lstseq']-lstmin]
            c1 = clouds1[idx1, lc['lstseq']-lstmin]
            
            plt.imshow(siny0, aspect='auto', interpolation='nearest', cmap=plotting.viridis, vmin=-.1, vmax=.1)
            plt.axvline(idx2[0], c='k', lw=2)
            plt.colorbar()
            plt.show()
            
            ax1 = plt.subplot(521)
            plt.plot(mag, '.')
            plt.plot(t0+ipx0+c0, '.')
            
            plt.subplot(522, sharex=ax1, sharey=ax1)
            plt.plot(mag, '.')
            plt.plot(t1+ipx1+c1, '.')
            plt.ylim(np.nanmedian(mag) - 1., np.nanmedian(mag)+1)
            
            ax2 = plt.subplot(523, sharex=ax1)
            plt.plot(mag-ipx0-c0, '.')
            plt.plot(t0, '.')
            
            plt.subplot(524, sharex=ax1, sharey=ax2)
            plt.plot(mag-ipx1-c1, '.')
            plt.plot(t1, '.')
            plt.ylim(np.nanmedian(mag) - 1., np.nanmedian(mag)+1)
            
            ax3 = plt.subplot(525, sharex=ax1)
            plt.plot(mag-t0-c0, '.')
            plt.plot(ipx0, '.')
            
            plt.subplot(526, sharex=ax1, sharey=ax3)
            plt.plot(mag-t1-c1, '.')
            plt.plot(ipx1, '.')
            plt.ylim(-1, 1)
            
            ax4 = plt.subplot(527, sharex=ax1)
            plt.plot(mag-t0-ipx0, '.')
            plt.plot(c0, '.')
            
            plt.subplot(528, sharex=ax1, sharey=ax4)
            plt.plot(mag-t1-ipx1, '.')
            plt.plot(c1, '.')
            plt.ylim(-1, 1)
            
            ax5 = plt.subplot(529, sharex=ax1)
            plt.plot(mag - t0-ipx0-c0, '.')
            
            print np.nanstd(mag - t0 - ipx0 - c0)
            
            plt.subplot(5,2,10, sharex=ax1, sharey=ax5)
            plt.plot(mag - t1-ipx1-c1, '.')
            plt.ylim(-1, 1)
            
            print np.nanstd(mag - t1 - ipx1 - c1)
            
            plt.show()
            plt.close()
            
        
    exit()
    
    fig = plt.figure(figsize=(16, 10))
    
    vmin = np.nanpercentile(trans0, 1)
    vmax = np.nanpercentile(trans0, 99)
    
    ax = plt.subplot(221, aspect='equal')
    plot_Polar(pgcam, trans0, vmin=vmin, vmax=vmax, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(222, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pgcam, trans1, vmin=vmin, vmax=vmax, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(223, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pgcam, trans0 - trans1, vmin=-.05, vmax=.05, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.tight_layout()
    plt.show()
    
    fig = plt.figure(figsize=(16, 10))
    
    ax = plt.subplot(221, aspect='equal')
    plot_Polar(pg, sinx0, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(222, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pg, sinx1, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(223, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pg, sinx0 - sinx1, vmin=-.05, vmax=.05, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.tight_layout()
    plt.show()
    
    fig = plt.figure(figsize=(16, 10))
    
    ax = plt.subplot(221, aspect='equal')
    plot_Polar(pg, cosx0, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(222, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pg, cosx1, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(223, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pg, cosx0 - cosx1, vmin=-.05, vmax=.05, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.tight_layout()
    plt.show()
    
    fig = plt.figure(figsize=(16, 10))
    
    ax = plt.subplot(221, aspect='equal')
    plot_Polar(pg, siny0, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(222, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pg, siny1, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(223, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pg, siny0 - siny1, vmin=-.05, vmax=.05, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.tight_layout()
    plt.show()
    
    fig = plt.figure(figsize=(16, 10))
    
    ax = plt.subplot(221, aspect='equal')
    plot_Polar(pg, cosy0, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(222, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pg, cosy1, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(223, aspect='equal', sharex=ax, sharey=ax)
    plot_Polar(pg, cosy0 - cosy1, vmin=-.05, vmax=.05, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.tight_layout()
    plt.show()
    exit()
    
    return
    
def main(args):
    
    compare()
    
    #plot_trans('/data2/talens/sys0_201506ALPW_23e13f.hdf5')
    #plot_intrapix()
    #plot_clouds()
    #plot_clouds2()
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
