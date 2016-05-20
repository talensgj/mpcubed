#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'
rcParams['axes.labelsize'] = 'x-large'
rcParams['image.interpolation'] = 'none'
rcParams['image.origin'] = 'lower'
rcParams['axes.titlesize'] = 'xx-large'

from package import plotting
from package import IO
from package import misc

def read_clouds(sysfile):
    
    from pea_grid import PolarEAGrid
    
    with h5py.File(sysfile, 'r') as f:
    
        grp = f['data/clouds']
        idx = grp['idx'].value
        lstseq = grp['lstseq'].value
        clouds = grp['clouds'].value
        nobs = grp['nobs'].value
        
        try: sigma = grp['sigma'].value
        except: sigma = None
        
        nx = grp.attrs['nx']
        lstmin = grp.attrs['lstmin']
        lstmax = grp.attrs['lstmax']
        lstlen = grp.attrs['lstlen']

    lstseq = lstseq - lstmin

    hg = PolarEAGrid(nx)
    
    tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
    tmp[idx, lstseq] = clouds
    clouds = tmp
    
    tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
    tmp[idx, lstseq] = nobs
    nobs = tmp
    
    if sigma is not None:
        tmp = np.full((hg.npix, lstlen), fill_value = np.nan)
        tmp[idx, lstseq] = sigma
        sigma = tmp

    return hg, clouds, sigma, nobs, lstmin, lstmax

def plot_clouds():
    
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
    
def _hadec2xy(ha, dec):
        
    from mascara import observer
        
    # Read the pointing from the file.
    #f = IO.SysFile(self.sysfile)
    #alt0, az0, th0, x0, y0 = f.read_pointing()
    
    # Initialize the coordinate transformations.
    site = observer.Site('LaPalma')
    cam = observer.Camera('east')
    #cam = observer.Camera(altitude=alt0, azimuth=az0, orientation=th0, Xo=x0, Yo=y0, nx=4008, ny=2672)
    
    tmp = ha.shape
    ha, dec = ha.ravel(), dec.ravel()
    
    # Perfrom the coordinate transformations.
    alt, az = site.hadec2altaz(ha, dec, degree=True)
    phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
    x, y = cam.PhiThe2XY(phi, theta)
    
    x, y = x.reshape(tmp), y.reshape(tmp)
    
    return x, y
    
def _add_hadecgrid():
        
        from mascara import observer
        
        # Read the pointing from the file.
        #f = IO.SysFile(self.sysfile)
        #alt0, az0, th0, x0, y0 = f.read_pointing()
        
        # Initialize the coordinate transformations.
        site = observer.Site('LaPalma')
        cam = observer.Camera('east')
        #cam = observer.Camera(altitude=alt0, azimuth=az0, orientation=th0, Xo=x0, Yo=y0, nx=4008, ny=2672)
        
        # Add lines of constant declination.
        ha = np.linspace(0, 360, 360)
        dec = np.linspace(-80, 80, 17)
        ha, dec = np.meshgrid(ha, dec)
        
        tmp = ha.shape
        
        ha, dec = ha.ravel(), dec.ravel()
        
        alt, az = site.hadec2altaz(ha, dec, degree=True)
        phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
        x, y = cam.PhiThe2XY(phi, theta)
        
        x, y = x.reshape(tmp), y.reshape(tmp)
        
        here = (x > -50) & (x < 4008+50) & (y > -50) & (y < 2672+50)
        x[~here] = np.nan
        y[~here] = np.nan
        
        plt.plot(x.T, y.T, c='k')
        
        # Add lines of constant hour angle.
        ha = np.linspace(0, 345, 24)
        dec = np.linspace(-80, 80, 160)
        ha, dec = np.meshgrid(ha, dec)
        
        tmp = ha.shape
        
        ha, dec = ha.ravel(), dec.ravel()
        
        alt, az = site.hadec2altaz(ha, dec, degree=True)
        phi, theta, goodpoint = cam.Hor2PhiThe(alt, az)
        x, y = cam.PhiThe2XY(phi, theta)
        
        x, y = x.reshape(tmp), y.reshape(tmp)
        
        here = (x > -50) & (x < 4008+50) & (y > -50) & (y < 2672+50)
        x[~here] = np.nan
        y[~here] = np.nan
        
        plt.plot(x, y, c='k')
        
        return
    
def plot_Polar(grid, skymap, **kwargs):
    
    xedges = grid.xedges
    yedges = grid.yedges
    
    xedges, yedges = np.meshgrid(xedges, yedges)
    xedges, yedges = _hadec2xy(xedges, yedges)
    
    skymap = skymap[1:-1,1:-1].T
    skymap = np.ma.masked_invalid(skymap)
    
    im = plt.pcolormesh(xedges, yedges, skymap, **kwargs)
    _add_hadecgrid()
    
    return im  
    
def plot_Healpix(grid, skymap, size=400, **kwargs):
    
    xedges = np.linspace(0, 360, size)
    yedges = np.linspace(-90, 90, size)
    
    x = (xedges[:-1] + xedges[1:])/2
    y = (yedges[:-1] + yedges[1:])/2
    
    x, y = np.meshgrid(x, y)
    
    idx = grid.radec2idx(x, y)
    
    data = skymap[idx]
    data = np.ma.masked_invalid(data)
    
    im = plt.pcolormesh(xedges*np.pi/180-np.pi, yedges*np.pi/180, data, **kwargs)
    
    return im 
  
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
    
def plot_PolarEA(grid, skymap, size=50, **kwargs):
    
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
    
    #plt.pcolormesh(xedges*np.pi/180-np.pi, yedges*np.pi/180, data, **kwargs)
    
    xedges, yedges = np.meshgrid(xedges, yedges)
    xedges, yedges = _hadec2xy(xedges, yedges)
    im = plt.pcolormesh(xedges, yedges, data, **kwargs)
    _add_hadecgrid()
    
    return im
   
def new_sysplot(filename):
    
    # Read the transmission and intrapixel variations.
    f = IO.SysFile(filename)
    pgcam, trans, nobs = f.read_trans()
    pgipx, sinx, cosx, siny, cosy, nobs = f.read_intrapix()
    
    # Plot the transmission.
    fig = plt.figure(figsize=(14,9))
    
    plt.suptitle('Transmission.', size='xx-large')
    
    gs = gridspec.GridSpec(2, 2, width_ratios = [15,.5], height_ratios = [1,10])
    
    ax = plt.subplot(gs[1,0], aspect='equal')
    
    vmin = np.nanpercentile(trans, 1)
    vmax = np.nanpercentile(trans, 99)
    delta = vmax - vmin
    vmin = vmin - .01*delta
    vmax = vmax + .01*delta
    
    im = plot_Polar(pgcam, trans, vmin=vmin, vmax=vmax, cmap=plotting.viridis)
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    cax = plt.subplot(gs[1,1])
    cb = plt.colorbar(im, cax=cax)
    cb.ax.invert_yaxis()
    cb.set_label(r'$\Delta m$')
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    # Plot the intrapixel variations.
    fig = plt.figure(figsize=(16, 10))

    plt.suptitle('Intrapixel variations.', size='xx-large')

    gs = gridspec.GridSpec(3, 3, width_ratios = [15,15,.5], height_ratios = [1,10,10])
    
    plt.subplot(gs[1,0], aspect='equal')
    plt.title(r'$\sin(2\pi x)$')
    im = plot_Polar(pgipx, sinx, vmin=-.1, vmax=.1, cmap=plotting.viridis) 
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[1,1], aspect='equal')
    plt.title(r'$\cos(2\pi x)$')
    im = plot_Polar(pgipx, cosx, vmin=-.1, vmax=.1, cmap=plotting.viridis) 
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[2,0], aspect='equal')
    plt.title(r'$\sin(2\pi y)$')
    im = plot_Polar(pgipx, siny, vmin=-.1, vmax=.1, cmap=plotting.viridis) 
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(gs[2,1], aspect='equal')
    plt.title(r'$\cos(2\pi y)$')
    im = plot_Polar(pgipx, cosy, vmin=-.1, vmax=.1, cmap=plotting.viridis) 
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    cax = plt.subplot(gs[1:,2])
    cb = plt.colorbar(im, cax=cax)
    cb.set_label('Amplitude')
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    return
    
def diff_trans(file0, file1):
    
    f = IO.SysFile(file0)
    pgcam, trans0, nobs = f.read_trans()
    
    f = IO.SysFile(file1)
    pgcam, trans1, nobs = f.read_trans()
    
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
    plt.close()
    
    return
    
def diff_ipx(file0, file1):
    
    f = IO.SysFile(file0)
    pg, sinx0, cosx0, siny0, cosy0, nobs = f.read_intrapix()
    
    f = IO.SysFile(file1)
    pg, sinx1, cosx1, siny1, cosy1, nobs = f.read_intrapix()
    
    fig = plt.figure(figsize=(16, 10))
    
    plt.subplot(3,4,1, aspect='equal')
    plot_Polar(pg, sinx0, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,2, aspect='equal')
    plot_Polar(pg, cosx0, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,3, aspect='equal')
    plot_Polar(pg, siny0, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,4, aspect='equal')
    plot_Polar(pg, cosy0, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,5, aspect='equal')
    plot_Polar(pg, sinx1, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,6, aspect='equal')
    plot_Polar(pg, cosx1, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,7, aspect='equal')
    plot_Polar(pg, siny1, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,8, aspect='equal')
    plot_Polar(pg, cosy1, vmin=-.1, vmax=.1, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,9, aspect='equal')
    plot_Polar(pg, sinx0 - sinx1, vmin=-.05, vmax=.05, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,10, aspect='equal')
    plot_Polar(pg, cosx0 - cosx1, vmin=-.05, vmax=.05, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,11, aspect='equal')
    plot_Polar(pg, siny0 - siny1, vmin=-.05, vmax=.05, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.subplot(3,4,12, aspect='equal')
    plot_Polar(pg, cosy0 - cosy1, vmin=-.05, vmax=.05, cmap=plotting.viridis)
    plt.colorbar()
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    return
    
def main(args):
    
    #filelist = glob.glob('/data2/talens/2015Q2_pea/LP?/sys0_pea_ra*')
    
    #for filename in filelist:
        #systematics = plotting.SysPlot(filename)
        #systematics.plot_trans(display=False, savefig=True)
        #systematics.plot_intrapix(display=False, savefig=True)
    
    #new_sysplot('/data2/talens/2015Q2_vmag/LPE/sys0_vmag_201506BLPE.hdf5')
    diff_trans('/data2/talens/2015Q2_pea/LPE/sys0_pea_ra_201506ALPE.hdf5', '/data2/talens/2015Q2_pea/LPE/sys0_pea_ha_201506ALPE.hdf5')
    #diff_ipx('/data2/talens/2015Q2_vmag/LPE/sys0_vmag_201506ALPE.hdf5', '/data2/talens/2015Q2_pea/LPE/sys0_pea_ha_201506ALPE.hdf5')
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
