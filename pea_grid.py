# -*- coding: utf-8 -*-

import h5py
import numpy as np

import healpy

import matplotlib.pyplot as plt

import mascara

def polar_eqarea(ra, dec, n, N1=1):
    
    # Find the rings.
    yedges = np.linspace(-90, 90, n+1)    
    ring = np.searchsorted(yedges, dec, 'right')
    
    # Clean up the out-of-bounds indices.
    ring[ring == 0] = 1
    ring[ring == (n+1)] = n
    ring = ring - 1 
    
    # Compute the number of cells in each ring.
    rid = np.arange(n)
    N = (np.cos(np.pi/n*rid) - np.cos(np.pi/n*(rid + 1)))/(1 - np.cos(np.pi/n))
    N = np.rint(N)
    
    cell = np.zeros(len(ra))
    for i in range(n):

        select = (ring == i)
        
        # Find the cells.
        xedges = np.linspace(0, 360, N[i]+1)
        cell_ = np.searchsorted(xedges, ra[select], 'right')

        # Clean up the out-of-bounds indices.
        cell_[cell_ == 0] = 1
        cell_[cell_ == (N[i]+1)] = N[i]
        cell_ = cell_ - 1        
        
        cell[select] = cell_  
        
    return ring, cell, N      

def polar_eqarea_caps(ra, dec, n, N1=1):
    
    # Find the rings.
    yedges = np.linspace(-90.+90./(n+1), 90.-90./(n+1), n+1)
    yedges = np.hstack([-90., yedges, 90.])
    ring = np.searchsorted(yedges, dec, 'right')
    
    # Clean up the out-of-bounds indices.
    ring[ring == 0] = 1
    ring[ring == (n+3)] = n + 2
    ring = ring - 1 
    
    # Compute the number of cells in each ring.
    rid = np.arange(n)
    N = (np.cos(np.pi*(rid+.5)/(n+1)) - np.cos(np.pi*(rid + 1.5)/(n+1)))/(1 - np.cos(np.pi*.5/(n+1)))
    N = np.rint(N).astype('int')
    N = np.hstack([1, N, 1])
    
    cell = np.zeros(len(ra), dtype='int')
    for i in range(n+2):

        select = (ring == i)
        
        # Find the cells.
        xedges = np.linspace(0, 360, N[i]+1)
        cell_ = np.searchsorted(xedges, ra[select], 'right')

        # Clean up the out-of-bounds indices.
        cell_[cell_ == 0] = 1
        cell_[cell_ == (N[i]+1)] = N[i]
        cell_ = cell_ - 1        
        
        cell[select] = cell_  
        
    return ring, cell, N  

def test():

    ra = np.linspace(0, 360, 1001)
    dec = np.linspace(-90, 90, 1001)

    ra_gr = (ra[:-1] + ra[1:])/2
    dec_gr = (dec[:-1] + dec[1:])/2
    ra_gr, dec_gr = np.meshgrid(ra_gr, dec_gr)
    ra_gr = np.ravel(ra_gr)
    dec_gr = np.ravel(dec_gr)

    ring, cell, N = polar_eqarea_caps(ra_gr, dec_gr, 59)
    print np.amin(ring), np.amax(ring), len(N)
    idx = (np.cumsum(N) - N)[ring] + cell

    cell = cell.reshape((1000,1000))
    ring = ring.reshape((1000,1000))
    idx = idx.reshape((1000, 1000))

    #plt.subplot(111)
    #plt.pcolormesh(ra, dec, idx, cmap=plt.cm.Greys)
    #plt.colorbar()
    #plt.xlim(0, 360)
    #plt.ylim(-90, 90)
    #plt.show()

    plt.subplot(111, projection='hammer')
    plt.pcolormesh((ra-180)*np.pi/180, dec*np.pi/180, idx%2, cmap=plt.cm.Greys)
    plt.colorbar()
    plt.show()
    
    return
    
def histogram():
    
    filename = '/data2/talens/2015Q2/LPC/fLC_201506ALPC.hdf5'
    with h5py.File(filename, 'r') as f:
        
        grp = f['header_table']
        ra = grp['ra'].value
        dec = grp['dec'].value
        
    xedges = np.linspace(0, 360, 1001)
    yedges = np.linspace(-90, 90, 1001)

    x = (xedges[:-1] + xedges[1:])/2
    y = (yedges[:-1] + yedges[1:])/2
        
    x, y = np.meshgrid(x, y)
    x = np.ravel(x)
    y = np.ravel(y)
        
    ring, cell, N = polar_eqarea_caps(ra, dec, 23)
    idx = (np.cumsum(N) - N)[ring] + cell
    npbin1 = np.bincount(idx.astype('int'), minlength=np.sum(N).astype('int'))
    ring, cell, N = polar_eqarea_caps(x, y, 23)
    idx = (np.cumsum(N) - N)[ring] + cell
    idx = idx.astype('int')
    map1 = npbin1[idx]
    map1 = map1.reshape(1000, 1000)

    print np.sum(N)

    idx = healpy.ang2pix(8, (90. - dec)*np.pi/180., ra*np.pi/180.)
    npbin2 = np.bincount(idx, minlength=healpy.nside2npix(8))
    idx = healpy.ang2pix(8, (90. - y)*np.pi/180., x*np.pi/180.)
    map2 = npbin2[idx]
    map2 = map2.reshape(1000, 1000)

    plt.subplot(221)
    plt.pcolormesh(xedges, yedges, map1, cmap=plt.cm.Greys_r, vmin=0, vmax=250)
    plt.colorbar()
    #plt.plot(ra, dec, '.', c='r')
    plt.xlim(0, 360)
    plt.ylim(-90, 90)

    plt.subplot(222)
    plt.pcolormesh(xedges, yedges, map2, cmap=plt.cm.Greys_r, vmin=0, vmax=250)
    plt.colorbar()
    #plt.plot(ra, dec, '.', c='r')
    plt.xlim(0, 360)
    plt.ylim(-90, 90)

    ax3 = plt.subplot(223)
    plt.hist(npbin1[npbin1>0], bins=np.linspace(0, 250, 51))

    plt.subplot(224, sharey=ax3)
    plt.hist(npbin2[npbin2>0], bins=np.linspace(0, 250, 51))

    plt.show()

    return
  
def grid_on_CCD():
    
    xedges = np.linspace(0, 4008, 3*167+1)
    yedges = np.linspace(0, 2672, 2*167+1)

    x = (xedges[:-1] + xedges[1:])/2
    y = (yedges[:-1] + yedges[1:])/2

    x, y = np.meshgrid(x, y)
    x = np.ravel(x)
    y = np.ravel(y)

    # Create instances of site and camera.
    site = mascara.observer.Site('LaPalma')
    cam = mascara.observer.Camera('north')

    # Perfrom the coordinate transformations.
    phi, the = cam.XY2PhiThe(x, y)
    alt, az = cam.PhiThe2Hor(phi, the)
    ha, dec = site.altaz2hadec(alt, az)

    ring, cell, N = polar_eqarea_caps(ha, dec, 23)
    idx1 = (np.cumsum(N) - N)[ring] + cell
    idx1 = idx1.astype('int')
    idx1 = idx1.reshape(2*167, 3*167)

    idx2 = healpy.ang2pix(8, (90. - dec)*np.pi/180., ha*np.pi/180.)
    idx2 = idx2.reshape(2*167, 3*167)

    plt.subplot(211, aspect='equal')
    plt.pcolormesh(xedges, yedges, idx2, cmap=plt.cm.Greys)
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)

    plt.subplot(212, aspect='equal')
    plt.pcolormesh(xedges, yedges, idx1, cmap=plt.cm.Greys)
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)

    plt.show()

    return




    
    
    
    
