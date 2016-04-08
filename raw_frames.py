#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob

import h5py
import numpy as np

import matplotlib.pyplot as plt

import healpy
from astropy.io import fits

from mascara.observer import Site, Camera
from mascara.funcs.utils import loadStellarCatalog

def main():
    
    filelist = glob.glob('/data2/talens/20160130LPC/*LPC.fits')
    filelist = np.sort(filelist)
    filelist = filelist[::50]
    
    ra, dec, mag, bmag, sID, sptype = loadStellarCatalog(Vmag=[4, 7.5])
    
    mysite = Site('LaPalma')
    mycamera = Camera('central')
    
    first=True
    for filename in filelist:
    
        hdulist = fits.open(filename)
        jd = hdulist[0].header['JD']
        jd = np.float64(jd)
        image = hdulist[0].data
    
        if first:
        
            vmin = np.percentile(image, 1)
            vmax = np.percentile(image, 99)
            
            alt, az, vis = mysite.apparent_stars(ra, dec, jd)
            vmag = mag[vis]
        
            #mycamera.find_pointing(image, alt, az, vmag)
            mycamera.solve_astrometry(image, alt, az, vmag)
        
            xedges = np.linspace(0, 4008, 2004+1)
            yedges = np.linspace(0, 2672, 1336+1)
            
            x = (xedges[:-1] + xedges[1:])/2
            y = (yedges[:-1] + yedges[1:])/2
            
            x, y = np.meshgrid(x, y)
            x = x.ravel()
            y = y.ravel()
            
            alt, az = mycamera.projection_on_sky(x, y)
            
            first = False
    
        image = image[0::2] + image[1::2]
        image = image[:,0::2] + image[:,1::2]
        image = image/4.
        image = image.ravel()
        
        ra, dec = mysite.Horizontal2Equatorial(alt, az, jd)
        ra = np.mod(ra + 180, 360)
        
        idx = healpy.ang2pix(1024, (90. - dec)*np.pi/180., ra*np.pi/180.)
        img = np.bincount(idx, image, minlength=healpy.nside2npix(1024))/np.bincount(idx, minlength=healpy.nside2npix(1024))
        
        # Plot the data.
        healpy.mollview(img, cmap=plt.cm.Greys, xsize=4000, min=vmin, max=vmax)
        plt.show()
        plt.close()
    
    return

if __name__ == '__main__':
   main()
