#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import pandas as pd

class StarCat(object):
    
    def __init__(self, catalogue='/data3/talens/Catalogue/I_280B.hdf5'):
        
        cat = dict()
        
        with h5py.File(catalogue, 'r') as f:
        
            grp = f['data']
            
            ascc = grp['ASCC'].value 
            cat['ra'] = grp['RAhour'].value
            cat['dec'] = grp['DEdeg'].value
            cat['plx'] = grp['Plx'].value
            cat['sptype'] = grp['SpType'].value
            cat['vmag'] = grp['Vmag'].value
            cat['bmag'] = grp['Bmag'].value
            cat['tyc1'] = grp['TYC1'].value
            cat['tyc2'] = grp['TYC2'].value
            cat['tyc3'] = grp['TYC3'].value
            cat['hd'] = grp['HD'].value
    
        cat['plx'] = cat['plx']/1e2 # [mas]
        cat['vmag'] = cat['vmag']/1e3 # [mag]
        cat['bmag'] = cat['bmag']/1e3 # [mag]
    
        ascc = ascc.astype('|S32')
        self.cat = pd.DataFrame(cat, index=ascc)
    
        return
        
    def get_star(self, ascc0):
        
        return self.cat.loc[ascc0]
