#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import h5py
import numpy as np

class fLCfile(object):
    """ Read data from fLC files.
    
    Attributes:
        fLCfile (str): The full path to the file.
    
    """
    
    def __init__(self, fLCfile):
        """ Initialize a reader of fLC files.
        
        Args:
            fLCfile (str): The full path to the file.
        
        """
        
        self.fLCfile = fLCfile
    
        return
        
    def read_global(self):
        """ Read the global group.
        
        Returns:
            data: A list of attribute (key, value) pairs.
        
        """
        
        with h5py.File(self.fLCfile, 'r') as f:
            data = f['global'].attrs.items()
        
        return data
        
    def read_header(self, fields, ascc=None):
        """ Read the specified fields from the header_table.
        
        Args:
            fields (str): List of fieldnames.
            ascc (str): List of ascc numbers for which the fields must be
                returned. Default is None.
        
        Returns:
            data: A list of data corresponding to the fields.
        
        """
        
        if ascc is None:
            
            # If no ascc are given read data for all stars.
            with h5py.File(self.fLCfile, 'r') as f:
                data = [f['header_table/'+field].value for field in fields]
            
        else:
            nstars = len(ascc)
            nfields = len(fields)
            
            # Create arrays.
            data = [np.zeros(nstars) for _ in range(nfields)]
            
            # Read the data for the specified stars.
            with h5py.File(self.fLCfile, 'r') as f:
                
                all_ascc = f['header_table/ascc'].value
                here = np.in1d(all_ascc, ascc)
                
                for j in range(nfields):
                    data[j] = f['header_table/'+fields[j]].value[here]
                    
        return data
    
    def read_data(self, fields, ascc=None, nobs=None):
        """ Read the specified fields from the lightcurves.
            
            Args:
                fields (str): List of fieldnames.
                ascc (str): List of ascc numbers for which the fields must be
                    returned. Default is None.
                nobs (int): Array containning the number of points in each
                    ligthcurve to be returned must correspond to ascc.
                    Default is None.
                    
            Returns:
                data: A lsit of data corresponding to the fields.
                
        """
        
        if ascc is None:
            # If no stars are specified read stars and nobs from the header.
            ascc, nobs = self.read_header(['ascc', 'nobs'])
        elif nobs is None:
            # If stars are specified but no nobs read nobs from the header.
            nobs, = self.read_header(['nobs'], ascc)

        nobs = nobs.astype('int')            
            
        nstars = len(ascc)
        nfields = len(fields)
        npoints = np.sum(nobs)
        select = np.append(0, np.cumsum(nobs))
        
        # Create arrays.
        data = [np.zeros(npoints) for _ in range(nfields)]
        
        # Read the data.
        with h5py.File(self.fLCfile, 'r') as f:
            for i in range(nstars):
                lc = f['data/'+ascc[i]].value
                for j in range(nfields):
                    data[j][select[i]:select[i+1]] = lc[fields[j]]
        
        return data
           
    def read_star(self, ascc):
        """ Read the lightcurve of the specified star.
        
        Args:
            ascc (str): The star to read.
        
        Returns:
            data: The full lightcurve of the star as a recarray.
        
        """
        
        with h5py.File(self.fLCfile, 'r') as f:
            data = f['data/' + ascc].value
        
        return data
        
    
