#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import numpy as np

from numpy.lib.recfunctions import stack_arrays

from time import time

from collections import namedtuple

Quality = namedtuple('Quality', 'npoints npars chisq niter') 

class test():
    
    def __init__(self, maxiter=10):
        
        for self.niter in range(maxiter):
            self.func()
                
    def func(self):
        print self.niter

f = test()
