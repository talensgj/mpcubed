#!/usr/bin/env python
# -*- coding: utf-8 -*-

def func(x):
    return x-1

a = 0.
b = 2.

for i in range(50):
    
    f = func((a+b)/2)
    
    if f < 0:
        a = (a+b)/2
    else:
        b = (a+b)/2
    
    print b, a
    
    
    


