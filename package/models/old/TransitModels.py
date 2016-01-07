#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import numpy as np
from occultquad import occultquad

def BoxModel(t, t0, P, q, delta):
    
    reduced_phase = (t-t0)/P%1
    
    model = np.zeros(len(t))

    here = (reduced_phase<q/2.)|(reduced_phase>(1-q/2.))
    model[here] -= delta
    
    return model
    
def BoxLikeModel(t, t0, P, eta, delta, c):
    
    phase = (t-t0)/P
    
    tp = P*np.sin(np.pi*phase)/(np.pi*eta)
    
    mu = 0.5*delta*(-np.tanh(c*(tp+.5))+np.tanh(c*(tp-.5)))
    
    return mu

def _get_E(M, e):
    
    # Iteratively solve for the eccentric anomaly.
    E = M
    eps = 1e-7
    while (np.abs(E-e*np.sin(E)-M)>eps):
        E = E - (E-e*np.sin(E)-M)/(1.-e*np.cos(E))
        
    return E

# Quick vectorization of _get_E.
get_E = np.vectorize(_get_E)
    
def SkySeparation(t, t_p, P, aRs, e, i, w):
    
    # Mean anomaly.
    M = 2.*np.pi*(t-t_p)/P
    
    # Eccentric anomaly.
    if e<1e-7:
        print 'Taking e=0.'
        e = 0.
        E = M
    else:
        E = get_E(M, e)
    
    # True anomaly
    v = 2.*np.arctan2(np.sqrt(1.+e)*np.sin(E/2.), np.sqrt(1.-e)*np.cos(E/2.))
    
    # Separation.
    r = aRs*(1.-e*e)/(1.+e*np.cos(v))
    d = r*np.sqrt(1.-np.sin(w+v)*np.sin(w+v)*np.sin(i)*np.sin(i))
    
    return d
    
def MA_Model(t, t_p, P, p, aRs, e, i, w, u1, u2):
    
    z = SkySeparation(t, t_p, P, aRs, e, i, w)
    
    model, limb_model = occultquad(z, u1, u2, p)
    
    return limb_model


def main():
    
    import matplotlib.pyplot as plt
    
    time = np.linspace(0,10,100000)
    
    # Box model params.
    t0 = 2.
    P = 2.23
    q = .1
    delta = .005
    
    # Box-like model params.
    eta = q*P
    c = 50.
    
    # Mandel & Agol assuming cicrcular orbit. 
    tp = t0+P/4.
    p = np.sqrt(delta)
    aRs = 1./(np.pi*q)
    e = 0.
    i = np.pi/2.
    w = 0.
    u1 = 0.
    u2 = 0.
    
    
    box_model = BoxModel(time, t0, P, q, delta)
    boxlike_model = BoxLikeModel(time, t0, P, eta, delta, c)
    MA_model = MA_Model(time, tp, P, p, aRs, e, i, w, u1, u2)
    
    plt.plot(time, box_model)
    plt.plot(time, boxlike_model)
    plt.plot(time, MA_model-1)
    plt.show()
	
    return 0

if __name__ == '__main__':
	main()

