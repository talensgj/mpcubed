#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
import numpy as np
 
def softbox(t, P, Tp, delta, eta, c=10., m=0.):
    
    x = P*np.sin(np.pi*(t - Tp)/P)/(np.pi*eta)
    model = .5*delta*(np.tanh(c*(x - .5)) - np.tanh(c*(x + .5)))
    
    return model + m
    
def softbox_der(t, P, Tp, delta, eta, c=10., m=0.):
    
    x = P*np.sin(np.pi*(t - Tp)/P)/(np.pi*eta)
    mod = .5*delta*(np.tanh(c*(x - .5)) - np.tanh(c*(x + .5)))
    tmp = .5*delta*c*(np.cosh(c*(x - .5))**-2 - np.cosh(c*(x + .5))**-2)
    
    der = np.zeros((len(t), 6))
    der[:,0] = tmp*(x/P - np.cos(np.pi*(t - Tp)/P)*(t - Tp)/(eta*P))
    der[:,1] = -tmp*np.cos(np.pi*(t - Tp)/P)/eta
    der[:,2] = mod/delta
    der[:,3] = -tmp*x/eta
    der[:,4] = .5*delta*((x-.5)*np.cosh(c*(x - .5))**-2 - (x+.5)*np.cosh(c*(x + .5))**-2)
    der[:,5] = 1.
 
    return der
    
    
def softbox_chisq(pars, t, y, weights):
    
    mod = softbox(t, *pars)
    chisq = np.sum(weights*(y - mod)**2)
    
    return chisq
    
def softbox_chisq_der(pars, t, y, weights):
    
    mod = softbox(t, *pars)
    der = softbox_der(t, *pars)
    
    if (len(pars) == 4):
        der = der[:,:4]
    if (len(pars) == 5):
        der = der[:,:5]
    
    chisq_der = -2*np.sum(weights*(y - mod)*der.T, axis=1)
        
    return chisq_der
    
    
def main(args):
    
    import matplotlib.pyplot as plt
    from scipy import optimize
    
    t = np.linspace(0, 30, 5000)
    pars = [3., .5, .02, .1, 5.]
    y = softbox(t, *pars)
    y = y + .002*np.random.randn(5000)
    weights = np.ones(y.shape)/(.002**2)
    
    res = optimize.minimize(softbox_chisq, [2.9, 0.48, 0.05, 0.2, 5.], args = (t, y, weights), method='BFGS', jac=softbox_chisq_der, options={'maxiter':5000, 'disp':True})
    
    print res.x
    
    fit = softbox(t, *res.x)
    
    plt.plot(t, y)
    plt.plot(t, fit)
    plt.show()
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
