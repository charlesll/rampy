# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 18:30:44 2014

@author: Charles LE LOSQ
Carnegie Institution of Washington D.C.
February 2014
This function calculates the increase/decrease rate of each frequency in a spectrum
Your data must have the same frequency axis of course...
Initial fitting function is a second order polynomial
Change that following your needs
"""
import numpy as np
from scipy.optimize import curve_fit

def fun(x,a,b,c):
        y = a + b*x + c*x*x
        return y

def taux(x,spectres):
    
    # we need an organized function before calling the curve_fit algorithm
    freq = spectres[:,0]
    # output array
    taux = np.zeros((len(freq),4));
    taux[:,0] = freq
    
    # We look a each frequency, we sort y data and fit them with a second order polynomial
    for i in range(len(freq)):
        y = spectres[i,1::]
        popt, pcov = curve_fit(fun,x,y,[0.5e-3,0.5e-4,1e-6])
        taux[i,1:len(x)]=popt
        
    return taux

