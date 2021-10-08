#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
from scipy.special import erfc

############ SIMPLE MATHEMATICAL FUNCTIONS ###########
def funlog(x,a,b,c,d):
    """log baseline function

    a * ln(-b *(x-c)) - d*x**2
    """
    return a*np.log(-b*(x-c))-d*x**2

def funexp(x,a,b,c):
    """exponential baseline function

    a*exp(b*(x-c))
    """
    return a*np.exp(b*(x-c))

def poly2(x,a,b,c):
    """returns a + b*x + c*x*x"""
    return a + b*x + c*x*x

def linear(x,a,b):
    """returns a + b*x"""
    return a + b*x

def linear0(x,a):
    """returns a*x"""
    return a*x

def constant(x,a):
    """returns a constant value

    Parameters
    ----------
    x : 1D array

    Returns
    -------
    y : 1D array
        array filled with a values
    """
    return np.zeros(len(x))+a

########### SPECIFIC GAUSSIAN/VOIGTR FUNCTIONS FOR USING WITH SCIPY OPTIMIZATION PROTOCOLS

def multigaussian(x,params):
    """old attempt to have a multigaussian function, do not use. Will be removed soon.

    """
    taille = len(params)
    y = np.zeros(len(x),4)
    for i in range(taille):
        y[:,taille+1] = params[taille,0]*np.exp(-np.log(2)*((x-params[taille,1])/params[taille,2])**2)
    y[:,0] = y[:,1]+y[:,2]+y[:,3]
    return y


def gauss_lsq(params,x):
    """predicts a sum of gaussian peaks with parameters params

    Parameters
    ----------
    params : 1D array
        an array of the parameters of the peaks. The number of peaks is assumed to be equal to len(params)/3.
        In this array, list intensities first, then all peak positions, then all peak half width at half maximum.
    x : 1D array
        x axis

    Returns
    -------
    y : 1D array
        y values at position x
    """
    nbpic = int(len(params)/3)
    a = np.zeros((1,nbpic))
    b = np.zeros((1,nbpic))
    c = np.zeros((1,nbpic))
    y = np.zeros((len(x),nbpic))
    for n in range(nbpic):
        m = 2*n # little trick for correct indexation
        a[0,n] = params[n+m]
        b[0,n] = params[n+m+1]
        c[0,n] = params[n+m+2]
        y[:,n] = a[0,n]*np.exp(-np.log(2)*((x[:]-b[0,n])/c[0,n])**2)
    ytot = sum(y,1)

    return ytot

def gauss_lsq_lfix(params,x):
    """predicts a sum of gaussian peaks with parameters params

    Assumes that all peaks share the same HWHM.

    Parameters
    ----------
    params : 1D array
        an array of the parameters of the peaks. The number of peaks is assumed to be equal to len(params)/3.
        In this array, list intensities first, then all peak positions, then the last element is the peaks' half width at half maximum.
    x : 1D array
        x axis

    Returns
    -------
    y : 1D array
        y values at position x
    """
    nbpic = int((len(params)-2)/2)
    a = np.zeros((1,nbpic))
    b = np.zeros((1,nbpic))
    c = np.zeros((1,2)) # FWMH fixed and in first position of params
    c[0,0] = params[0]
    c[0,1] = params[1]
    b[0,:] = params[2:(nbpic+2)]
    a[0,:] = params[nbpic+2:(2*nbpic+2)]
    y = np.zeros((len(x),nbpic))
    for n in range(nbpic):
        if n == 0:
            y[:,n] = a[0,n]*np.exp(-np.log(2)*((x[:]-b[0,n])/c[0,0])**2)
        else:
            y[:,n] = a[0,n]*np.exp(-np.log(2)*((x[:]-b[0,n])/c[0,1])**2)
    ytot = sum(y,1)

    return ytot

########### SPECIFIC FUNCTIONS FOR CHEMICAL DIFFUSION



def diffshort(x, t, C0, C1, D):
    """1D equation for the diffusion into a semi-infinite slab, see Crank 1975

    Parameters
    ----------
    C0 : float
        the concentration in the core
    C1 : float
        the concentration at the border
    D : float
        the diffusion coefficient in log10 unit, m^2.s^-1
    x : 1D array
        the profil length in meters
    t : float
        time in seconds

    Returns
    -------
    Cx : 1D array
        concentration at x
    """

    Cx = (C1 - C0) * erfc(x / (2. * np.sqrt((10**D)*t))) + C0

    return Cx

def difffull(x1, x2, t, C0, C1, D):
    """Equation for the diffusion into a full slab, see Crank 1975

    Here we assume the profil to have 2 surfaces of contact on each side

    Parameters
    ----------
    C0 : float
        the concentration in the core
    C1 : float
        the concentration at the border
    D : float
        the diffusion coefficient in log10 unit, m^2.s^-1
    x1 and x2 : float
        the profil lengths from beginning and end respectively, in meters
    t : float
        time in seconds
    """
    x = (x2-x1)
    Cx = (C1 - C0) * ( erfc(x / (2. * np.sqrt((10**D)*t))) +  erfc((x) / (2. * np.sqrt((10**D)*t))))+ C0

    return Cx
