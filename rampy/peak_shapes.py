#!/usr/bin/env python
#-*- coding: utf-8 -*-
import numpy as np

def gaussian(x,amp,freq,HWHM): # for spectral fit
    """compute a Gaussian peak

    Inputs
    ------
    x : ndarray
        the positions at which the signal should be sampled
    amp : float or ndarray with size equal to x.shape
        amplitude
    freq : float or ndarray with size equal to x.shape
        frequency/position of the Gaussian component
    HWHM : float or ndarray with size equal to x.shape
        half-width at half-maximum

    Returns
    -------
    out : ndarray
        the signal

    Remarks
    -------
    Formula is amp*np.exp(-np.log(2)*((x-freq)/HWHM)**2)
    """
    return amp*np.exp(-np.log(2)*((x-freq)/HWHM)**2)
	
def lorentzian(x,amp,freq,HWHM):
    """compute a Lorentzian peak

    Inputs
    ------
    x : ndarray
        the positions at which the signal should be sampled
    amp : float or ndarray with size equal to x.shape
        amplitude
    freq : float or ndarray with size equal to x.shape
        frequency/position of the Gaussian component
    HWHM : float or ndarray with size equal to x.shape
        half-width at half-maximum

    Returns
    -------
    out : ndarray
        the signal

    Remarks
    -------
    Formula is amp/(1+((x-freq)/HWHM)**2)
    """
    return amp/(1+((x-freq)/HWHM)**2)

def pseudovoigt(x,amp,freq,HWHM,L_ratio):
    """compute a pseudo-Voigt peak
    Inputs
    ------
    x : ndarray
        the positions at which the signal should be sampled. Can be provided as vector, nx1 or nxm array.
    amp : float or ndarray with size equal to x.shape
        amplitude
    freq : float or ndarray with size equal to x.shape
        frequency/position of the Gaussian component
    HWHM : float or ndarray with size equal to x.shape
        half-width at half-maximum
    L_ratio : float or ndarray with size equal to x.shape
        ratio pf the Lorentzian component, should be between 0 and 1 (included)

    Returns
    -------
    out : ndarray of size equal to x.shape
        the signal

    Remarks
    -------
    Formula is (1-L_ratio)*gaussian(amp,freq,HWHM) + L_ratio*lorentzian(amp,freq,HWHM)
    """
    try:
        if (L_ratio.any()>1) or (L_ratio.any()<0): # if entries are lists/arrays
            raise ValueError("L_ratio should be comprised between 0 and 1")
    except:
        if (L_ratio>1) or (L_ratio<0): # if entries are floats
            raise ValueError("L_ratio should be comprised between 0 and 1")

    return L_ratio*lorentzian(x,amp,freq,HWHM) + (1-L_ratio)*gaussian(x,amp,freq,HWHM)

def pearson7(x,a0,a1,a2,a3):
    """compute a Peason7 peak

    Inputs
    ------
    x : ndarray
        the positions at which the signal should be sampled
    a0, a1, a2, a3 : float or ndarrays of size equal to x.shape
        parameters of the Pearson7 equation

    Returns
    -------
    out : ndarray
        the signal

    Remarks
    -------
    Formula is a0 / ( (1.0 + ((x-a1)/a2)**2.0 * (2.0**(1.0/a3) -1.0))**a3 )
    """
    return a0 / ( (1.0 + ((x-a1)/a2)**2.0 * (2.0**(1.0/a3) -1.0))**a3 )
	
def create_gauss():
    def gauss(x,amp,freq,HWHM,bcg,slope):
        return amp*np.exp(-np.log(2)*((x-freq)/HWHM)**2)+slope*x+bcg
    return gauss

def create_lorenz():
    def lorenz(x,amp,freq,HWHM,bcg,slope):
        return amp/(1+((x-freq)/HWHM)**2)+slope*x+bcg
    return lorenz
"""
def create_pseudovoigt():
	def pseudovoigt(x,amp,freq,HWHM,L_ratio):
            try:
				if (L_ratio.any()>1) or (L_ratio.any()<0):
					raise ValueError("L_ratio should be comprised between 0 and 1")
			except:
				if (L_ratio.any()>1) or (L_ratio.any()<0):
					raise ValueError("L_ratio should be comprised between 0 and 1")

		return L_ratio*lorentzian(x,amp,freq,HWHM) + (1-L_ratio)*gaussian(x,amp,freq,HWHM)
    return pseudovoigt
"""
