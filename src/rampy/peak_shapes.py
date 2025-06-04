# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import numpy as np
from scipy.optimize import curve_fit
import plotly.graph_objs as go
from plotly.subplots import make_subplots

def gaussian(x: np.ndarray, amp, freq, HWHM) -> np.ndarray:
    """Computes a Gaussian peak.

    Args:
        x (np.ndarray): x axis
        amp (float or np.ndarray): Peak amplitude
        freq (float or np.ndarray): Peak position
        HWHM (float or np.ndarray): Peak half-width at half-maximum.

    Returns:
        np.ndarray: calculated peak.

    Notes:
        Formula used: \( \text{amp} \cdot \exp(-\log(2) \cdot ((x - \text{freq}) / \text{HWHM})^2) \).

    Examples:
        You can create a single peak like:

        >>> x = np.linspace(0, 10, 1000)
        >>> y = rampy.gaussian(x, 10., 3., 1.0)

        You can also create an array with several peaks in columns using arrays for the peak parameters:

        >>> x = np.linspace(0, 10, 1000)
        >>> y = rampy.gaussian(x.reshape(-1, 1), np.array([[10, 10.]]), np.array([[3, 7.]]), np.array([[1., 1.]]))

        In this case, `y` will be an array of shape (len(x), 2) with one peak per columns. Your peaks can even share parameters, here the HWHM:

        >>> y = rampy.gaussian(x.reshape(-1, 1), np.array([[10, 10.]]), np.array([[3, 7.]]), 2.0)
        
        The composite signal is simply given by

        >>> y_sum = np.sum(y, axis=1) 

    """
    return amp * np.exp(-np.log(2) * ((x - freq) / HWHM)**2)

def lorentzian(x: np.ndarray, amp, freq, HWHM) -> np.ndarray:
    """Computes a Lorentzian peak.

    Args:
        x (np.ndarray): x axis
        amp (float or np.ndarray): Peak amplitude
        freq (float or np.ndarray): Peak position
        HWHM (float or np.ndarray): Peak half-width at half-maximum.

    Returns:
        np.ndarray: calculated peak.

    Notes:
        Formula used: \( \text{amp} / (1 + ((x - \text{freq}) / \text{HWHM})^2) \).

    Examples:
        You can create a single peak like:

        >>> x = np.linspace(0, 10, 1000)
        >>> y = rampy.lorentzian(x, 10., 3., 1.0)

        You can also create an array with several peaks in columns using arrays for the peak parameters:

        >>> x = np.linspace(0, 10, 1000)
        >>> y = rampy.lorentzian(x.reshape(-1, 1), np.array([[10, 10.]]), np.array([[3, 7.]]), np.array([[1., 1.]]))

        In this case, `y` will be an array of shape (len(x), 2) with one peak per columns. Your peaks can even share parameters, here the HWHM:

        >>> y = rampy.lorentzian(x.reshape(-1, 1), np.array([[10, 10.]]), np.array([[3, 7.]]), 2.0)
        
        The composite signal is simply given by

        >>> y_sum = np.sum(y, axis=1) 

    """
    return amp / (1.0 + ((x - freq) / HWHM)**2)

def pseudovoigt(x: np.ndarray, amp, freq, HWHM, L_ratio) -> np.ndarray:
    """Computes a pseudo-Voigt peak.

    Args:
        x (np.ndarray): x axis
        amp (float or np.ndarray): Peak amplitude
        freq (float or np.ndarray): Peak position
        HWHM (float or np.ndarray): Peak half-width at half-maximum.
        L_ratio (float): Ratio of the Lorentzian component, must be between 0 and 1 (inclusive).

    Returns:
        np.ndarray: The computed pseudo-Voigt signal.

    Raises:
        ValueError: If `L_ratio` is not between 0 and 1.

    Notes:
        Formula used: \( (1 - L_{\text{ratio}}) \cdot \text{gaussian}(x, \text{amp}, \text{freq}, \text{HWHM}) + 
                      L_{\text{ratio}} \cdot \text{lorentzian}(x, \text{amp}, \text{freq}, \text{HWHM}) \).

    Examples:
        You can create a single peak like:

        >>> x = np.linspace(0, 10, 1000)
        >>> y = rampy.pseudovoigt(x, 10., 3., 1.0, 0.5)

        You can also create an array with several peaks in columns using arrays for the peak parameters:

        >>> x = np.linspace(0, 10, 1000)
        >>> y = rampy.pseudovoigt(x.reshape(-1, 1), np.array([[10, 10]]), np.array([[3, 7]]), np.array([[1., 1.]]), np.array([[0.5, 0.5]]))

        In this case, `y` will be an array of shape (len(x), 2) with one peak per columns. Your peaks can even share parameters, here the L_ratio:

        >>> y = rampy.pseudovoigt(x.reshape(-1, 1), np.array([[10, 10]]), np.array([[3, 7]]), np.array([[1., 1.]]), 0.5)
        
        The composite signal is simply given by

        >>> y_sum = np.sum(y, axis=1) 

    """
    if not 0 <= L_ratio <= 1:
        raise ValueError("L_ratio must be between 0 and 1.")
    
    return L_ratio * lorentzian(x, amp, freq, HWHM) + (1 - L_ratio) * gaussian(x, amp, freq, HWHM)

def pearson7(x, a0, a1, a2, a3):
    """Computes a Pearson VII peak.

    Args:
        x (np.ndarray): Positions at which the signal should be sampled.
        x (np.ndarray): x axis
        a0 (float or np.ndarray): Peak amplitude
        a1 (float or np.ndarray): Peak position
        a2 (float or np.ndarray): Peak width
        a3 (float or np.ndarray): Peak shape parameter.

    Returns:
        np.ndarray: The computed Pearson VII signal.

    Notes:
        Formula used: \( a_0 / ((1 + ((x - a_1) / a_2)^2 \cdot (2^{(1/a_3)} - 1))^{a_3}) \).

    Examples:
        You can create a single peak like:

        >>> x = np.linspace(0, 10, 1000)
        >>> y = rampy.pearson7(x, 10., 3., 1.0, 0.5)

        You can also create an array with several peaks in columns using arrays for the peak parameters:

        >>> x = np.linspace(0, 10, 1000)
        >>> y = rampy.pearson7(x.reshape(-1, 1), np.array([[10, 10]]), np.array([[3, 7]]), np.array([[1., 1.]]), np.array([[0.5, 0.5]]))

        In this case, `y` will be an array of shape (len(x), 2) with one peak per columns. Your peaks can even share parameters, here the a3:

        >>> y = rampy.pearson7(x.reshape(-1, 1), np.array([[10, 10]]), np.array([[3, 7]]), np.array([[1., 1.]]), 0.5)
        
        The composite signal is simply given by

        >>> y_sum = np.sum(y, axis=1) 
    """
    return a0 / ( (1.0 + ((x-a1)/a2)**2.0 * (2.0**(1.0/a3) -1.0))**a3 )

"""
def create_gauss():
    def gauss(x,amp,freq,HWHM,bcg,slope):
        return amp*np.exp(-np.log(2)*((x-freq)/HWHM)**2)+slope*x+bcg
    return gauss

def create_lorenz():
    def lorenz(x,amp,freq,HWHM,bcg,slope):
        return amp/(1+((x-freq)/HWHM)**2)+slope*x+bcg
    return lorenz

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
