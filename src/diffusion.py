# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 15:18:12 2015

@author: Charles Le Losq
Research School of Earth Science

Functions for diffusion experiments
"""

from scipy.special import erfc
import numpy as np

def diffshort(x, t, C0, C1, D):
    """
    Simple equation for the diffusion into a semi-infinite slab, see Crank 1975
    C0 is the concentration in the core
    C1 is the concentration at the border
    D is the diffusion coefficient in log10 unit, m^2.s^-1
    x is the profil length in meters
    t is the time in seconds
    """
        
    Cx = (C1 - C0) * erfc(x / (2. * np.sqrt((10**D)*t))) + C0
    
    return Cx
    
def difffull(x1, x2, t, C0, C1, D):
    """
    Simple equation for the diffusion into a semi-infinite slab, see Crank 1975
    C0 is the concentration in the core
    C1 is the concentration at the border
    D is the diffusion coefficient in log10 unit, m^2.s^-1
    x1 and x2 are the profil lengths from beginning and end respectively, in meters
    t is the time in seconds
    """
        
    Cx = (C1 - C0) * ( erfc(x / (2. * np.sqrt((10**D)*t))) +  erfc((x2-x1) / (2. * np.sqrt((10**D)*t))))+ C0
    
    return Cx
