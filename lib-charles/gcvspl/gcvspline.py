# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 16:50:59 2014

@author: charleslelosq

Python Wrapper for the gcvspl.f spline function
It assumes that you already have the good path in your pathway
gcvspl.so should exist in this path...
"""

import gcvspl
import numpy as np

def gcvspline(x,y,ese,VAL,**options):
    """
    
    c, wk, ier = gcvspline(x,y,f,**options)
    
    Function for using the wrapper of gcvspl.f
    Input:
    x (1D array) are the independent variables
    y (2D array as (len(x),1)) are the observations
    (we assume here that you want to use this 
    spline only on 1 dataset... see gcvspl.f if not)
    ese are the errors on y
    VAL (double) depends on the mode, see gcvspl.f for details (=ese**2 for mode3)
    options:
    NC is the number of output spline coefficients, NC >= len(y) 
        default: NC = len(y)
    splorder is the half order of the required B-splines. The values 
    splorder = 1,2,3,4 correspond to linear, cubic, quintic, and 
    heptic splines, respectively.
        default: splorder = 2 (cubic) 
    splmode is the Optimization mode switch:
        default: splmode = 2 (General Cross Validated)
                       splmode = 1: Prior given value for p in VAL
                                 (VAL.ge.ZERO). This is the fastest
                                 use of GCVSPL, since no iteration
                                 is performed in p.
                       splmode = 2: Generalized cross validation.
                       splmode = 3: True predicted mean-squared error,
                                 with prior given variance in VAL.
                       splmode = 4: Prior given number of degrees of
                                 freedom in VAL (ZERO.le.VAL.le.N-M).
                       splmode  < 0: It is assumed that the contents of
                                 X, W, M, N, and WK have not been
                                 modified since the previous invoca-
                                 tion of GCVSPL. If MD < -1, WK(4)
                                 is used as an initial estimate for
                                 the smoothing parameter p.  At the
                                 first call to GCVSPL, MD must be > 0.
                       Other values for |MD|, and inappropriate values
                       for VAL will result in an error condition, or
                       cause a default value for VAL to be selected.
                       After return from MD.ne.1, the same number of
                       degrees of freedom can be obtained, for identical
                       weight factors and knot positions, by selecting
                       |MD|=1, and by copying the value of p from WK(4)
                       into VAL. In this way, no iterative optimization
                       is required when processing other data in Y.     
   
    see gcvspl.f for more info about the two last parameter
    
    Output:
    c: the spline
    wk: work vector, see gcvspl.f
    ier: error parameter. 
        ier = 0: Normal exit 
        ier = 1: M.le.0 .or. N.lt.2*M
        ier = 2: Knot sequence is not strictly
                 increasing, or some weight
                 factor is not positive.
        ier = 3: Wrong mode parameter or value.
    
    Wrapped by C. Le Losq. with using f2py. See the howto associated for wrapping.
    
    """
    # We will do the effort in the following to keep the same name as in gcvspl.f    
    # The following stuffs are useless since gcvspl will handle that
    # by itself. However the wrapper assumes the following statement:
    # dimensions = shape(y)
    # K = dimensions[1] # number of dataset, always 1 for our case
    # N = dimensions[0] # number of independent variable
    # NY = N # Number of observation, N = NY

    wx = 1/(ese**2) # relative variance of observations
    wy = np.zeros((1))+1 # systematic errors... not used so put them to 1
    
    # Default values
    if options.get("splorder") == None:
        splorder = 2
    else:
        splorder = options.get("splorder")   
        
    if options.get("splmode") == None:
        splmode = 2
    else:
        splmode = options.get("splmode")  
        
    if options.get("NC") == None:
        NC = len(y)
    else:
        NC = options.get("NC")      

    c, wk, ier = gcvspl.gcvspl(x,y,wx,wy,splorder,splmode,VAL,NC)
    
    return c, wk, ier

def splderivative(xfull,xparse,cparse,**options):
    """
    Wrapper to the SPLDER function of gcvspl.f
    Powerful for interpolation purpose
    
    xfull (1D array) contains the entire x range where the spline has 
    to be evaluated
    xparse (1D array) contains the x values of interpolation regions   
    WARNING xparse[0] <= xfull[0] <= xparse[n] 
    cparse (1D array) is the evaluated spline coefficients returned by gcvspl
    for xparse
    
    options:    
        splineorder (integer): is the spline order, 
        default: splineorder = 2 (cubic)
        L (integer): see gcvspl.f for details, default: L = 1
        IDER: the Derivative order required, with 0.le.IDER 
        and IDER.le.2*M. If IDER.eq.0, the function
        value is returned; otherwise, the IDER-th
        derivative of the spline is returned.
        
    see the splder function in gcvspl.f for more information
    """    
   
    # Default values
    if options.get("splineorder") == None:
        splineorder = 2
    else:
        splineorder = options.get("splineorder")
        
    if options.get("L") == None:
        L = 1
    else:
        L = options.get("L")
        
    if options.get("IDER") == None:
        IDER = 0
    else:
        IDER = options.get("IDER")
        
    q = np.zeros((2*splineorder))  # working array
    ycalc = np.zeros((len(xfull))) # Output array
    
    # we loop other xfull to create the output values
    for i in range(len(xfull)):
        ycalc[i] = gcvspl.splder(IDER,splineorder,xfull[i],xparse,cparse,L)
    
    return ycalc
