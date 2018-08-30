# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.spatial import ConvexHull
import scipy.sparse as sparse
from numpy.linalg import norm
from sklearn import preprocessing

import rampy

def get_portion_interest(x,y,bir):
    """Extracts the signals indicated in the bir.

    Parameters
    ----------

    x: the x axis

    y: the y values

    bir: the regions x values where the signal needs to be extracted, must be a n x 2 dimension array

    Returns
    -------
    yafit
        a 2 columns x-y array containing the signals in the bir.

    """
    birlen = np.array(bir.shape[0])

    sp = np.transpose(np.vstack((x.reshape(-1),y.reshape(-1))))
    ### selection of bir data
    for i in range(birlen):
        if i == 0:
            yafit = sp[np.where((sp[:,0]> bir[i,0]) & (sp[:,0] < bir[i,1]))]
        else:
            je = sp[np.where((sp[:,0]> bir[i,0]) & (sp[:,0] < bir[i,1]))]
            yafit = np.vstack((yafit,je))

    return yafit

def baseline(x_input,y_input,bir,method, **kwargs):
    """Allows subtracting a baseline under a x y spectrum.

    Parameters
    ----------
    x_input : ndarray
        x values.
    y_input : ndarray
        y values.
    bir : ndarray
        Contain the regions of interest, organised per line. 
        For instance, roi = np.array([[100., 200.],[500.,600.]]) will 
        define roi between 100 and 200 as well as between 500 and 600,.
    methods
        "poly": polynomial fitting, with splinesmooth the degree of the polynomial.
        "unispline": spline with the UnivariateSpline function of Scipy, splinesmooth is 
                     the spline smoothing factor (assume equal weight in the present case);
        "gcvspline": spline with the gcvspl.f algorythm, really robust. 
                     Spectra must have x, y, ese in it, and splinesmooth is the smoothing factor;
                     For gcvspline, if ese are not provided we assume ese = sqrt(y). 
                     Requires the installation of gcvspline with a "pip install gcvspline" call prior to use;
        "exp": exponential background;
        "log": logarythmic background;
        "rubberband": rubberband baseline fitting;
        "als": automatic least square fitting following Eilers and Boelens 2005;
        "arPLS": automatic baseline fit using the algorithm from Baek et al. 2015 
                 Baseline correction using asymmetrically reweighted penalized least squares smoothing, Analyst 140: 250-257.

    kwargs
    ------
    polynomial_order : Int
        The degree of the polynomial (0 for a constant), default = 1.
    s : Float
        spline smoothing coefficient for the unispline and gcvspline algorithms.
    lam : Float
        float, the lambda smoothness parameter for the ALS and ArPLS algorithms. Typical values are between 10**2 to 10**9, default = 10**5.
    p : Float
        float, for the ALS algorithm, advised value between 0.001 to 0.1, default = 0.01.
    niter : Int
        number of iteration of the ALS algorithm, default = 10.
    p0_exp : List
        containg the starting parameter for the exp baseline fit with curve_fit. Default = [1.,1.,1.].

    p0_log : List
        containg the starting parameter for the log baseline fit with curve_fit. Default = [1.,1.,1.,1.].

    Returns
    -------
    out1 : ndarray
        Contain the corrected signal.
    out2 : ndarray
        Contain the baseline.

    """
    # we get the signals in the bir
    yafit_unscaled = get_portion_interest(x_input,y_input,bir)

    # signal standard standardization with sklearn
    # this helps for polynomial fitting
    X_scaler = preprocessing.StandardScaler().fit(x_input.reshape(-1, 1))
    Y_scaler = preprocessing.StandardScaler().fit(y_input.reshape(-1, 1))

    # transformation
    x = X_scaler.transform(x_input.reshape(-1, 1))
    y = Y_scaler.transform(y_input.reshape(-1, 1))

    yafit = np.copy(yafit_unscaled)
    yafit[:,0] = X_scaler.transform(yafit_unscaled[:,0].reshape(-1, 1))[:,0]
    yafit[:,1] = Y_scaler.transform(yafit_unscaled[:,1].reshape(-1, 1))[:,0]

    y = y.reshape(len(y_input))

    if method == 'poly':

        # optional parameters
        poly_order = kwargs.get('polynomial_order',1)

        coeffs = np.polyfit(yafit[:,0],yafit[:,1],poly_order)

        baseline_fitted = np.polyval(coeffs,x)

    elif method == 'unispline':

        # optional parameters
        splinesmooth = kwargs.get('s',2.0)

        # fit of the baseline
        coeffs = UnivariateSpline(yafit[:,0],yafit[:,1], s=splinesmooth)

        baseline_fitted = coeffs(x)

    elif method == 'gcvspline':

        try:
            from gcvspline import gcvspline, splderivative
        except ImportError:
            print('ERROR: Install gcvspline to use this mode (needs a working FORTRAN compiler).')
            
        # optional parameters
        splinesmooth = kwargs.get('s',2.0)

        # Spline baseline with mode 1 of gcvspl.f, see gcvspline documentation
        c, wk, ier = gcvspline(yafit[:,0],yafit[:,1],np.sqrt(np.abs(yafit[:,1])),splinesmooth,splmode = 1) # gcvspl with mode 1 and smooth factor

        baseline_fitted = splderivative(x,yafit[:,0],c)

    elif method == 'gaussian':
        ### Baseline is of the type y = a*exp(-log(2)*((x-b)/c)**2)
        # optional parameters
        p0_gauss = kwargs.get('p0_gaussian',[1.,1.,1.])
        ## fit of the baseline
        coeffs, pcov = curve_fit(rampy.gaussian,yafit[:,0],yafit[:,1],p0 = p0_gauss)

        baseline_fitted = rampy.gaussian(x,coeffs[0],coeffs[1],coeffs[2])

    elif method == 'exp':
        ### Baseline is of the type y = a*exp(b*(x-xo))
        # optional parameters
        p0_exp = kwargs.get('p0_exp',[1.,1.,1.])
        ## fit of the baseline
        coeffs, pcov = curve_fit(rampy.funexp,yafit[:,0],yafit[:,1],p0 = p0_exp)

        baseline_fitted = rampy.funexp(x,coeffs[0],coeffs[1],coeffs[2])

    elif method == 'log':
        ### Baseline is of the type y = a*exp(b*(x-xo))
        # optional parameters
        p0_log = kwargs.get('p0_log',[1.,1.,1.,1.])
        ## fit of the baseline
        coeffs, pcov = curve_fit(rampy.funlog,yafit[:,0],yafit[:,1],p0 = p0_log)

        baseline_fitted = rampy.funlog(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3])

    elif method == 'rubberband':
        # code from this stack-exchange forum
        #https://dsp.stackexchange.com/questions/2725/how-to-perform-a-rubberband-correction-on-spectroscopic-data

        # Find the convex hull
        v = ConvexHull(np.array([x, y])).vertices

        # Rotate convex hull vertices until they start from the lowest one
        v = np.roll(v, -v.argmin())
        # Leave only the ascending part
        v = v[:v.argmax()]

        # Create baseline using linear interpolation between vertices
        baseline_fitted = np.interp(x, x[v], y[v])

    elif method == 'als':
        # Matlab code in Eilers et Boelens 2005
        # Python addaptation found on stackoverflow: https://stackoverflow.com/questions/29156532/python-baseline-correction-library

        # optional parameters
        lam = kwargs.get('lam',1.0*10**5)
        p = kwargs.get('p',0.01)
        niter = kwargs.get('niter',10)

        # starting the algorithm
        L = len(y)
        D = sparse.csc_matrix(np.diff(np.eye(L), 2))
        w = np.ones(L)
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = sparse.linalg.spsolve(Z, w*y)
            w = p * (y > z) + (1-p) * (y < z)

        baseline_fitted = z

    elif method == 'arPLS':
        # Adaptation of the Matlab code in Baek et al 2015

        # optional parameters
        lam = kwargs.get('lam',1.0*10**5)
        ratio = kwargs.get('ratio',0.01)

        N = len(y)
        D = sparse.csc_matrix(np.diff(np.eye(N), 2))
        w = np.ones(N)

        while True:
            W = sparse.spdiags(w, 0, N, N)
            Z = W + lam * D.dot(D.transpose())
            z = sparse.linalg.spsolve(Z, w*y)
            d = y - z
            # make d- and get w^t with m and s
            dn = d[d<0]
            m = np.mean(dn)
            s = np.std(dn)
            wt = 1.0/(1 + np.exp( 2* (d-(2*s-m))/s ) )
            # check exit condition and backup
            if norm(w-wt)/norm(w) < ratio:
                break
            w = wt

        baseline_fitted = z

    return y_input.reshape(-1,1)-Y_scaler.inverse_transform(baseline_fitted.reshape(-1, 1)), Y_scaler.inverse_transform(baseline_fitted.reshape(-1, 1))
    #return y_corrected, baseline_fitted

    #
