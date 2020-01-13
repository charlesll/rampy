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
    x : ndarray
        the x axis
    y : ndarray
        the y values
    bir : n x 2 array
        the x values of regions where the signal needs to be extracted,
        must be a n x 2 dimension array, where n is the number of regions to extract
        and column 0 contains the low bounds, column 1 the high ones.

    Returns
    -------
    yafit : ndarray
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
        define roi between 100 and 200 as well as between 500 and 600.
        Note: This is NOT used by the "als" and "arPLS" algorithms, but still is a requirement when calling the function.
        bir and method probably will become args in a futur iteration of rampy to solve this.
    method : str
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
        "als": (automatic) baseline least square fitting following Eilers and Boelens 2005;
        "arPLS": (automatic) Baseline correction using asymmetrically reweighted penalized least squares smoothing. Baek et al. 2015, Analyst 140: 250-257;
        'drPLS': (automatic) Baseline correction method based on doubly reweighted penalized least squares. Xu et al., Applied Optics 58(14):3913-3920.

    kwargs
    ------
    polynomial_order : Int
        The degree of the polynomial (0 for a constant), default = 1.
    s : Float
        spline smoothing coefficient for the unispline and gcvspline algorithms.
    lam : Float
        float, the lambda smoothness parameter for the ALS, ArPLS and drPLS algorithms. Typical values are between 10**2 to 10**9, default = 10**5 for ALS and ArPLS and default = 10**6 for drPLS.
    p : Float
        float, for the ALS algorithm, advised value between 0.001 to 0.1, default = 0.01.
    ratio : float
        ratio parameter of the arPLS and drPLS algorithm. default = 0.01 for arPLS and 0.001 for drPLS.
    niter : Int
        number of iteration of the ALS and drPLS algorithm, default = 10 for ALS and default = 100 for drPLS.
    eta : Float
        roughness parameter for the drPLS algorithm, is between 0 and 1, default = 0.5
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

    elif method == 'drPLS':
        #according to Applied Optics, 2019, 58, 3913-3920.

        #optional parameters
        niter = kwargs.get('niter',100)
        lam = kwargs.get('lam',1000000)
        eta = kwargs.get('eta',0.5)
        ratio = kwargs.get('ratio',0.001)

        #optional smoothing in the next line, currently commented out
        #y = np.around(savgol_filter(raw_data,19,2,deriv=0,axis=1),decimals=6)

        L = len(y)

        D = sparse.diags([1,-2,1],[0,-1,-2],shape=(L,L-2),format='csr')
        D = D.dot(D.transpose())
        D_1 = sparse.diags([-1,1],[0,-1],shape=(L,L-1),format='csr')
        D_1 = D_1.dot(D_1.transpose())

        w_0 = np.ones(L)
        I_n = sparse.diags(w_0,format='csr')

        #this is the code for the fitting procedure
        w = w_0
        W = sparse.diags(w,format='csr')
        Z = w_0

        for jj in range(int(niter)):
            W.setdiag(w)
            Z_prev = Z
            Z = sparse.linalg.spsolve(W + D_1 + lam * (I_n - eta*W) * D,W*y,permc_spec='NATURAL')
            if np.linalg.norm(Z - Z_prev) > ratio:
                d = y - Z
                d_negative = d[d<0]
                sigma_negative = np.std(d_negative)
                mean_negative = np.mean(d_negative)
                w = 0.5 * (1 - np.exp(jj) * (d - (-mean_negative + 2*sigma_negative))/sigma_negative / (1 + np.abs(np.exp(jj) * (d - (-mean_negative + 2*sigma_negative))/sigma_negative)))
            else:
                break
        #end of fitting procedure

        baseline_fitted = Z
    else:
        raise ValueError("method not found, check you entered the right name.")

    return y_input.reshape(-1,1)-Y_scaler.inverse_transform(baseline_fitted.reshape(-1, 1)), Y_scaler.inverse_transform(baseline_fitted.reshape(-1, 1))
