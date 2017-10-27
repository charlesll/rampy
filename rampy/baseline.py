import numpy as np
from gcvspline import gcvspline, splderivative 
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.spatial import ConvexHull
import scipy.sparse as sparse
from numpy.linalg import norm
from sklearn import preprocessing

def get_portion_interest(x,y,bir):
    """
    This function extracts the signals indicated in the bir.
    
    Inputs
    ------
    
    x: the x axis
    
    y: the y values
    
    bir: the regions x values where the signal needs to be extracted, must be a n x 2 dimension array
    
    Outputs
    -------
    
    an 2 columns x-y array containing the signals in the bir.
    
    """
    birlen = np.array(bir.shape[0])
    
    sp = np.transpose(np.vstack((x,y)))
    ### selection of bir data
    for i in range(birlen):
        if i == 0:
            yafit = sp[np.where((sp[:,0]> bir[i,0]) & (sp[:,0] < bir[i,1]))] 
        else:
            je = sp[np.where((sp[:,0]> bir[i,0]) & (sp[:,0] < bir[i,1]))] 
            yafit = np.vstack((yafit,je))
            
    return yafit
    

def baseline(x_input,y_input,bir,method, **kwargs):
    """
    This function allows subtracting a baseline under the spectra
    spectre is a spectrum or an array of spectra constructed with the spectrarray function
    bir contains the Background Interpolation Regions, it must be a n x 2 dimension array
    
    Inputs
    ------
    
        X: Array with x values.
    
        Y: array with y values.
    
        bir: an Array containing the regions of interest, organised per line. for instance, roi = np.array([[100., 200.],[500.,600.]]) will define roi between 100 and 200 as well as between 500 and 600,.
    
        methods:
    
            "poly": polynomial fitting, with splinesmooth the degree of the polynomial.
    
            "unispline": spline with the UnivariateSpline function of Scipy, splinesmooth is the spline smoothing factor (assume equal weight in the present case);
    
            "gcvspline": spline with the gcvspl.f algorythm, really robust. Spectra must have x, y, ese in it, and splinesmooth is the smoothing factor;
            for gcvspline, if ese are not provided we assume ese = sqrt(y);
    
            'exp': baseline of type 
    
            'log': baseline of type a*np.log(-b*(x-c))-d*x**2, parameters adjusted by scipy.optimize.curve_fit
    
            "rubberband": rubberband baseline fitting
    
            "als": automatic least square fitting following Eilers and Boelens 2005.
    
            "arPLS": automatic baseline fit using the algorithm from Baek et al. 2015 Baseline correction using asymmetrically reweighted penalized least suqares smoothing, Analyst 140: 250-257.
    
    Options
    -------
    
        Options for the ALS algorithm are:
    
        polynomial_order: integer, the degree of the polynomial (0 for a constant), default = 1
        
        s: float, spline smoothing coefficient for the unispline and gcvspline algorithms
    
        lam: float, between 10**2 to 10**9, default = 10**5
    
        p: float, advised value between 0.001 to 0.1, default = 0.01
        niter: number of iteration, default = 10
    
        Options for the ALS algorithm are:
    
        lam: float, the lambda smoothness parameter. Typical values are between 10**2 to 10**9, default = 10**5
    
        ratio = float, the termination ratio. Set it high to start and lower it. Default = 0.01.
    
    Outputs
    -------
    
        out1: an 2 columns x-y array containing the corrected signal
    
        out2: an 2 columns x-y array containing the baseline
    
        coefs: contains spline coefficients.
    
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

    y = y.reshape(len(y),)
    
    if method == 'poly':
        
        # optional parameters
        poly_order = kwargs.get('polynomial_order',1)
        
        coeffs = np.polyfit(yafit[:,0],yafit[:,1],poly_order)
                
        baseline_fitted = np.polyval(coeffs,x)
        y_corrected = y-baseline_fitted
        
    elif method == 'unispline':
        
        # optional parameters
        splinesmooth = kwargs.get('s',2.0)
        
        # fit of the baseline
        coeffs = UnivariateSpline(yafit[:,0],yafit[:,1], s=splinesmooth)
            
        baseline_fitted = coeffs(x)
        y_corrected = y-baseline_fitted
        
    elif method == 'gcvspline':
        
        # optional parameters
        splinesmooth = kwargs.get('s',2.0)
        
        # Spline baseline with mode 1 of gcvspl.f, see gcvspline documentation
        c, wk, ier = gcvspline(yafit[:,0],yafit[:,1],np.sqrt(np.abs(yafit[:,1])),splinesmooth,splmode = 1) # gcvspl with mode 1 and smooth factor
        
        baseline_fitted = splderivative(x,yafit[:,0],c)       
        y_corrected = y-baseline_fitted
            
    elif method == 'exp':
        ### Baseline is of the type y = a*exp(b*(x-xo))
        # optional parameters
        p0_exp = kwargs.get('p0_exp',[1.,1.,1.])
        ## fit of the baseline
        coeffs, pcov = curve_fit(funexp,yafit[:,0],yafit[:,1],p0 = p0_exp)
        
        baseline_fitted = funexp(x,coeffs[0],coeffs[1],coeffs[2])
        y_corrected = y-baseline_fitted
    
    elif method == 'log':
        ### Baseline is of the type y = a*exp(b*(x-xo))
        # optional parameters
        p0_log = kwargs.get('p0_log',[1.,1.,1.,1.])
        ## fit of the baseline
        coeffs, pcov = curve_fit(funlog,yafit[:,0],yafit[:,1],p0 = p0_exp)
        
        baseline_fitted = funlog(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3])
        y_corrected = y-baseline_fitted
            
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
        y_corrected = y-baseline_fitted
        
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
        y_corrected = y-baseline_fitted
        
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
        y_corrected = y-baseline_fitted

    return Y_scaler.inverse_transform(y_corrected.reshape(-1, 1)), Y_scaler.inverse_transform(baseline_fitted.reshape(-1, 1))
    #return y_corrected, baseline_fitted
    
    #