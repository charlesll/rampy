import numpy as np
from gcvspline import gcvspline, splderivative 
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.spatial import ConvexHull
import scipy.sparse as sparse
from numpy.linalg import norm

def baseline(spectre,bir,method, **kwargs):
    """
    This function allows subtracting a baseline under the spectra
    spectre is a spectrum or an array of spectra constructed with the spectrarray function
    bir contains the Background Interpolation Regions, it must be a n x 2 dimensiona rray
    
    Inputs
    ------
    
        Spectre: Array with 2 or more columns. First column contain x axis, subsequent columns contains y values. If using gcvspline, only treat a spectrum at a call, third column can contain known ese.
    
        bir: an Array containing the regions of interest, organised per line. for instance, roi = np.array([[100., 200.],[500.,600.]]) will define roi between 100 and 200 as well as between 500 and 600,.
    
        methods:
    
    "linear": linear baseline, with spectre = array[x y];
    "hori': constant baseline, fitted at the minimum in the provided region of spectra. Splinesmooth in this case is the 1/2 extent of the region where the mean minimum is calculated;
    "unispline": spline with the UnivariateSpline function of Scipy, splinesmooth is the spline smoothing factor (assume equal weight in the present case);
    "gcvspline": spline with the gcvspl.f algorythm, really robust. Spectra must have x, y, ese in it, and splinesmooth is the smoothing factor;
    for gcvspline, if ese are not provided we assume ese = sqrt(y);
    "poly": polynomial fitting, with splinesmooth the degree of the polynomial.
    "rubberband": rubberband baseline fitting
    "als": automatic least square fitting following Eilers and Boelens 2005.
    "arPLS": automatic baseline fit using the algorithm from Baek et al. 2015 Baseline correction using asymmetrically reweighted penalized least suqares smoothing, Analyst 140: 250-257.
    
    Options
    -------
    
        Options for the ALS algorithm are:
    
        lambda: float, between 10**2 to 10**9, default = 10**5
        p: float, advised value between 0.001 to 0.1, default = 0.01
        niter: number of iteration, default = 10
    
        Options for the ALS algorithm are:
    
        lam: float, the lambda smoothness parameter. Typical values are between 10**2 to 10**9, default = 10**5
        ratio = float, the termination ratio. Set it high to start and lower it. Default = 0.01.
    
        Options for the unispline and gcvspline algorithms are:
    
        s: smoothing coefficient, the higher the more smooth.
    
    Outputs
    -------
    
        out1: an 2 columns x-y array containing the corrected signal
    
        out2: an 2 columns x-y array containing the baseline
    
        coefs: contains spline coefficients.
    
    """
    
    # we already say what is the output array
    out1 = np.zeros(spectre.shape) # matrix for corrected spectra
    out2 = np.zeros(spectre.shape) # matrix with the baselines    
    x = spectre[:,0] # x axis
    out1[:,0] = x[:]
    out2[:,0] = x[:]
 
    nbsp = np.array(spectre.shape[1]-1)
    birlen = np.array(bir.shape[0])
    
    if method == 'linear':
        
        taux = np.zeros([nbsp,bir.shape[1]]) # thrid output matrix with coefficients
           
        ### selection of bir data
        for i in range(birlen):
            if i == 0:
                yafit = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))] 
            else:
                je = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))]
                yafit = np.concatenate((yafit,je),axis=0)
        
        # fit of the baseline
        for i in range(nbsp):
            popt, pcov = curve_fit(fun1,yafit[:,0],yafit[:,i+1],p0=[1,1])
            taux[i,:]=popt
            out1[:,i+1] = spectre[:,i+1] - (popt[0] + popt[1]*x)
            out2[:,i+1] = (popt[0] + popt[1]*x)
            
        return out1, out2, taux
        
    elif method =='hori':
        # we take the data in the region of interest        
        yafit = spectre[np.where((spectre[:,0]> bir[0,0]) & (spectre[:,0] < bir[0,1]))] 
        
        idx = np.where(yafit[:,1] == np.min(yafit[:,1])) #we search where the minimum is in yafit
        yafit2 = yafit[np.where((yafit[:,0] > (yafit[idx[0],0]-splinesmooth)) & (yafit[:,0] < (yafit[idx[0],0]+splinesmooth)))]
        constant = np.mean(yafit2[:,1]) #and we calculate the mean horizontal baseline in the region with an extent defined by splinesmooth
        
        out1[:,1] = spectre[:,1] - constant
        out2[:,1] = constant
        
        return out1, out2
        
    elif method == 'unispline':
        
        # optional parameters
        splinesmooth = kwargs.get('s',2.0)
        
        ### selection of bir data
        for i in range(birlen):
            if i == 0:
                yafit = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))] 
            else:
                je = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))]
                yafit = np.concatenate((yafit,je),axis=0)
        # fit of the baseline
        for i in range(nbsp):
            coeffs = UnivariateSpline(yafit[:,0],yafit[:,i+1], s=splinesmooth)
            out2[:,i+1] = coeffs(x)
            out1[:,i+1] = spectre[:,i+1]-out2[:,i+1]
        
        return out1, out2
        
    elif method == 'gcvspline':
        ## WARNING THEIR IS THE ERROR HERE IN THE SPECTRE MATRIX, not the case for the other functions
        ## ONLY TREAT ONE SPECTRA AT A TIME       
        
        # optional parameters
        splinesmooth = kwargs.get('s',2.0)
        
        ### selection of bir data
        for i in range(birlen):
            if i == 0:
                yafit = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))] 
            else:
                je = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))]
                yafit = np.concatenate((yafit,je),axis=0)
        
        # Spline baseline with mode 1 of gcvspl.f, see gcvspline documentation
        xdata = yafit[:,0]
        ydata = np.zeros((len(xdata),1))
        ydata[:,0] = yafit[:,1]
        test = np.shape(yafit)
        if test[1] > 2:
            ese = yafit[:,2]
        else:
            ese = np.sqrt(np.abs(yafit[:,1]))
        VAL = 1.0
        c, wk, ier = gcvspline(xdata,ydata,ese,splinesmooth,splmode = 1) # gcvspl with mode 1 and smooth factor
        out2[:,1] = splderivative(x,xdata,c)       
        out1[:,1] = spectre[:,1]-out2[:,1]
        
        return out1, out2
        
    elif method == 'poly':
        ## Polynomial baseline
        ### selection of bir data
        for i in range(birlen):
            if i == 0:
                yafit = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))] 
            else:
                je = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))]
                yafit = np.concatenate((yafit,je),axis=0)
        # fit of the baseline
        # Note that splineorder will serve for contraining also the polynomial order
        for i in range(nbsp):
            coeffs = np.polyfit(yafit[:,0],yafit[:,i+1],splinesmooth)
            out2[:,i+1] = np.polyval(coeffs,x)
            out1[:,i+1] = spectre[:,i+1]-out2[:,i+1]
            
        return out1, out2, coeffs
            
    elif method == 'log':
        ### Baseline is of the type y = a*exp(b*(x-xo))
        ### selection of bir data
        for i in range(birlen):
            if i == 0:
                yafit = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))] 
            else:
                je = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))]
                yafit = np.concatenate((yafit,je),axis=0)
        ## fit of the baseline
        for i in range(nbsp):
            coeffs, pcov = curve_fit(funlog,yafit[:,0],yafit[:,i+1],p0 = splinesmooth)
            out1[:,i+1] = spectre[:,i+1] - funlog(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3])
            out2[:,i+1] = funlog(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3])
            
        return out1, out2, coeffs
            
    elif method == 'rubberband':
        # code from this stack-exchange forum
        #https://dsp.stackexchange.com/questions/2725/how-to-perform-a-rubberband-correction-on-spectroscopic-data
        
        x = spectre[:,0]
        y = spectre[:,1]
        
        # Find the convex hull
        v = ConvexHull(np.array([x, y])).vertices
            
        # Rotate convex hull vertices until they start from the lowest one
        v = np.roll(v, -v.argmin())
        # Leave only the ascending part
        v = v[:v.argmax()]

        # Create baseline using linear interpolation between vertices
        out2[:,1] = np.interp(x, x[v], y[v])
        out1[:,1] = y - out2[:,1]
        
        return out1, out2
        
    elif method == 'als':
        # Matlab code in Eilers et Boelens 2005
        # Python addaptation found on stackoverflow: https://stackoverflow.com/questions/29156532/python-baseline-correction-library
        
        # setting y
        y = spectre[:,1]
        
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
        
        out2[:,1] = z
        out1[:,1] = y - z
        
        return out1, out2
        
    elif method == 'arPLS':
        # Adaptation of the Matlab code in Baek et al 2015 
    
        # setting y
        y = spectre[:,1]
        
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

        out2[:,1] = z
        out1[:,1] = y-z

        return out1, out2
