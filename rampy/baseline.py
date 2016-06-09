from dependencies import gcvspline
import numpy as np
from scipy import interpolate 
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

def baseline(spectre,bir,method,splinesmooth):
    """
    This function allows subtracting a baseline under the spectra
    spectre is a spectrum or an array of spectra constructed with the spectrarray function
    bir contains the Background Interpolation Regions, it must be a n x 2 dimensiona rray
    
    methods:
    "linear": linear baseline, with spectre = array[x y];
    "hori': constant baseline, fitted at the minimum in the provided region of spectra. Splinesmooth in this case is the 1/2 extent of the region where the mean minimum is calculated;
    "unispline": spline with the UnivariateSpline function of Scipy, splinesmooth is the spline smoothing factor (assume equal weight in the present case);
    "gcvspline": spline with the gcvspl.f algorythm, really robust. Spectra must have x, y, ese in it, and splinesmooth is the smoothing factor;
    for gcvspline, if ese are not provided we assume ese = sqrt(y);
    "poly": polynomial fitting, with splinesmooth the degree of the polynomial.
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
        coeffs = None        
        
    elif method == 'unispline':
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
    elif method == 'gcvspline':
        ## WARNING THEIR IS THE ERROR HERE IN THE SPECTRE MATRIX, not the case for the other functions
        ## ONLY TREAT ONE SPECTRA AT A TIME       
        ### selection of bir data
        for i in range(birlen):
            if i == 0:
                yafit = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))] 
            else:
                je = spectre[np.where((spectre[:,0]> bir[i,0]) & (spectre[:,0] < bir[i,1]))]
                yafit = np.concatenate((yafit,je),axis=0)
        
        # Spline baseline with mode 3 of gcvspl.f
        xdata = yafit[:,0]
        ydata = np.zeros((len(xdata),1))
        ydata[:,0] = yafit[:,1]
        test = np.shape(yafit)
        if test[1] > 2:
            ese = yafit[:,2]
        else:
            ese = np.sqrt(np.abs(yafit[:,1]))
        VAL = ese**2
        c, wk, ier = gcvspline.gcvspline(xdata,ydata,splinesmooth*ese,VAL,splmode = 3) # gcvspl with mode 3 and smooth factor
        out2[:,1] = gcvspline.splderivative(x,xdata,c)       
        out1[:,1] = spectre[:,1]-out2[:,1]
        coeffs = None
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

