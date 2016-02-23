# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:42:21 2014

@author: Charles LE LOSQ
Carnegie Institution of Washington D.C.
February 2014

This is a module with several function for helping importing/treating spectra

"""
import sys
sys.path.append("/Users/charles/Tresors/charles/Software/Rampy/lib-charles/gcvspl/")

import numpy as np
from scipy import interpolate 
from scipy import signal
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import gcvspline


############ SIMPLE MATHEMATICAL FUNCTIONS ###########
def gaussian(x,amp,freq,HWHM): # for spectral fit
    gauss = amp*np.exp(-np.log(2)*((x-freq)/HWHM)**2)
    return gauss
 
def gaussianarea(Amplitude,HWHM,**options): 
    """
    Return the area of a gaussian with inpu amplitude and HWHM.
    Options are eseAmplitude (None by default) and eseHWHM (None by default)
    """
    area = np.sqrt(np.pi/np.log(2))*Amplitude*HWHM
    if options.get("eseAmplitude") != None:
        eseAmplitude = options.get("eseAmplitude")
        if options.get("eseHWHM") != None:
            eseHWHM = options.get("eseHWHM")
            esearea = np.sqrt((np.pi/np.log(2)*HWHM)**2 * eseAmplitude**2 + (np.pi/np.log(2)*Amplitude)**2 * eseHWHM**2)
    else:
        esearea = None
    
    return area, esearea
    
def pseudovoigt(x,amp,freq,HWHM,LGratio): # for spectral fit
    pv1 = LGratio*(amp/(1+((x-freq)/HWHM)**2)) + (1-LGratio)*(amp*np.exp(-np.log(2)*((x-freq)/HWHM)**2))
    return pv1

def funlog(x,a,b,c,d):
    y = a*np.log(-b*(x-c))-d*x**2
    return y
    
def funexp(x,a,b,c):
    y = a*np.exp(b*(x-c))
    return y

def fun2(x,a,b,c):
    y = a + b*x + c*x*x
    return y

def fun1(x,a,b):
    y = a + b*x
    return y
    
def fun(x,a):
    y = a*x
    return y

################## SPECIFIC FUNCTIONS FOR TREATMENT OF SPECTRA

def spectrarray(name,sh,sf,x,interpmethod):
    """The function needs the x general axis, an array containing the name of files, and a precision concenrning the interpolation method
    Spectra are interpolated to have linear and constant x values and be able to be in a same array
    Interpoation methods proposed are the np.interp function or an interpolation using the interpolate 
    function of scipy, based on spline functions
    """
    for i in range(len(name)):
        rawspectre = np.genfromtxt(name[i],skip_header=sh, skip_footer=sf)
        rawspectre = rawspectre[~np.isnan(rawspectre).any(1)] # on verifie qu on a pas de nan
        # on echantillonne le spectre pour avoir les memes x
        if interpmethod == 'np.interp':
            y = np.interp(x,rawspectre[:,0],rawspectre[:,1]) # resampling 
        elif interpmethod == 'scipy.interp':
            tck = interpolate.splrep(rawspectre[:,0],rawspectre[:,1],s=0)
            y = interpolate.splev(x,tck,der=0)
        else:
            print('Error, choose the good interp method name: np.interp or scipy.interp')
        
        # Now we construct the output matrix
        # 1st column is the x axis
        # then others are the spectra in the order provided in the list of names input array
        if i == 0:
            out = np.zeros((len(x),len(name)+1))
            out[:,0]=x
            out[:,i+1]=y
        else:
            out[:,i+1]=y
            
    return out

#### FILTERING OF DATA WITH A BUTTERWORTH FILTER 
   
def spectrafilter(spectre,filtertype,fq,numtaps,columns):
    """
    This function filters the HF of the spectra
    It use a low pass butterworth filter
    y is an array of spectra constructed with the previous spectraarray function
    filtertype contains the string defining which type of filter you want
    Choose between low, high, bandstop, bandpass
    fq is the frequency of the periodic signal you try to erase
    if it is a bandpass or band stop, f must be an array containing the cutoff frequencies
    columns is an array defining which columns you want to treat
    first column of y must contain the frequency
    """
    
    # we already say what is the output array
    out = np.zeros(spectre.shape)
    
    # Butterworth band stop filter caracteristics
    a = spectre[1,0] - spectre[0,0]
    samplerate = 1/a  #Hertz
    nyq_rate = samplerate/2 # frequence Nyquist
    cutf = fq # cutoff frequency
    #bandwidth = 0.005 # largeur filtre, for band pass/stop filters
    numtaps = 1 # ordre du filtre...
    
    for i in range(len(columns)):
        y = spectre[:,columns[i]]
        if (filtertype == 'low') or (filtertype == 'high'):
            b, a = signal.butter(numtaps, [(cutf/nyq_rate)], btype = filtertype)
            out[:,columns[i]] = signal.filtfilt(b, a, y) # filter with phase shift correction
        else:
            b, a = signal.butter(numtaps, [(cutf[0]/nyq_rate),(cutf[1]/nyq_rate)], btype = filtertype)
            out[:,columns[i]] = signal.filtfilt(b, a, y) # filter with phase shift correction

    # Note forgetting to register the x axis        
    out[:,0] = spectre[:,0]
    
    return out

#### OFFSETTING DATA FROM THE SPECTRARRAY FUNCTION

def spectraoffset(spectre,offsets):
    """
    This function allows to Y offset your spectra for representing them or correcting them
    with an horizontal baseline for instance
    spectre is an array constructed with the spectrarray function
    offsets is an array constructed with numpy and containing the coefficient for the offset to apply to spectre   
    """
    
    # we already say what is the output array
    out = np.zeros(spectre.shape)
    # and we offset the spectra
    for i in range(len(offsets)):
        out[:,i+1] = spectre[:,i+1] + offsets[i]
    
    return out
    
#### SPLINES, LINEAR OR POLYNOMIAL BASELINES    
    
def linbaseline(spectre,bir,method,splinesmooth):
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

def smooth(x,window_len=11,window="hamming"):
    """smooth the data using a window with requested size.
    
    Taken from http://wiki.scipy.org/Cookbook/SignalSmooth
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)] #Not working pretty well....

def spectrataux(x,spectres):
    """
    Created on Wed Feb 12 18:30:44 2014
    
    @author: Charles LE LOSQ
    Carnegie Institution of Washington D.C.
    February 2014
    This function calculates the increase/decrease rate of each frequency in a spectrum
    Your data must have the same frequency axis of course...
    Initial fitting function is a second order polynomial
    Change that following your needs
    """
    # we need an organized function before calling the curve_fit algorithm
    freq = spectres[:,0]
    # output array
    taux = np.zeros((len(freq),4));
    taux[:,0] = freq[:]
    
    # We look a each frequency, we sort y data and fit them with a second order polynomial
    for i in range(len(freq)):
        y = spectres[i,1::]
        popt, pcov = curve_fit(fun2,x,y,[0.5e-3,0.5e-4,1e-6])
        taux[i,1:len(x)]=popt
        
    return taux
            
def longcorr(data,temp,wave): # input are a two column matrix of data, temperature in C, wave as the wavelength in cm-1
    """
    # Long Correction
    # Charles Le Losq
    # CIW Washington 2014
    
    # Long's correction of Raman spectra and normalisation
    # last rev. Oct 2010, converted to Matlab and then to Python
    # ensures strictly increasing values of wavenumber
    # calc. e.s.e. as Long cor. norm. sqrt(n_raw) 3d output col.
    # exp. format to avoid null ese.
    # program long3;
    
    #Original program in pascal from J. Roux, modified for Python C. Le Losq.
    
    # See Shucker and Gammon, Phys. Rev. Lett. 1970; Galeener and Sen, Phys. Rev. B 1978; Neuville and Mysen, GCA 1996; Le Losq et al. AM 2012 and GCA 2014 for equations and theory
    """
    h = 6.62606896*10**-34   # J.s    Plank constant
    k = 1.38066e-23;     # J/K    Boltzman
    c = 2.9979e8;        # m/s    Speed of light
    v = wave;            # cm-1   Excitating laser line
    nu0 = 1.0/v*10**9;    # conversion of v in m
    T = temp + 273.15;   # C input temperature in K

    x = data[:,0]
    y = data[:,1]
    
    
    # Calculate the error on data as sqrt(y). If y <= 0, then error = sqrt of absolute value of y.
    # for doing that we construct an artificial array y containing only positive values 
    
    ese = np.sqrt(np.absolute(data[:,1])) # we assume that errors on raw data are in sqrt(y)
    
    # For retrieving the good errors after long correction, one simple way is
    # to work with relative errors...
    error = ese/y;
    
    # then we proceed to the correction (Neuville and Mysen, 1996; Le Losq et
    # al., 2012)
    
    nu = 100.0*x; # cm-1 -> m-1 Raman shift
    rnu = nu0-nu; # nu0 is in m-1
    t0 = nu0*nu0*nu0*nu/rnu/rnu/rnu/rnu;
    t1 = -h*c*nu/k/T; # c in m/s  : t1 dimensionless
    t2 = 1 - np.exp(t1);
    longsp = y*t0*t2; # for y values
    #long2 = ese*t0*t2; # for errors, as comment as we use relative errors
    
    
    # normalized to max intensity
    # tried area with np.trapz but their is an issue for now
    norm = np.trapz(longsp,x)
    #norm = np.max(longsp)
    longsp = longsp/norm
    
    eselong = error*longsp
    
    spectreout = np.zeros((len(x),3))
    spectreout[:,0] = x
    spectreout[:,1] = longsp
    spectreout[:,2] = eselong
    
    return spectreout
  
 
########### SPECIFIC GAUSSIAN/VOIGTR FUNCTIONS FOR USING WITH SCIPY OPTIMIZATION PROTOCOLS
   
def multigaussian(x,params):
    taille = len(params)
    y = np.zeros(len(x),4)
    for i in range(taille):
        y[:,taille+1] = params[taille,0]*np.exp(-np.log(2)*((x-params[taille,1])/params[taille,2])**2)
    y[:,0] = y[:,1]+y[:,2]+y[:,3]
    return y
    

def gauss_lsq(params,x): 
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
    
