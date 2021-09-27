#!/usr/bin/env python
#-*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

import rampy as rp
from rampy import peak_shapes

class map():
    """treat maps of Raman spectra

    Parameters
    ----------
    file : str
        filename, including path

    spectrometer_type : str
        type of spectrometer, choose between "horiba" or "renishaw"
    """

    def __init__(self,file, spectrometer_type = "horiba"):
        self.spectrometer_type = spectrometer_type
        self.filename = file

        if spectrometer_type == "horiba":
            self.X, self.Y, self.w, self.I = read_horiba(file)
        elif spectrometer_type == "renishaw":
            self.X, self.Y, self.w, self.I = read_renishaw(file)
        else:
            print("Only renishaw or horiba are supported at the moment.")

        # here we check that the frequency axis w is increasing only:
        if self.w[0]>self.w[-1]:
            self.w = np.flip(self.w)
            self.I = np.flip(self.I,axis=0)

    def background(self, bir, method = "poly", **kwargs):
        """correct a background from the initial signal I on a map using rampy.baseline

        Parameters
        ----------
        bir : ndarray
            arrays of the backgroudn interpolation regions.
        method : string
            see rampy.baseline documentation for methods available. Default is polynomial

        All kwargs argument for rampy.baseline() will be forwarded and can be used there.

        Returns
        -------
        Background and corrected spectra area available at self.background and self.I_corrected
        """
        self.I_background = np.copy(self.I)
        self.I_corrected = np.copy(self.I)
        for i in range(len(self.X)):
            y_, bkg_ = rp.baseline(self.w, self.I[:,i], bir, method, **kwargs)
            self.I_corrected[:,i] = y_.ravel()
            self.I_background[:,i] = bkg_.ravel()

    def normalise(self, y, method="intensity"):
        """normalise the spectra to their max intensity, their area or min-max normalisation

        This uses the internals of rampy.normalise.

        Parameters
        ----------
        y : object intensities
            the intensities to normalise. For instance, if you want to normalised the background corrected I, pass self.I_corrected.
        method : string
            method used, choose between area, intensity, minmax

        Returns
        -------
        The normalised spectra are available at self.I_normalised
        """
        self.I_normalised = np.copy(self.I)
        for i in range(len(self.X)):
            y_ = rp.normalise(y[:,i].ravel(),x=self.w,method=method)
            self.I_normalised[:,i] = y_.ravel()

    def smooth(self, y, method="whittaker",**kwargs):
        """uses the smooth function of rampy to smooth the signals
        Parameters
        ----------
        y : object intensities
            the intensities to normalise. For instance, if you want to normalised the background corrected I, pass self.I_corrected.
        method : str
            Method for smoothing the signal;
            choose between savgol (Savitzky-Golay), GCVSmoothedNSpline, MSESmoothedNSpline, DOFSmoothedNSpline, whittaker, 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'.

        kwargs
        ======
        window_length : int
            The length of the filter window (i.e. the number of coefficients). window_length must be a positive odd integer.
        polyorder : int
            The order of the polynomial used to fit the samples. polyorder must be less than window_length.
        Lambda : float
            smoothing parameter of the Whittaker filter described in Eilers (2003). The higher the smoother the fit.
        d : int
            d parameter in Whittaker filter, see Eilers (2003).
        ese_y : ndarray
            errors associated with y (for the gcvspline algorithms)

        Returns
        -------
        self.y_smoothed : ndarray
            the smoothed signal for the map

        """
        self.y_smoothed = np.copy(self.I)
        for i in range(len(self.X)):
            y_ = rp.smooth(self.w, y[:,i],method=method, **kwargs)
            self.y_smoothed[:,i] = y_.ravel()

    def centroid(self, y, region_to_investigate):
        """calculate the centroid in a given region of interest

        Parameters
        ----------
        y : object intensities
            the intensities to normalise. For instance, pass self.normalised for performing the calculation on normalised spectra.
        region_to_investigate : 1x2 array
            the x values of regions where the centroid will be measured

        Returns
        -------
        self.centroid_position : ndarray
	    	centroid position for the map
        """
        self.centroid_position = np.copy(self.X)
        for i in range(len(self.X)):
            sp_ = rp.get_portion_interest(self.w, y[:,i], region_to_investigate)
            self.centroid_position[i] = rp.centroid(sp_[:,0], sp_[:,1])

    def intensity(self, y, region_to_investigate):
        """get the maximum intensity in the region to investigate.

        The intensity maximum is estimated from a simple np.max() search.
        Do not forget to smooth the signal if necessary prior to using this.
        
        Parameters
        ----------
        y : object intensities
            the intensities to consider. For instance, pass self.normalised for performing the calculation on normalised spectra.
        region_to_investigate : 1x2 array
            the x values of regions where the intensity will be measured

        Returns
        -------
        self.I_max : ndarray
            Intensity maximum
        """
        self.I_max = np.copy(self.X)
        for i in range(len(self.X)):
            sp_ = rp.get_portion_interest(self.w, y[:,i], region_to_investigate)
            self.I_max[i] = np.max(sp_[:,1])
            
    def area(self, y, region_to_investigate):
        """get the area under the curve in the region to investigate.

        The area is calculated by trapezoidal integration, using np.trapz()        
        Do not forget to smooth the signal if necessary prior to using this.
        
        Parameters
        ----------
        y : object intensities
            the intensities to consider. For instance, pass self.normalised for performing the calculation on normalised spectra.
        region_to_investigate : 1x2 array
            the x values of regions where the area will be measured

        Returns
        -------
        self.A : ndarray
                maximum to make a nice plot
        """
        self.A = np.copy(self.X)
        for i in range(len(self.X)):
            sp_ = rp.get_portion_interest(self.w, y[:,i], region_to_investigate)
            self.A[i] = np.trapz(sp_[:,1],sp_[:,0])

    def intensity_ratio(self, y, region_to_investigate):
        """get the intensity ratio between two regions of interest.

        The intensity maxima are estimated from a simple np.max() search.
        Do not forget to smooth the signals if necessary prior to using this.
        Parameters
        ----------
        y : object intensities
            the intensities to consider. For instance, pass self.normalised for performing the calculation on normalised spectra.
        region_to_investigate : 2x2 array
            the x values of regions where the intensity ratios will be measured. The two lines record the two regions of interest.

        Returns
        -------
        self.I_ratio : ndarray
			Intensity ratio
        """
        self.I_ratio = np.copy(self.X)
        I_max1 = np.copy(self.X)
        I_max2 = np.copy(self.X)
        for i in range(len(self.X)):
            sp_1 = rp.get_portion_interest(self.w, y[:,i], region_to_investigate[0,:].reshape(1,2))
            sp_2 = rp.get_portion_interest(self.w, y[:,i], region_to_investigate[1,:].reshape(1,2))
            I_max1 = np.max(sp_1[:,1])
            I_max2 = np.max(sp_2[:,1])
            self.I_ratio[i] = I_max1/I_max2
            
            
    def area_ratio(self, y, region_to_investigate):
        """get the area ratio between two regions of interest.

        The areas are calculated by trapezoidal integration, using np.trapz()        
        Do not forget to smooth the signals if necessary prior to using this.
        
        Parameters
        ----------
        y : object intensities
            the intensities to consider. For instance, pass self.normalised for performing the calculation on normalised spectra.
        region_to_investigate : 1x2 array
            the x values of regions where the areas will be measured. The two lines record the two regions of interest.

        Returns
        -------
        self.A_ratio : ndarray
			Area ratio = area region 1 / area region 2
        """
        self.A_ratio = np.copy(self.X)
        A_max1 = np.copy(self.X)
        A_max2 = np.copy(self.X)
        for i in range(len(self.X)):
            sp_1 = rp.get_portion_interest(self.w, y[:,i], region_to_investigate[0,:].reshape(1,2))
            sp_2 = rp.get_portion_interest(self.w, y[:,i], region_to_investigate[1,:].reshape(1,2))
            A_max1 = np.trapz(sp_1[:,1],sp_1[:,0])
            A_max2 = np.trapz(sp_2[:,1],sp_2[:,0])
            self.A_ratio[i] = A_max1/A_max2

def read_renishaw(file):
    """read Renishaw csv maps

    Parameters
    ==========
    file : str
        filename, including path

    Returns
    -------
    X : m by n array
        X positions
    Y : m by n array
        Y position
    lambdas : m by n array
        Raman shift
    intensities : m by n array
        Intensities
    """
    df=pd.read_csv(file,names=['x','y','lambda','int'],sep='\t')
    table=df.loc[:,'int'].values
    lambdas=df.loc[:,'lambda'].values
    lambda_0=lambdas[0]
    lambdas_one=lambdas[: (np.argwhere(lambdas==lambda_0)[1])[0]]

    X=df.iloc[::lambdas_one.shape[0],0].values
    Y=df.iloc[::lambdas_one.shape[0],1].values

    intensities=np.transpose(np.reshape(table,(X.shape[0],lambdas_one.shape[0])))
    lambdas=np.transpose(np.reshape(lambdas,(X.shape[0],lambdas_one.shape[0])))

    return X, Y, lambdas_one, intensities

def read_horiba(file):
    """read Horiba csv maps

    Parameters
    ==========
    file : str
        filename, including path

    Returns
    -------
    X : m by n array
        X positions
    Y : m by n array
        Y position
    lambdas : m by n array
        Raman shift
    intensities : m by n array
        Intensities
    """

    df = pd.read_csv(file,sep='\t')
    intensities = df.iloc[:,2:].values
    lambdas = df.columns.values[2:].values.astype(float)
    X = df.iloc[:,0].values
    Y = df.iloc[:,1].values
    return X, Y, lambdas, intensities.T

def peak(X, Y, lambdas, intensities, function, Xrange, amp, Xmean, sigma, y0, A):
    """
    #fitting
    """
    if function=='gauss':
        fun=peak_shapes.create_gauss()
    elif function=='lorenz':
        fun=peak_shapes.create_lorenz()
    results=np.empty(5)
    for d in np.ndindex(intensities.shape[1]):
        try:
            popt, pcov = curve_fit(fun, lambdas[Xrange[0]:Xrange[1]],
                                   np.squeeze(intensities[Xrange[0]:Xrange[1],d]),
                                   p0=(amp,Xmean,sigma,y0,A))
        except RuntimeError:
            print("Error - curve_fit failed")
        results=np.vstack((results,popt))
    #maps
    if X[0]!=X[1]:
        n_X0=np.argwhere(X!=X[0])[0,0] # while main axis in x
        n_X1=int(X.shape[0]/n_X0)
    else:
        n_X1=np.argwhere(Y!=Y[0])[0,0] # while main axis in y
        n_X0=int(Y.shape[0]/n_X1)


    rmap=np.empty([n_X0,n_X1])

    for d in np.ndindex(results.shape[1]):
        rmap=np.dstack((rmap,results[1:,d].reshape(n_X0,n_X1)))

    return results, rmap
