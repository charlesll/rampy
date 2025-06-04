# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

import rampy as rp
from rampy import peak_shapes

class maps():
    """Class for handling and analyzing Raman spectral maps.

    This class provides methods for reading, processing, and analyzing
    1D or 2D Raman spectral maps from Horiba or Renishaw spectrometers.

    Args:
        file (str): Path to the map file.
        spectrometer_type (str, optional): Type of spectrometer. 
            Must be either 'horiba' or 'renishaw'. Defaults to 'horiba'.
        map_type (str, optional): Type of map. Must be either '2D' or '1D'. Defaults to '2D'.
        despiking (bool, optional): Whether to perform despiking on the spectra. Defaults to False.
        neigh (int, optional): Number of neighboring points for despiking. Defaults to 4.
        threshold (float, optional): Threshold for spike detection. Defaults to 3.

    Raises:
        ValueError: If `map_type` is not '1D' or '2D'.
    """

    def __init__(self, file, 
                 spectrometer_type = "horiba", 
                 map_type = "2D", 
                 despiking=False, 
                 neigh=4, 
                 threshold = 3):
        self.spectrometer_type = spectrometer_type
        self.filename = file

        if map_type == "2D":
            if spectrometer_type == "horiba":
                self.X, self.Y, self.w, self.I = read_horiba(file, map_type = "2D")
            elif spectrometer_type == "renishaw":
                self.X, self.Y, self.w, self.I = read_renishaw(file)
            else:
                print("Only renishaw or horiba are supported for 2D maps at the moment.")
        elif map_type == "1D":
            if spectrometer_type == "horiba":
                self.X, self.w, self.I = read_horiba(file, map_type = "1D")
            else:
                print("Only horiba line maps are supported at the moment.")
        else:
            raise ValueError("map_type should be set to either '1D' or '2D'")

        # here we check that the frequency axis w is increasing only:
        if self.w[0]>self.w[-1]:
            self.w = np.flip(self.w)
            self.I = np.flip(self.I,axis=0)

        # here we remove spikes if the user wants it
        if despiking == True:
            for i in range(self.I.shape[1]):
                self.I = despiking(self.w, self.I[:,i], neigh=neigh, threshold=threshold)

    def background(self, bir, method = "poly", **kwargs):
        """Corrects background from the signal using rampy.baseline.

        Applies baseline correction to each spectrum in the map using the specified method.

        Args:
            bir (ndarray): Regions of interest for baseline fitting.
            method (str, optional): Baseline correction method. Defaults to 'poly'.
                Supported methods include:
                    - 'poly': Polynomial fitting.
                    - 'unispline': Spline fitting (UnivariateSpline).
                    - 'gcvspline': GCVSpline fitting.
                    - 'gaussian': Gaussian function fitting.
                    - 'exp': Exponential fitting.
                    - 'log': Logarithmic fitting.
                    - 'rubberband': Rubberband correction.
                    - 'als': Asymmetric Least Squares.
                    - 'arPLS': Asymmetrically Reweighted Penalized Least Squares.
                    - 'drPLS': Doubly Reweighted Penalized Least Squares.
                    - 'whittaker': Whittaker smoothing.
                    - 'GP': Gaussian process method.
            **kwargs: Additional method-specific parameters, such as:
                roi (ndarray): Regions of interest (default: full range).
                polynomial_order (int): Polynomial degree for 'poly'.
                s (float): Spline smoothing parameter.
                ese_y (ndarray): Errors for 'gcvspline'.
                lam (float): Smoothness parameter for 'als', 'arPLS', etc.
                p (float): Weighting parameter for 'als'.
                ratio (float): Convergence ratio for 'arPLS', 'drPLS'.
                niter (int): Number of iterations for 'als', 'drPLS'.
                eta (float): Roughness for 'drPLS'.
                p0_gaussian (list): Initial params for 'gaussian'.
                p0_exp (list): Initial params for 'exp'.
                p0_log (list): Initial params for 'log'.

        Returns:
            None. Sets `self.I_background` and `self.I_corrected` with the background and background-corrected spectra.
        """
        self.I_background = np.copy(self.I)
        self.I_corrected = np.copy(self.I)
        for i in range(len(self.X)):
            y_, bkg_ = rp.baseline(self.w, self.I[:,i], method, roi=bir, **kwargs)
            self.I_corrected[:,i] = y_.ravel()
            self.I_background[:,i] = bkg_.ravel()

    def normalise(self, y, method="intensity"):
        """Normalises the spectra using rampy.normalise.

        Args:
            y (ndarray): Intensities to normalise (e.g., `self.I_corrected`).
            method (str, optional): Normalisation method. One of 'area', 'intensity', or 'minmax'. Defaults to 'intensity'.

        Returns:
            None. Sets `self.I_normalised` with the normalised spectra.
        """
        self.I_normalised = np.copy(self.I)
        for i in range(len(self.X)):
            y_ = rp.normalise(y[:,i].ravel(),x=self.w,method=method)
            self.I_normalised[:,i] = y_.ravel()

    def smooth(self, y, method="whittaker",**kwargs):
        """Smooths spectra using rampy.smooth.

        Args:
            y (ndarray): Intensities to smooth (e.g., `self.I_corrected`).
            method (str, optional): Smoothing method. Supported methods include:
                'savgol', 'GCVSmoothedNSpline', 'MSESmoothedNSpline', 'DOFSmoothedNSpline',
                'whittaker', 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'.
                Defaults to 'whittaker'.
            **kwargs: Additional method-specific parameters, such as:
                window_length (int): Length of the filter window (for Savitzky-Golay).
                polyorder (int): Polynomial order (for Savitzky-Golay).
                Lambda (float): Smoothing parameter for Whittaker filter.
                d (int): Parameter for Whittaker filter.
                ese_y (ndarray): Errors for spline algorithms.

        Returns:
            None. Sets `self.y_smoothed` with the smoothed spectra.
        """
        self.y_smoothed = np.copy(self.I)
        for i in range(len(self.X)):
            y_ = rp.smooth(self.w, y[:,i], method=method, **kwargs)
            self.y_smoothed[:,i] = y_.ravel()

    def centroid(self, y, region_to_investigate):
        """Calculates the centroid in a given region of interest.

        Args:
            y (ndarray): Intensities to analyze (e.g., `self.I_normalised`).
            region_to_investigate (array-like): [min, max] values of the region where the centroid is measured.

        Returns:
            None. Sets `self.centroid_position` with the centroid positions for the map.
        """
        self.centroid_position = np.copy(self.X)
        for i in range(len(self.X)):
            sp_ = rp.extract_signal(self.w, y[:,i], region_to_investigate)
            self.centroid_position[i] = rp.centroid(sp_[:,0], sp_[:,1])

    def intensity(self, y, region_to_investigate):
        """Finds the maximum intensity in a specified region.

        The maximum is estimated using `np.max` within the region.

        Args:
            y (ndarray): Intensities to analyze (e.g., `self.I_normalised`).
            region_to_investigate (array-like): [min, max] values of the region to search for the maximum.

        Returns:
            None. Sets `self.I_max` with the intensity maxima for the map.
        """
        self.I_max = np.copy(self.X)
        for i in range(len(self.X)):
            sp_ = rp.extract_signal(self.w, y[:,i], region_to_investigate)
            self.I_max[i] = np.max(sp_[:,1])

    def area(self, y, region_to_investigate):
        """Calculates the area under the curve in a specified region.

        The area is computed via trapezoidal integration (`np.trapz`).

        Args:
            y (ndarray): Intensities to analyze (e.g., `self.I_normalised`).
            region_to_investigate (array-like): [min, max] values of the region to integrate.

        Returns:
            None. Sets `self.A` with the integrated areas for the map.
        """
        self.A = np.copy(self.X)
        for i in range(len(self.X)):
            sp_ = rp.extract_signal(self.w, y[:,i], region_to_investigate)
            self.A[i] = np.trapz(sp_[:,1],sp_[:,0])

    def intensity_ratio(self, y, region_to_investigate):
        """Calculates the intensity ratio between two regions of interest.

        The ratio is computed as the maximum intensity in region 1 divided by the maximum in region 2.

        Args:
            y (ndarray): Intensities to analyze (e.g., `self.I_normalised`).
            region_to_investigate (ndarray): 2x2 array, each row [min, max] for the two regions.

        Returns:
            None. Sets `self.I_ratio` with the intensity ratios for the map.
        """
        self.I_ratio = np.copy(self.X)
        I_max1 = np.copy(self.X)
        I_max2 = np.copy(self.X)
        for i in range(len(self.X)):
            sp_1 = rp.extract_signal(self.w, y[:,i], region_to_investigate[0,:].reshape(1,2))
            sp_2 = rp.extract_signal(self.w, y[:,i], region_to_investigate[1,:].reshape(1,2))
            I_max1 = np.max(sp_1[:,1])
            I_max2 = np.max(sp_2[:,1])
            self.I_ratio[i] = I_max1/I_max2


    def area_ratio(self, y, region_to_investigate):
        """Calculates the area ratio between two regions of interest.

        The ratio is computed as the integrated area in region 1 divided by the area in region 2.

        Args:
            y (ndarray): Intensities to analyze (e.g., `self.I_normalised`).
            region_to_investigate (ndarray): 2x2 array, each row [min, max] for the two regions.

        Returns:
            None. Sets `self.A_ratio` with the area ratios for the map.
        """
        self.A_ratio = np.copy(self.X)
        A_max1 = np.copy(self.X)
        A_max2 = np.copy(self.X)
        for i in range(len(self.X)):
            sp_1 = rp.extract_signal(self.w, y[:,i], region_to_investigate[0,:].reshape(1,2))
            sp_2 = rp.extract_signal(self.w, y[:,i], region_to_investigate[1,:].reshape(1,2))
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

def read_horiba(file, map_type="2D"):
    """read Horiba csv maps (1D, 2D)
    
    Parameters
    ==========
    file : str
        filename, including path
    map_type : str
        1D map (line) or 2D map, default: 2D
        
    Returns
    -------
    X : m by n array
        X positions
    Y : m by n array
        Y position
    lambdas : n array
        Raman shift
    intensities : m by n array
        Intensities
    """

    df = pd.read_csv(file,sep='\t')
    
    if map_type == "2D":
        intensities = df.iloc[:,2:].values
        lambdas = df.columns.values[2:].astype(float)
        X = df.iloc[:,0].values
        Y = df.iloc[:,1].values
        return X, Y, lambdas, intensities.T
    elif map_type == "1D":
        intensities = df.iloc[:,1:].values
        lambdas = df.columns.values[1:].astype(float)
        X = df.iloc[:,0].values
        return X, lambdas, intensities.T
    else:
        raise ValueError("Not implemented, set map_type to '1D' or '2D'")

def peak(X, Y, lambdas, intensities, function, Xrange, amp, Xmean, sigma, y0, A):
    """to fit peaks in a map. Work in progress.
    
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
