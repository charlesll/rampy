# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simpson

import rampy

def spectrarray(name: np.ndarray, sh: int, sf: int, x: np.ndarray) -> np.ndarray:
    """Constructs a unified array of spectra with a common x-axis.

    This function reads spectral data from multiple files, resamples the y-values to match 
    a common x-axis, and combines them into a single array.

    Args:
        name (np.ndarray): Array of file names containing the spectra.
        sh (int): Number of header lines to skip in each file.
        sf (int): Number of footer lines to skip in each file.
        x (np.ndarray): The common x-axis values to which all spectra will be resampled.

    Returns:
        np.ndarray: A 2D array where the first column contains the common x-axis values 
        and subsequent columns contain the resampled y-values for each spectrum.

    Raises:
        ValueError: If any file contains invalid or missing data.

    Notes:
        - The function uses `np.genfromtxt` to read spectral data and `resample` for interpolation.
        - NaN values in the input data are automatically removed.

    Example:

        >>> import numpy as np
        >>> import rampy as rp
        >>> filenames = np.array(["spectrum1.txt", "spectrum2.txt"])
        >>> x_common = np.linspace(100, 2000, 500)
        >>> spectra_array = rp.spectrarray(filenames, sh=1, sf=1, x=x_common)
    """
    for i in range(len(name)):
        rawspectre = np.genfromtxt(name[i], skip_header=sh, skip_footer=sf)
        rawspectre = rawspectre[~np.isnan(rawspectre).any(1)] # check for nan

        y = resample(rawspectre[:,0],rawspectre[:,1],x) # resample the signal

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

def spectrataux(spectres: np.ndarray) -> np.ndarray:
    """Calculates the rate of change for each frequency in a set of spectra.

    This function fits a second-order polynomial to the y-values at each frequency across 
    multiple spectra and calculates the polynomial coefficients.

    Args:
        spectres (np.ndarray): A 2D array where the first column contains the common x-axis 
            (frequencies) and subsequent columns contain y-values for multiple spectra.

    Returns:
        np.ndarray: A 2D array where the first column contains the frequencies and subsequent 
        columns contain the polynomial coefficients for each frequency.

    Raises:
        ValueError: If curve fitting fails or input data is invalid.

    Example:

        >>> import numpy as np
        >>> import rampy as rp
        >>> spectres = np.array([[100, 10, 12], [200, 20, 22], [300, 30, 32]])
        >>> taux = rp.spectrataux(spectres)
    """
    from scipy.optimize import curve_fit

    # Extract frequencies (x-axis) and initialize output array
    freq = spectres[:, 0]  # First column contains frequencies
    taux = np.zeros((len(freq), 4))  # Output array with frequencies + coefficients
    taux[:, 0] = freq[:]  # First column is frequencies

    # Define second-order polynomial model
    def poly2(x, a, b, c):
        return a * x**2 + b * x + c

    # Fit polynomial coefficients for each frequency
    for i in range(len(freq)):
        y = spectres[i, 1:]  # Extract y-values at this frequency across all spectra
        try:
            popt, _ = curve_fit(poly2, np.arange(len(y)), y, p0=[0.5e-3, 0.5e-4, 1e-6])  # Fit coefficients
            taux[i, 1:] = popt  # Store coefficients in output array
        except RuntimeError as e:
            raise ValueError(f"Curve fitting failed at frequency {freq[i]}: {e}")

    return taux

def spectraoffset(spectre, oft):
    """Applies vertical offsets to spectra.

    This function offsets the y-values of each spectrum by a specified amount.

    Args:
        spectre (np.ndarray): A 2D array where rows represent x-axis values 
            and columns represent spectra (first column is x-axis).
        oft (np.ndarray): A 1D array specifying offset values for each spectrum.

    Returns:
        np.ndarray: A 2D array with the same shape as `spectre`, where specified columns 
        have been vertically offset.

    Raises:
        ValueError: If the length of `oft` does not match the number of spectra.

    Example:

        >>> import numpy as np
        >>> import rampy as rp
        >>> spectre = np.array([[100, 10], [200, 20], [300, 30]])
        >>> offsets = np.array([5])
        >>> offset_spectra = rp.spectraoffset(spectre, offsets)
    """

    out = np.zeros(spectre.shape) # we already say what is the output array
    for i in range(len(oft)): # and we offset the spectra
        out[:,i+1] = spectre[:,i+1] + oft[i]
    return out

#
# Simple data treatment functions
#
def shiftsp(sp: np.ndarray, shift: float) -> np.ndarray:
    """Shifts the x-axis values of a spectrum by a given amount.

    Args:
        sp (np.ndarray): A 2D array where the first column contains x-axis values 
            (e.g., frequency or wavenumber) and subsequent columns contain y-values.
        shift (float): The amount by which to shift the x-axis values. Negative values 
            shift left; positive values shift right.

    Returns:
        np.ndarray: The input array with shifted x-axis values.

    Example:

        >>> import numpy as np
        >>> import rampy as rp
        >>> sp = np.array([[100, 10], [200, 20], [300, 30]])
        >>> shifted_sp = rp.shiftsp(sp, shift=50)
    """
    sp[:,0] = sp[:,0] - shift
    return sp

def flipsp(sp: np.ndarray) -> np.ndarray:
    """Sorts or flips a spectrum along its row dimension based on x-axis values.

    Args:
        sp (np.ndarray): A 2D array where the first column contains x-axis values 
            and subsequent columns contain y-values.

    Returns:
        np.ndarray: The input array sorted in ascending order based on the first column.

    Notes:
        - Uses `np.argsort` to ensure sorting regardless of initial order.

    Example:

        >>> import numpy as np
        >>> import rampy as rp
        >>> sp = np.array([[300, 30], [100, 10], [200, 20]])
        >>> sorted_sp = rp.flipsp(sp)
    """
    return sp[sp[:, 0].argsort()] # we actually use argsort to sort the array in ascending order

def resample(x: np.ndarray, y: np.ndarray, x_new: np.ndarray, fill_value="extrapolate", **kwargs) -> np.ndarray:
    """
    Resamples a y signal along new x-axis values using interpolation.

    Args:
        x (np.ndarray): Original x-axis values.
        y (np.ndarray): Original y-axis values corresponding to `x`.
        x_new (np.ndarray): New x-axis values for resampling.
        fill_value (array-like or (array-like, array_like) or “extrapolate”, optional): behavior of the interpolation for requested points outside of the data range. See [scipy help for details](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html). Default is '“extrapolate”'.
        **kwargs: Additional arguments passed to `scipy.interpolate.interp1d`.

            - kind (str or int): Type of interpolation ('linear', 'cubic', etc.). Default is 'linear'.
            - bounds_error (bool): If True, a ValueError is raised any time interpolation is attempted on a value outside of the range of x (where extrapolation is necessary). If False, out of bounds values are assigned fill_value. By default, an error is raised unless fill_value="extrapolate".

    Returns:
        np.ndarray: Resampled y-values corresponding to `x_new`.

    Example:
    
        >>> import numpy as np
        >>> import rampy as rp
        >>> original_x = np.array([100, 200, 300])
        >>> original_y = np.array([10, 20, 30])
        >>> new_x = np.linspace(100, 300, 5)
        >>> resampled_y = rp.resample(original_x, original_y, new_x)
    """
    f = interp1d(x, y, fill_value=fill_value, **kwargs)
    return f(x_new)

def normalise(y: np.ndarray, x : np.ndarray = None, method: str = "intensity") -> np.ndarray:
    """Normalizes the y signal(s) using specified methods.

    This function normalizes the input y signal(s) based on the chosen method: 
    by area under the curve, maximum intensity, or min-max scaling.

    Args:
        y (np.ndarray): A 2D array of shape (m values, n samples) containing the y values to normalize.
        x (np.ndarray, optional): A 2D array of shape (m values, n samples) containing the x values 
            corresponding to `y`. Required for area normalization. Default is None.
        method (str): The normalization method to use. Options are:
            - 'area': Normalize by the area under the curve.
            - 'intensity': Normalize by the maximum intensity.
            - 'minmax': Normalize using min-max scaling.

    Returns:
        np.ndarray: A 2D array of normalized y signals with the same shape as the input `y`.

    Raises:
        ValueError: If `x` is not provided when using the 'area' normalization method.
        NotImplementedError: If an invalid normalization method is specified.

    Example:
        >>> import numpy as np
        >>> import rampy as rp
        >>> x = np.linspace(0, 10, 100)
        >>> y = rp.gaussian(x, 10., 50., 2.0)
        >>> y_norm = rp.normalise(y, x=x, method="area")
    """

    if method == "area":
        try:
            y = y/simpson(y, x=x, axis=0)
        except:
            raise ValueError("Input array of x values for area normalisation")
    elif method == "intensity":
        y = y/np.max(y, axis=0)

    elif method == "minmax":
        y = (y-np.min(y, axis=0))/(np.max(y, axis=0)-np.min(y, axis=0))
    else:
        raise NotImplementedError("Wrong method name, choose between area, intensity and minmax.")

    return y

def centroid(x: np.ndarray, y: np.ndarray, smoothing: bool = False, **kwargs) -> np.ndarray:
    """Calculates the centroid of y signal(s).

    The centroid is calculated as \( \text{centroid} = \sum(y / \sum(y) \cdot x) \).

    Args:
        x (np.ndarray): A 2D array of shape (m values, n samples) containing the x-axis values.
        y (np.ndarray): A 2D array of shape (m values, n samples) containing the y-axis values.
        smoothing (bool): Whether to smooth the signals before calculating centroids. 
            If True, smoothing is applied using `rampy.smooth` with additional arguments passed via `kwargs`.

    Returns:
        np.ndarray: A 1D array of centroids for each sample in `y`.

    Example:

        >>> import numpy as np
        >>> import rampy as rp
        >>> x = np.linspace(0, 10, 100).reshape(-1, 1)
        >>> y = rp.gaussian(x, 10., 50., 2.0)
        >>> centroids = rp.centroid(x, y)
    """


    y_ = y.copy()

    if smoothing == True:
        for i in range(x.shape[1]):
            y_[:,i] = rampy.smooth(x[:,i],y[:,i],**kwargs)

    return np.sum(y_/np.sum(y_,axis=0)*x,axis=0)

def despiking(x: np.ndarray, y: np.ndarray, neigh:int = 4, threshold:int = 3) -> np.ndarray:
    """Removes spikes from a 1D signal using a threshold-based approach.

    This function identifies spikes in a signal by comparing local residuals to a threshold 
    based on the root mean square error (RMSE). Spikes are replaced with the mean of neighboring points.

    Args:
        x (np.ndarray): A 1D array containing the x-axis values of the signal.
        y (np.ndarray): A 1D array containing the y-axis values of the signal to despike.
        neigh (int): The number of neighboring points to use for calculating average values 
            during despiking and for smoothing. Default is 4.
        threshold (int): The multiplier of RMSE used to identify spikes. Default is 3.

    Returns:
        np.ndarray: A 1D array of the despiked signal.

    Example:

        >>> import numpy as np
        >>> import rampy as rp
        >>> x = np.linspace(0, 10, 100)
        >>> y = rp.gaussian(x, 10., 50., 2.0)
        >>> y_despiked = rp.despiking(x, y)
    """
    # we make sure we work with vectors
    y = y.reshape(-1)
    y_out = y.copy() # So we don’t overwrite y for i in np.arange(len(spikes)):
    
    y_smo = rampy.smooth(x.reshape(-1), y, method="savgol", window_length=neigh)
    rmse_local = np.sqrt((y-y_smo)**2)
    rmse_mean = np.sqrt(np.mean((y-y_smo)**2))

    # if the point is at more than 3 sigmas, we consider it as an outlier
    spikes = rmse_local > threshold* rmse_mean

    for i in range(len(y)):
        if spikes[i] != False: # If we have an spike in position i
            
            # we must be careful avoiding the boundaries
            low_i = i-neigh
            high_i = i+1+neigh
        
            if low_i < 0:
                low_i = 0
            if high_i > len(y):
                high_i = len(y)
            
            w = np.arange(low_i,high_i) # we select 2 m + 1 points around our spike
            w2 = w[spikes[w] == 0] # From such interval, we choose the ones which are not spikes
            y_out[i] = np.mean(y[w2]) # and we average their values
 
    return y_out

def invcm_to_nm(x_inv_cm: np.ndarray, laser_nm: float = 532.0) -> np.ndarray:
    """Converts Raman shifts from inverse centimeters (cm⁻¹) to nanometers (nm).

    Args:
        x_inv_cm (float or np.ndarray): Raman shift(s) in inverse centimeters (cm⁻¹).
        laser_nm (float): Wavelength of the excitation laser in nanometers. Default is 532.0 nm.

    Returns:
        float or np.ndarray: Wavelength(s) in nanometers corresponding to the Raman shifts.

    Example:

        >>> import rampy as rp
        >>> x_inv_cm = 1000.0
        >>> wavelength_nm = rp.invcm_to_nm(x_inv_cm)
    """
    laser_inv_cm = 1.0 / (laser_nm * 1.0e-7)
    x_inv_cm_absolute = laser_inv_cm - x_inv_cm
    return 1.0 / (x_inv_cm_absolute * 1.0e-7)

def nm_to_invcm(x: np.ndarray, laser_nm: float = 532.0) -> np.ndarray:
    """Converts wavelengths from nanometers (nm) to Raman shifts in inverse centimeters (cm⁻¹).

    Args:
        x (float or np.ndarray): Wavelength(s) in nanometers.
        laser_nm (float): Wavelength of the excitation laser in nanometers. Default is 532.0 nm.

    Returns:
        float or np.ndarray: Raman shift(s) in inverse centimeters (cm⁻¹).

    Example:

        >>> import rampy as rp
        >>> wavelength_nm = 600
        >>> shift_inv_cm = nm_to_invcm(wavelength_nm)
    """
    x_inv_cm = 1.0/(x*1.0e-7)
    laser_inv_cm = 1.0/(laser_nm*1.0e-7)
    return laser_inv_cm-x_inv_cm

def convert_x_units(x: np.ndarray, laser_nm: float = 532.0, way: str = "nm_to_cm-1") -> np.ndarray:
    """Converts between nanometers and inverse centimeters for Raman spectroscopy.

    Args:
        x (np.ndarray): Array of x-axis values to convert.
        laser_nm (float): Wavelength of the excitation laser in nanometers. Default is 532.0 nm.
        way (str): Conversion direction. Options are:
            - "nm_to_cm-1": Convert from nanometers to inverse centimeters.
            - "cm-1_to_nm": Convert from inverse centimeters to nanometers.

    Returns:
        np.ndarray: Converted x-axis values.

    Raises:
        ValueError: If an invalid conversion direction is specified.

    Example:
        Convert from nanometers to inverse centimeters:

        >>> import rampy as rp
        >>> x_nm = np.array([600.0])
        >>> x_cm_1 = rp.convert_x_units(x_nm)

        Convert from inverse centimeters to nanometers:

        >>> x_cm_1 = np.array([1000.0])
        >>> x_nm = rp.convert_x_units(x_cm_1, way="cm-1_to_nm")
    """

    if way == "nm_to_cm-1":
        return nm_to_invcm(x, laser_nm = laser_nm)
    elif way == "cm-1_to_nm":
        return invcm_to_nm(x, laser_nm = laser_nm)
    else:
        ValueError("way should be set to either nm_to_cm-1 or cm-1_to_nm")