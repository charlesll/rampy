# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import numpy as np
import scipy
from scipy import signal
from scipy.interpolate import make_smoothing_spline
import plotly.graph_objects as go

def smooth(x: np.ndarray, y: np.ndarray, method: str = "whittaker", **kwargs) -> np.ndarray:
    """Smooths the provided signal using various smoothing algorithms.

    This function applies different smoothing techniques to the input signal, 
    including Whittaker smoothing, Savitzky-Golay filtering, spline-based methods, 
    and window-based filters. Each method is designed to reduce noise while preserving 
    signal features.

    Args:
        x (np.ndarray): The x-axis values (e.g., time or wavelength). For window-based methods, 
          these values should be equally spaced.
        y (np.ndarray): The y-axis values to be smoothed.
          method (str): The smoothing method to apply. Default is "whittaker".
            Available options:
                - "whittaker": Whittaker smoother (Eilers 2003).
                - "savgol": Savitzky-Golay filter.
                - "GCVSmoothedNSpline": Spline with generalized cross-validation.
                - "MSESmoothedNSpline": Spline with mean squared error criterion 
                  (requires the gcvspline library).
                - "DOFSmoothedNSpline": Spline with degrees of freedom criterion 
                  (requires the gcvspline library).
                - "flat": Moving average window.
                - "hanning", "hamming", "bartlett", "blackman": Various window filters.
        **kwargs: Additional parameters specific to the chosen method.
            - window_length (int): Length of the filter window for window-based methods 
              and Savitzky-Golay filter. Must be a positive odd integer. Default is 5.
            - polyorder (int): Polynomial order for Savitzky-Golay filter. Must be less 
              than `window_length`. Default is 2.
            - Lambda (float): Smoothing parameter for the Whittaker filter. Higher values 
              produce smoother fits. Default is \(10^5\).
            - d (int): Difference order in the Whittaker filter. Default is 2.
            - ese_y (float or np.ndarray): Errors associated with y values for spline methods. 
              Default is 1.0.

    Returns:
        np.ndarray: The smoothed signal sampled on `x`.

    Raises:
        ValueError: If the input vector is smaller than the window size or if an invalid 
            smoothing method is specified.
        ImportError: If `gcvspline` is not installed when using MSESmoothedNSpline or DOFSmoothedNSpline.

    Notes:
        - The "whittaker" method implements the perfect smoother described by Eilers (2003).
        - The "savgol" method uses `scipy.signal.savgol_filter()`.
        - Spline methods require `scipy >= 1.10.0` or the `gcvspline` package.
        - Window-based methods are implemented following the SciPy Cookbook.

    Examples:
        Smooth a noisy signal using Savitzky-Golay filtering:

        >>> import numpy as np
        >>> x = np.linspace(0, 100, 101)
        >>> y = np.sin(x / 10) + np.random.normal(0, 0.1, size=x.size)
        >>> smoothed_signal = smooth(x, y, method="savgol", window_length=7, polyorder=3)

        Apply Whittaker smoothing:

        >>> smoothed_signal = smooth(x, y, method="whittaker", Lambda=1e6)

        Use a moving average window:

        >>> smoothed_signal = smooth(x, y, method="flat", window_length=5)
    """


    window_len = kwargs.get("window_length",5)
    polyorder = kwargs.get("polyorder",2)
    lam = kwargs.get("Lambda",10.0**5)
    d = kwargs.get("d",2)
    ese_y = kwargs.get("ese_y",1.0)

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if not method in ['GCVSmoothedNSpline','MSESmoothedNSpline','DOFSmoothedNSpline','whittaker','savgol',
                      'flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Method should be on 'GCVSmoothedNSpline','MSESmoothedNSpline','DOFSmoothedNSpline','whittaker','savgol','flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    if method == "GCVSmoothedNSpline":
            w = 1.0 / (np.ones((y.shape[0],1)) * ese_y) # errors
            flt = make_smoothing_spline(x.reshape(-1),y.reshape(-1),w.reshape(-1))
            return flt(x)
    
    if (method == "MSESmoothedNSpline") or (method == "DOFSmoothedNSpline"): # gcvspline methods

        try: # we test if gcvspline is installed
            import gcvspline
        except ImportError:
            print('ERROR: Install gcvspline to use this smoothing mode (needs a working FORTRAN compiler).')

        w = 1.0 / (np.ones((y.shape[0],1)) * ese_y) # errors

        #if method == "GCVSmoothedNSpline":
            #flt = gcvspline.GCVSmoothedNSpline(x.reshape(-1),y.reshape(-1),w.reshape(-1))
        if method == "MSESmoothedNSpline":
            flt = gcvspline.MSESmoothedNSpline(x.reshape(-1),y.reshape(-1),w.reshape(-1))
        elif method == "DOFSmoothedNSpline":
            flt = gcvspline.DOFSmoothedNSpline(x.reshape(-1),y.reshape(-1),w.reshape(-1))
        return flt(x)

    elif method == "whittaker": # whittaker smoother
        return whittaker(y,Lambda=lam,d=d)
    
    elif method == "savgol": # Savtisky-Golay filter
        return scipy.signal.savgol_filter(y, window_len, polyorder)
    
    elif method in frozenset(('flat', 'hanning', 'hamming', 'bartlett', 'blackman')): # various window filters, from https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html?highlight=smooth

        s=np.r_[y[window_len-1:0:-1],y,y[-2:-window_len-1:-1]]
        if method == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:
            w=getattr(np,method)(window_len)

        y_filt = np.convolve(w/w.sum(),s,mode='valid')
        shift = int((len(y_filt) - len(y))/2)

        return y_filt[shift:-shift]

def whittaker(y: np.ndarray, weights: np.ndarray = None, **kwargs) -> np.ndarray:
    """Smooths a signal using the Whittaker smoother.

    This function applies Whittaker smoothing to reduce noise in a signal while preserving 
    its features. It uses penalized least squares optimization with a specified smoothing 
    coefficient.

    Args:
        y (np.ndarray): An array containing the values to smooth. The data should be equally spaced.
        weights (np.ndarray, optional): Weights for the signal. Regions with weight 0 will not 
        be smoothed. Default is uniform weights across all points.
        **kwargs: Additional parameters for Whittaker smoothing.
            - Lambda (float): The smoothing coefficient; higher values produce smoother results. 
              Default is \(10^5\).

    Returns:
        np.ndarray: An array containing the smoothed values.

    Raises:
        ValueError: If `y` or `weights` are invalid.

    Notes:
        - This implementation follows the description provided by Eilers (2003).

    References:
        Eilers, P.H.C., 2003. A Perfect Smoother. Anal. Chem. 75, 3631–3636.

    Examples:
        Smooth a noisy signal using default weights and Lambda:

        >>> import numpy as np
        >>> y = np.sin(np.linspace(0, 10, 100)) + np.random.normal(0, 0.1, size=100)
        >>> smoothed_signal = whittaker(y)

        Apply custom weights and Lambda:

        >>> weights = np.ones(100)
        >>> weights[40:60] = 0  # Do not smooth this region
        >>> smoothed_signal = whittaker(y, weights=weights, Lambda=1e6)
    """

    lam = kwargs.get('Lambda', 1.0 * 10**5)
    L = len(y)
    D = scipy.sparse.csc_matrix(np.diff(np.eye(L), 2))

    if weights is None:
        weights = np.ones(L)

    W = scipy.sparse.spdiags(weights, 0, L, L)
    Z = W + lam * D.dot(D.transpose())
    z = scipy.sparse.linalg.spsolve(Z, weights * y)

    return z 

def spectrafilter(spectre : np.ndarray, filtertype: str, fq: int, numtaps: int, columns: np.ndarray) -> np.ndarray:
    """Applies a Butterworth filter to specific columns of spectral data.

    This function filters specific frequencies in the provided spectral data using a Butterworth filter. 
    It supports low-pass, high-pass, band-pass, and band-stop filtering. The input spectra must be 
    provided as an array where the first column represents the x-axis (e.g., wavelength or time), 
    and subsequent columns represent the y-axis values of multiple spectra.

    Args:
        spectre (np.ndarray): A 2D array of spectral data. The first column contains x-axis values 
            (e.g., time or wavelength), and subsequent columns contain y-axis values for multiple spectra.
        filtertype (str): The type of filter to apply. Choose from:
            - 'low': Low-pass filter.
            - 'high': High-pass filter.
            - 'bandpass': Band-pass filter.
            - 'bandstop': Band-stop filter.
        fq (float or np.ndarray): The cutoff frequency/frequencies for the filter. For 'low' or 'high' filters, 
            this is a single float value. For 'bandpass' or 'bandstop', this must be a sequence of two values 
            `[low_cutoff, high_cutoff]`.
        numtaps (int): The order of the Butterworth filter. Higher values result in sharper transitions but may 
            introduce more ringing artifacts.
        columns (np.ndarray): An array specifying which columns in the input `spectre` to apply the filter to. 
            Each value should correspond to a valid column index.

    Returns:
        np.ndarray: A 2D array with the same shape as `spectre`, where the specified columns have been filtered.

    Raises:
        ValueError: If an invalid `filtertype` is provided or if `fq` is improperly specified for 'bandpass' or 
            'bandstop' filters.

    Notes:
        - The x-axis values (first column) are copied directly to the output without modification.
        - The cutoff frequencies (`fq`) are normalized based on the Nyquist frequency, which is calculated from 
          the sampling rate inferred from the x-axis spacing.
        - The Butterworth filter is applied using zero-phase filtering (`scipy.signal.filtfilt`) to avoid phase distortion.

    Examples:
        Apply a low-pass filter to a single spectrum:

        >>> import numpy as np
        >>> import rampy
        >>> x = np.linspace(0, 10, 1000)  # Time axis
        >>> y = np.sin(2 * np.pi * 5 * x) + np.sin(2 * np.pi * 50 * x)  # Signal with two frequencies
        >>> spectre = np.column_stack((x, y))
        >>> filtered_spectre = rampy.spectrafilter(spectre, filtertype='low', fq=10, numtaps=4, columns=np.array([1]))

        Apply a band-pass filter to multiple spectra:

        >>> spectre = np.column_stack((x, y, y * 0.5))  # Two spectra
        >>> filtered_spectre = spectrafilter(spectre, filtertype='bandpass', fq=[5, 15], numtaps=4, columns=np.array([1, 2]))

    References:
        - Butterworth Filter Design: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
        - Zero-Phase Filtering: https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
    """
    # Input validation
    if filtertype not in ['low', 'high', 'bandstop', 'bandpass']:
        raise ValueError(f"Invalid filtertype: {filtertype}. Choose from 'low', 'high', 'bandstop', 'bandpass'")
    
    if (filtertype in ['bandstop', 'bandpass']) and (not hasattr(fq, '__len__') or len(fq) != 2):
        raise ValueError(f"{filtertype} filter requires fq to be a sequence of two values [low, high]")
    
    # Create output array and copy x-axis
    out = np.zeros_like(spectre)
    out[:, 0] = spectre[:, 0]
    
    # Calculate sampling parameters
    sample_spacing = spectre[1, 0] - spectre[0, 0]
    samplerate = 1 / sample_spacing  # Hertz
    nyq_rate = samplerate / 2  # Nyquist frequency
    
    # Normalize cutoff frequency(ies)
    if filtertype in ['low', 'high']:
        normalized_cutoff = fq / nyq_rate
    else:  # bandpass or bandstop
        normalized_cutoff = [fq[0] / nyq_rate, fq[1] / nyq_rate]
    
    # Design filter once
    b, a = signal.butter(numtaps, normalized_cutoff, btype=filtertype)
    
    # Apply filter to each specified column
    for col in columns:
        out[:, col] = signal.filtfilt(b, a, spectre[:, col])
    
    return out