# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline, make_smoothing_spline
from scipy.spatial import ConvexHull
import scipy.sparse as sparse
from numpy.linalg import norm
from sklearn import preprocessing
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
import warnings
import rampy

def get_portion_interest(x: np.ndarray, y: np.ndarray, roi) -> np.ndarray:
    print('Warning: get_portion_interest will be renamed extract_signal in a next version. You can already use the new name!')
    return extract_signal(x, y, roi)

def extract_signal(x: np.ndarray, y: np.ndarray, roi) -> np.ndarray:
    """Extracts the signal from specified regions of interest (ROI) in the x-y data.

    This function selects and extracts portions of the input x-y data based on the 
    specified regions of interest (ROI) provided in `roi`. Each region is defined 
    by a lower and upper bound along the x-axis.

    Args:
        x (ndarray): The x-axis values (e.g., time, wavelength, or other independent variables).
        y (ndarray): The y-axis values corresponding to the x-axis values (e.g., signal intensity).
        roi (ndarray or list of lists): Regions of interest (ROI) where the signal should be extracted.
            Must be an n x 2 array or a list of lists, where `n` is the number of regions to extract.
            Each sublist or row should contain two elements:
              - The lower bound of the region (inclusive).
              - The upper bound of the region (inclusive).

            Example:
                - Array: `np.array([[10, 20], [50, 70]])`
                - List: `[[10, 20], [50, 70]]`

    Returns:
        ndarray: A 2-column array containing the extracted x-y signals from the specified regions.
            The first column contains the x values, and the second column contains the corresponding y values.

    Raises:
        ValueError: If `roi` is not a valid n x 2 array or list of lists, or if any region in `roi` 
            falls outside the range of `x`.

    Notes:
        - Overlapping regions in `roi` are not merged; they are processed as separate regions.
        - If no valid regions are found within `roi`, an empty array is returned.

    Examples:
        Extracting signal from two regions in an x-y dataset:

        >>> import numpy as np
        >>> x = np.linspace(0, 100, 101)
        >>> y = np.sin(x / 10) + np.random.normal(0, 0.1, size=x.size)
        >>> roi = [[10, 20], [50, 70]]
        >>> extracted_signal = extract_signal(x, y, roi)
        >>> print(extracted_signal)
    """
    
    # Check if roi is a numpy array or a list of lists, and convert
    if isinstance(roi, list):
        roi = np.array(roi)
    
    # check the size
    if roi.ndim != 2 or roi.shape[1] != 2:
        raise ValueError("roi must be an n x 2 array.")

    xy = np.column_stack((x, y))
    extracted_regions = [xy[(xy[:, 0] > bounds[0]) & (xy[:, 0] < bounds[1])] for bounds in roi]
    return np.vstack(extracted_regions)

def validate_input(x, y):
    """Validate the input arrays."""
    if not isinstance(x, np.ndarray) or not isinstance(y, np.ndarray):
        raise TypeError("x and y must be numpy arrays.")

    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape.")

def validate_roi(roi):
    # Check if roi is a numpy array or a list of lists
    if isinstance(roi, np.ndarray):
        if roi.ndim != 2 or roi.shape[1] != 2:
            raise ValueError("roi must be an n x 2 array.")
    elif isinstance(roi, list):
        if not all(isinstance(region, list) and len(region) == 2 for region in roi):
            raise ValueError("roi must be a list of lists, where each sublist contains two elements.")
    else:
        raise TypeError("roi must be either a numpy array or a list of lists.")

def standardize_data(x: np.ndarray, y: np.ndarray):
    """Standardize the data using sklearn's StandardScaler."""
    X_scaler = preprocessing.StandardScaler().fit(x.reshape(-1, 1))
    Y_scaler = preprocessing.StandardScaler().fit(y.reshape(-1, 1))
    x_scaled = X_scaler.transform(x.reshape(-1, 1)).flatten()
    y_scaled = Y_scaler.transform(y.reshape(-1, 1)).flatten()
    return x_scaled, y_scaled, X_scaler, Y_scaler

def baseline(x_input: np.ndarray, 
             y_input: np.ndarray, 
             roi = None,
             method: str = "poly", 
             **kwargs) -> tuple:
    """Subtracts a baseline from an x-y spectrum using various methods.

    This function performs baseline subtraction on spectroscopic data by fitting a model 
    to the background signal. It supports multiple correction methods, including polynomial 
    fitting, splines, and advanced penalized least squares techniques.

    Args:
        x_input (ndarray): The x-axis values (e.g., wavelength or wavenumber).
        y_input (ndarray): The y-axis values (e.g., intensity or absorbance).
        roi (ndarray): Regions of interest for baseline fitting. Default is the entire range of `x_input`.
        method (str, optional): The method used for baseline fitting. Default is "poly".
            Available options:
                - "poly": Polynomial fitting. Requires `polynomial_order`.
                - "unispline": Spline fitting using Scipy's `UnivariateSpline`. Requires `s`.
                - "gcvspline": Spline fitting using GCVSpline. Requires `s` and optionally `ese_y`.
                - "gaussian": Gaussian function fitting. Requires `p0_gaussian`.
                - "exp": Exponential background fitting. Requires `p0_exp`.
                - "log": Logarithmic background fitting. Requires `p0_log`.
                - "rubberband": Rubberband correction using convex hull interpolation.
                - "als": Asymmetric Least Squares correction. Requires `lam`, `p`, and `niter`.
                - "arPLS": Asymmetrically Reweighted Penalized Least Squares smoothing. Requires `lam` and `ratio`.
                - "drPLS": Doubly Reweighted Penalized Least Squares smoothing. Requires `niter`, `lam`, `eta`, and `ratio`.
                - "whittaker": Whittaker smoothing with weights applied to regions of interest. Requires `lam`.
                - "GP": Gaussian process method using a rational quadratic kernel.

        **kwargs: Additional parameters specific to the chosen method.

            - polynomial_order (int): Degree of the polynomial for the "poly" method. Default is 1.
            - s (float): Spline smoothing coefficient for "unispline" and "gcvspline". Default is 2.0.
            - ese_y (ndarray): Errors associated with y_input for the "gcvspline" method.
              Defaults to sqrt(abs(y_input)) if not provided.
            - lam (float): Smoothness parameter for methods like "als", "arPLS", and others. Default is 1e5.
            - p (float): Weighting parameter for ALS method. Recommended values are between 0.001 and 0.1. Default is 0.01.
            - ratio (float): Convergence ratio parameter for arPLS/drPLS methods. Default is 0.01 for arPLS and 0.001 for drPLS.
            - niter (int): Number of iterations for ALS/drPLS methods. Defaults are 10 for ALS and 100 for drPLS.
            - eta (float): Roughness parameter for drPLS method, between 0 and 1. Default is 0.5.
            - p0_gaussian (list): Initial parameters [a, b, c] for Gaussian fitting: 
              \(y = a \cdot \exp(-\log(2) \cdot ((x-b)/c)^2)\). Default is [1., 1., 1.].
            - p0_exp (list): Initial parameters [a, b, c] for exponential fitting: 
              \(y = a \cdot \exp(b \cdot (x-c))\). Default is [1., 1., 0.].
            - p0_log (list): Initial parameters [a, b, c, d] for logarithmic fitting: 
              \(y = a \cdot \log(-b \cdot (x-c)) - d \cdot x^2\). Default is [1., 1., 1., 1.].

    Returns:
        tuple:
            corrected_signal (ndarray): The signal after baseline subtraction.
            baseline_fitted (ndarray): The fitted baseline.

    Raises:
        ValueError: If the specified method is not recognized or invalid parameters are provided.

    Notes:
        The input data is standardized before fitting to improve numerical stability during optimization.
        The fitted baseline is transformed back to the original scale before subtraction.

    Examples:
        Example with polynomial baseline correction:

        >>> import numpy as np
        >>> x = np.linspace(50, 500, nb_points)
        >>> noise = 2.0 * np.random.normal(size=nb_points)
        >>> background = 5 * np.sin(x / 50) + 0.1 * x
        >>> peak = rampy.gaussian(x, 100.0, 250.0, 7.0)
        >>> y = peak + background + noise
        >>> roi = np.array([[0, 200], [300, 500]])
        >>> corrected_signal, baseline = baseline(x, y, method="poly", roi = roi, polynomial_order=5)

        Example with GCVSpline algorithm:

        >>> corrected_signal, baseline = baseline(x, y, method="gcvspline", roi=roi, s=2.0)

        Example with Whittaker smoothing:

        >>> corrected_signal, baseline = baseline(x, y, method="whittaker", roi=roi)

        Example with Gaussian process:

        >>> corrected_signal, baseline = baseline(x, y, method="GP", roi=roi)

        Example with rubberband correction:

        >>> corrected_signal, baseline = baseline(x, y, method="rubberband")

        Example with ALS algorithm:

        >>> corrected_signal, baseline = baseline(x, y, method="als", lam=1e5, p=0.01)
    """

    # first check the input
    validate_input(x_input, y_input)

    # then get the roi or set it to the entire x_input range
    if roi is not None:
        validate_roi(roi)
    else:
        roi = [[np.min(x_input), np.max(x_input)]]

    # now extract the signal and standardize the data
    yafit_unscaled = extract_signal(x_input, y_input, roi)
    x, y, X_scaler, Y_scaler = standardize_data(x_input, y_input)

    yafit = np.copy(yafit_unscaled)
    yafit[:, 0] = X_scaler.transform(yafit_unscaled[:, 0].reshape(-1, 1)).flatten()
    yafit[:, 1] = Y_scaler.transform(yafit_unscaled[:, 1].reshape(-1, 1)).flatten()
    roi_scaled = X_scaler.transform(np.array(roi).reshape(-1,1)).reshape(-1, 2)

    baseline_fitted = fit_baseline(x, y, roi_scaled, yafit, method, **kwargs)

    corrected_signal = y_input.reshape(-1, 1) - Y_scaler.inverse_transform(baseline_fitted.reshape(-1, 1))
    return corrected_signal, Y_scaler.inverse_transform(baseline_fitted.reshape(-1, 1))

def fit_baseline(x: np.ndarray, y: np.ndarray, roi_: np.ndarray, yafit: np.ndarray, method: str, **kwargs) -> np.ndarray:
    """Fit the baseline using the specified method."""
    if method == 'poly':
        poly_order = kwargs.get('polynomial_order', 1)
        coeffs = np.polyfit(yafit[:, 0], yafit[:, 1], poly_order)
        return np.polyval(coeffs, x)

    elif method == 'unispline':
        spline = UnivariateSpline(yafit[:, 0], 
                                  yafit[:, 1], 
                                  s=kwargs.get('s', 2.0))
        return spline(x)

    elif method == 'gcvspline':
        splinesmooth = kwargs.get('s', None)
        if splinesmooth is not None:
            splinesmooth /= 10
        ese = kwargs.get('ese_y', np.sqrt(np.abs(yafit[:, 1])))
        spline = make_smoothing_spline(yafit[:, 0], yafit[:, 1], w=ese, lam=splinesmooth)
        return spline(x)

    elif method == 'gaussian':
        p0_gauss = kwargs.get('p0_gaussian', [1., 1., 1.])
        coeffs, _ = curve_fit(rampy.gaussian, yafit[:, 0], yafit[:, 1], p0=p0_gauss)
        return rampy.gaussian(x, *coeffs)

    elif method == 'exp':
        p0_exp = kwargs.get('p0_exp', [1., 1., 0.])
        coeffs, _ = curve_fit(rampy.funexp, yafit[:, 0], yafit[:, 1], p0=p0_exp)
        return rampy.funexp(x, *coeffs)

    elif method == 'log':
        p0_log = kwargs.get('p0_log', [1., 1., 1., 1.])
        coeffs, _ = curve_fit(rampy.funlog, yafit[:, 0], yafit[:, 1], p0=p0_log)
        return rampy.funlog(x, *coeffs)

    elif method == 'rubberband':
        return rubberband_baseline(x, y)

    elif method == 'als':
        return als_baseline(y, 
                            kwargs.get('lam', 1.0e5), 
                            kwargs.get('p', 0.01), 
                            kwargs.get('niter', 10))

    elif method == 'arPLS':
        return arPLS_baseline(y, 
                              kwargs.get('lam', 1.0e5), 
                              kwargs.get('ratio', 0.01))

    elif method == 'drPLS':
        return drPLS_baseline(y, 
                              kwargs.get('niter', 100), 
                              kwargs.get('lam', 1.0e5), 
                              kwargs.get('eta', 0.5), 
                              kwargs.get('ratio', 0.001))
    
    elif method == 'whittaker':
        # Create weights array with 0s outside the regions of interest
        weights = np.zeros_like(x)
        for bounds in roi_:
            weights[(x >= bounds[0]) & (x <= bounds[1])] = 1

        return rampy.whittaker(y, 
                               weights=weights, 
                               Lambda=kwargs.get('lam', 1.0e5))

    elif method == "GP":
        # Gaussian process method with the rational quadratic kernel
        gp_scale = 1.0
        kernel = np.std(y) * RBF(length_scale=gp_scale*(np.max(yafit[:,0])-np.min(yafit[:,0])), 
                                length_scale_bounds=[0.001,100000])

        # declare a GP, with a noise alpha estimated from yafit
        gaussian_process = GaussianProcessRegressor(kernel=kernel, 
                                                    alpha=np.std(yafit[:,1]), # noise in data a priori
                                                    n_restarts_optimizer=1)

        # perform the fit
        gaussian_process.fit(yafit[:,0].reshape(-1,1), yafit[:,1])

        # get result
        return gaussian_process.predict(x.reshape(-1,1))

    else:
        raise ValueError("Method not found, check you entered the right name.")

def rubberband_baseline(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Perform rubberband baseline correction.

    Parameters
    ----------
    x : ndarray
        The x-axis values.
    y : ndarray
        The y values corresponding to x.

    Returns
    -------
    baseline_fitted : ndarray
        The fitted baseline.

    Raises
    ------
    ValueError
        If the convex hull does not have enough points for interpolation.
    """
    # Ensure x and y are numpy arrays
    x = np.asarray(x)
    y = np.asarray(y)

    # Find the convex hull
    points = np.column_stack((x, y))
    hull = ConvexHull(points)
    vertices = hull.vertices

    # Roll the vertices to start from the lowest x-value
    vertices = np.roll(vertices, -vertices.argmin())

    # Ensure there are enough points for interpolation
    if len(vertices) < 2:
        raise ValueError("Convex hull does not have enough points for interpolation.")

    # Leave only the ascending part
    vertices = vertices[:vertices.argmax() + 1]

    # Perform linear interpolation
    baseline_fitted = np.interp(x, x[vertices], y[vertices], left=y[vertices[0]], right=y[vertices[-1]])

    return baseline_fitted

def als_baseline(y: np.ndarray, lam: float, p: float, niter: int) -> np.ndarray:
    """Asymmetric Least Squares baseline fitting."""
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sparse.linalg.spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z

def arPLS_baseline(y: np.ndarray, lam: float, ratio: float) -> np.ndarray:
    """Asymmetrically Reweighted Penalized Least Squares baseline fitting."""
    N = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(N), 2))
    H = lam * D.dot(D.transpose())
    w = np.ones(N)
    while True:
        W = sparse.spdiags(w, 0, N, N)
        Z = W + H
        z = sparse.linalg.spsolve(Z, w * y)
        d = y - z
        dn = d[d < 0]
        m = np.mean(dn)
        s = np.std(dn)
        wt = 1.0 / (1 + np.exp(2 * (d - (2 * s - m)) / s))
        if norm(w - wt) / norm(w) < ratio:
            break
        w = wt
    return z

def drPLS_baseline(y: np.ndarray, niter: int, lam: float, eta: float, ratio: float) -> np.ndarray:
    """Doubly Reweighted Penalized Least Squares baseline fitting."""
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2), format='csr').dot(sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2), format='csr').transpose())
    D_1 = sparse.diags([-1, 1], [0, -1], shape=(L, L - 1), format='csr').dot(sparse.diags([-1, 1], [0, -1], shape=(L, L - 1), format='csr').transpose())
    w = np.ones(L)
    W = sparse.diags(w, format='csr')
    Z = w
    for _ in range(niter):
        W.setdiag(w)
        Z_prev = Z
        Z = sparse.linalg.spsolve(W + D_1 + lam * (sparse.diags(np.ones(L), format='csr') - eta * W) * D, W * y, permc_spec='NATURAL')
        if np.linalg.norm(Z - Z_prev) > ratio:
            d = y - Z
            d_negative = d[d < 0]
            sigma_negative = np.std(d_negative)
            mean_negative = np.mean(d_negative)
            w = 0.5 * (1 - np.exp(_) * (d - (-mean_negative + 2 * sigma_negative)) / sigma_negative / (1 + np.abs(np.exp(_) * (d - (-mean_negative + 2 * sigma_negative)) / sigma_negative)))
        else:
            break
    return Z