# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import numpy as np
from scipy.integrate import simpson
import rampy

def peakarea(shape: str, amp: float, HWHM: float, pos: float = None, 
             a3: float = None, L_ratio: float = None, ese_amp: float = None, 
             ese_HWHM: float = None) -> tuple[float, float | None]:
    """
    Computes the area of a peak given its shape and parameters.

    For Gaussian peaks, the area is calculated analytically. For other shapes 
    (Lorentzian, pseudo-Voigt, Pearson VII), the area is calculated using numerical 
    integration (trapezoidal rule).

    Args:
        shape (str): The shape of the peak. Must be one of 'gaussian', 'lorentzian', 
            'pseudovoigt', or 'pearson7'.
        amp (float, list of floats, or np.ndarray): Amplitude of the peak.
        HWHM (float, list of floats, or np.ndarray): Half-width at half-maximum of the peak.
        pos (float, list of floats, np.ndarray, optional): Peak position (required for non-Gaussian shapes). Default is None.
        a3 (float, list of floats, np.ndarray, optional): Shape parameter for Pearson VII peaks. Required if `shape='pearson7'`.
        L_ratio (float, list of floats, np.ndarray, optional): Lorentzian component ratio for pseudo-Voigt peaks. Must be between 0 and 1.
            Required if `shape='pseudovoigt'`. Default is None.
        ese_amp (float, list of floats, np.ndarray, optional): Standard deviation of the amplitude. Used to calculate error on area.
            Default is None.
        ese_HWHM (float, list of floats, np.ndarray, optional): Standard deviation of HWHM. Used to calculate error on area.
            Default is None.

    Returns:
        tuple[float, float | None]: A tuple containing:
            - area (float): The computed area of the peak.
            - ese_area (float or None): The error on the computed area if `ese_amp` and 
              `ese_HWHM` are provided; otherwise, None.

    Raises:
        ValueError: If required parameters are missing or invalid, or if the list of floats or the np.ndarray for the peak parameters have different sizes.
        NotImplementedError: If an unsupported shape is specified.

    Notes:
        - For Gaussian peaks, the formula used is:
          \( \text{area} = \sqrt{\pi / \ln(2)} \cdot \text{amp} \cdot \text{HWHM} \).
        - For other shapes, numerical integration is performed over a range of \( \pm 10 \cdot \text{HWHM} \).
    
    Example:
        Compute the area of a Gaussian peak:

        >>> amp = 10
        >>> HWHM = 2
        >>> area, ese_area = peakarea("gaussian", amp=amp, HWHM=HWHM)
        >>> print(area)

        Compute the area of a Lorentzian peak with numerical integration:

        >>> amp = 10
        >>> HWHM = 2
        >>> pos = 5
        >>> area, _ = peakarea("lorentzian", amp=amp, HWHM=HWHM, pos=pos)
        >>> print(area)

        Compute the area of several Lorentzian peaks with numerical integration (beware that in this case, the lists should have equal sizes):

        >>> amp = [10., 5., 3.]
        >>> HWHM = 2 # they share a common width
        >>> pos = [2., 5., 5.5]
        >>> area, _ = peakarea("lorentzian", amp=amp, HWHM=HWHM, pos=pos)
        >>> print(area)


    """
    # Validate input parameters
    if shape not in ["gaussian", "lorentzian", "pseudovoigt", "pearson7"]:
        raise NotImplementedError("Supported shapes are 'gaussian', 'lorentzian', 'pseudovoigt', and 'pearson7'.")
    
    # Ensure inputs are arrays for consistent handling
    amp = np.asarray(amp)
    HWHM = np.asarray(HWHM)
    pos = np.asarray(pos) if pos is not None else None
    a3 = np.asarray(a3) if a3 is not None else None
    L_ratio = np.asarray(L_ratio) if L_ratio is not None else None
    ese_amp = np.asarray(ese_amp) if ese_amp is not None else None
    ese_HWHM = np.asarray(ese_HWHM) if ese_HWHM is not None else None
    
    # Gaussian peak: Analytical calculation
    if shape == "gaussian":
        area = np.sqrt(np.pi / np.log(2)) * amp * HWHM
        ese_area = None
        
        # Calculate error on area if uncertainties are provided
        if ese_amp is not None and ese_HWHM is not None:
            ese_area = np.sqrt(
                ((np.sqrt(np.pi / np.log(2)) * HWHM) ** 2) * ese_amp**2 +
                ((np.sqrt(np.pi / np.log(2)) * amp) ** 2) * ese_HWHM**2
            )
        
        return area, ese_area

    # Non-Gaussian peaks: Numerical integration
    if pos is None:
        raise ValueError("Parameter 'pos' must be provided for non-Gaussian peaks.")

    # we create a broad X axis with the extreme values to cover all the peaks
    x_int = np.linspace(np.min(pos) - 10 * np.max(HWHM), np.max(pos) + 10 * np.max(HWHM), 10000).reshape(-1,1)

    if shape == "lorentzian":
        y_int = rampy.lorentzian(x_int, amp=amp, freq=pos, HWHM=HWHM)
    elif shape == "pseudovoigt":
        if L_ratio is None or not (0 <= L_ratio <= 1):
            raise ValueError("Parameter 'L_ratio' must be provided and between 0 and 1 for pseudo-Voigt peaks.")
        y_int = rampy.pseudovoigt(x_int, amp=amp, freq=pos, HWHM=HWHM, L_ratio=L_ratio)
    elif shape == "pearson7":
        if a3 is None:
            raise ValueError("Parameter 'a3' must be provided for Pearson VII peaks.")
        y_int = rampy.pearson7(x_int, a0=amp, a1=pos, a2=HWHM, a3=a3)

    # Calculate area using numerical integration
    area = simpson(y_int, x=x_int, axis=0)
    
    return area, None

def gaussianarea(amp,HWHM,**options):
    """returns the area of a Gaussian peak

    Args:
        amp (float or ndarray): amplitude of the peak
        HWHM (float or ndarray): half-width at half-maximum
        eseAmplitude (float or ndarray, optional): standard deviation on amp; Default = None
        eseHWHM (float or ndarray, optional): standard deviation on HWHM; Default = None

    Returns:
        area (float or ndarray): peak area
        esearea  (float or ndarray): error on peak area; will be None if no errors on amp and HWHM were provided.
    """
    area = np.sqrt(np.pi/np.log(2))*amp*HWHM
    if options.get("eseAmplitude") != None:
        eseAmplitude = options.get("eseAmplitude")
        if options.get("eseHWHM") != None:
            eseHWHM = options.get("eseHWHM")
            esearea = np.sqrt((np.pi/np.log(2)*HWHM)**2 * eseAmplitude**2 + (np.pi/np.log(2)*amp)**2 * eseHWHM**2)
    else:
        esearea = None

    return area, esearea
