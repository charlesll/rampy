# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import numpy as np
from scipy.constants import c, h, k

def tlcorrection(x: np.ndarray, y: np.ndarray, temperature: float, wavelength: float, **kwargs) -> tuple:
    """Corrects Raman spectra for temperature and excitation line effects.

    This function applies corrections to Raman spectra to account for temperature and laser 
    excitation wavelength effects. It supports multiple correction equations and normalization 
    methods, making it suitable for a variety of materials and experimental conditions.

    Args:
        x (np.ndarray): Raman shifts in cm⁻¹.
        y (np.ndarray): Intensity values (e.g., counts).
        temperature (float): Temperature in °C.
        wavelength (float): Wavelength of the laser that excited the sample, in nm.
        correction (str, optional): The correction equation to use. Options are:
            - 'long': Default equation from Galeener and Sen (1978) with a \(v_0^3\) coefficient correction.
            - 'galeener': Original equation from Galeener and Sen (1978), based on Shuker and Gammon (1970).
            - 'hehlen': Equation from Hehlen et al. (2010), preserving the Boson peak signal. Default is 'long'.
        normalisation (str, optional): Normalization method for the corrected data. Options are:
            - 'intensity': Normalize by maximum intensity.
            - 'area': Normalize by total area under the curve.
            - 'no': No normalization. Default is 'area'.
        density (float, optional): Density of the studied material in kg/m³, used only with the 'hehlen' equation. 
            Default is 2210.0 (density of silica).

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]:
            - `x` (np.ndarray): Raman shift values after correction.
            - `ycorr` (np.ndarray): Corrected intensity values.
            - `ese_corr` (np.ndarray): Propagated errors calculated as \(\sqrt{y}\) on raw intensities.

    Raises:
        ValueError: If an invalid correction or normalization method is specified.

    Notes:
        - The 'galeener' equation is a modification of Shuker and Gammon's formula to account for 
          \((v_0 - v)^4\) dependence of Raman intensity.
        - The 'long' equation includes a \(v_0^3\) coefficient to remove cubic meter dimensions, 
          as used in several studies like Mysen et al. (1982).
        - The 'hehlen' equation avoids signal suppression below 500 cm⁻¹, preserving features like 
          the Boson peak in glasses.

    References:
        - Galeener, F.L., & Sen, P.N. (1978). Theory of the first-order vibrational spectra of disordered solids. *Physical Review B*, 17(4), 1928–1933.
        - Hehlen, B. (2010). Inter-tetrahedra bond angle of permanently densified silicas extracted from their Raman spectra. *Journal of Physics: Condensed Matter*, 22(2), 025401.
        - Brooker, M.H., Nielsen, O.F., & Praestgaard, E. (1988). Assessment of correction procedures for reduction of Raman spectra. *Journal of Raman Spectroscopy*, 19(2), 71–78.
        - Mysen, B.O., Finger, L.W., Virgo, D., & Seifert, F.A. (1982). Curve-fitting of Raman spectra of silicate glasses. *American Mineralogist*, 67(7-8), 686–695.
        - Neuville, D.R., & Mysen, B.O. (1996). Role of aluminium in the silicate network: In situ high-temperature study of glasses and melts on the join SiO₂-NaAlO₂. *Geochimica et Cosmochimica Acta*, 60(9), 1727–1737.
        - Le Losq, C., Neuville, D.R., Moretti, R., & Roux, J. (2012). Determination of water content in silicate glasses using Raman spectrometry: Implications for the study of explosive volcanism. *American Mineralogist*, 97(5-6), 779–790.
        - Shuker, R., & Gammon, R.W. (1970). Raman-scattering selection-rule breaking and the density of states in amorphous materials. *Physical Review Letters*, 25(4), 222.

    Examples:
        Correct a simple spectrum using default parameters:

        >>> import numpy as np
        >>> x = np.array([100, 200, 300])  # Raman shifts in cm⁻¹
        >>> y = np.array([10, 20, 30])     # Intensity values
        >>> temperature = 25.0             # Temperature in °C
        >>> wavelength = 532.0             # Wavelength in nm
        >>> x_corr, y_corr, ese_corr = correct_spectra(x, y, temperature, wavelength)

        Use a specific correction equation and normalization method:

        >>> x_corr, y_corr, ese_corr = correct_spectra(
                x=x,
                y=y,
                temperature=25,
                wavelength=532,
                correction='hehlen',
                normalisation='intensity',
                density=2500
            )
    """
    # Input validation
    if not isinstance(x, np.ndarray) or not isinstance(y, np.ndarray):
        raise TypeError("x and y must be numpy arrays")
    
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")
        
    if x[-1] < x[0]:  # to raise an error if decreasing x values are provided
        raise ValueError('x values should be increasing.')

    # Getting the kwargs
    correction = kwargs.get('correction', 'long')
    normalisation = kwargs.get('normalisation', 'area')
    density = kwargs.get('density', 2210.0)

    # Physical constants
    #h = 6.626070040e-34    # J S    Planck constant from NIST
    #hb = 1.054571800e-34   # J S    Reduced Planck constant from NIST
    #k = 1.38064852e-23     # J K-1   Boltzmann constant from NIST
    #c = 299792458.         # m s-1   Speed of light from NIST
    
    nu0 = 1.0/wavelength*1e9     # nu0 laser is in m-1 (wavelength is in nm)
    T = temperature + 273.15     # K temperature
    # density is in kg m-3

    # Calculate the relative error on data as sqrt(y). If y <= 0, then error = abs(y).
    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        ese = np.sqrt(np.absolute(y))/np.absolute(y)  # relative errors
    ese[~np.isfinite(ese)] = 1.0  # Set to 1.0 where y is zero

    # Get the Raman shift in m-1
    nu = 100.0*x  # cm-1 -> m-1 Raman shift

    # CORRECTION
    if correction == 'long':
        # Formula used in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012)
        frequency = nu0**3*nu/((nu0-nu)**4)  # frequency correction; dimensionless
        boltzman = 1.0 - np.exp(-h*c*nu/(k*T))  # temperature correction with Boltzmann distribution
        ycorr = y*frequency*boltzman  # correction
    elif correction == 'galeener':
        # Formula from Galeener and Sen (1978) and Brooker et al. (1988)
        frequency = nu/((nu0-nu)**4)  # frequency correction; m^3
        boltzman = 1.0 - np.exp(-h*c*nu/(k*T))  # temperature correction
        ycorr = y*frequency*boltzman  # correction
    elif correction == 'hehlen':
        # Formula from Hehlen et al. 2010
        frequency = 1.0/(nu0**3*density)  # frequency + density correction; m/kg
        boltzman = 1.0 - np.exp(-h*c*nu/(k*T))  # dimensionless
        ycorr = nu*y*frequency*boltzman  # correction
    else:
        raise ValueError(f"Unknown correction method: {correction}. Choose 'long', 'galeener', or 'hehlen'.")

    # Normalization
    if normalisation == 'area':
        area = np.trapz(ycorr, x)
        if area != 0:
            ycorr = ycorr/area  # area normalization
        else:
            raise ValueError("Cannot normalize by area: integral is zero")
    elif normalisation == 'intensity':
        max_intensity = np.max(ycorr)
        if max_intensity != 0:
            ycorr = ycorr/max_intensity  # max. intensity normalization
        else:
            raise ValueError("Cannot normalize by intensity: maximum is zero")
    elif normalisation == 'no':
        pass  # No normalization
    else:
        raise ValueError(f"Unknown normalization method: {normalisation}. Choose 'area', 'intensity', or 'no'.")

    # Calculate error on corrected data
    ese_corr = ese * ycorr

    return x, ycorr, ese_corr
