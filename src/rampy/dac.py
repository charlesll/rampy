import numpy as np

def ruby_T_corr(T_K: float | np.ndarray) -> float | np.ndarray:
    """
    Computes the temperature correction for the Ruby pressure scale.

    This function calculates the temperature-dependent correction to the 
    wavelength of the R1 Ruby fluorescence line, based on Datchi et al. (2007), Eq. 2.

    Args:
        T_K (float or np.ndarray): Temperature in Kelvin.

    Returns:
        float or np.ndarray: Temperature correction in nanometers.

    References:
        Datchi, F., et al. (2007). Optical pressure sensors for high-pressure–high-temperature 
        studies in a diamond anvil cell. High Pressure Research, 27, 447–463.
    """
    d_T = T_K - 296.0
    return 0.00746 * d_T - 3.01e-6 * d_T**2 + 8.76e-9 * d_T**3

def borate_T_corr(T_K: float | np.ndarray) -> float | np.ndarray:
    """
    Computes the temperature correction for the SrB4O7:Sm2+ pressure scale.

    This function calculates the temperature-dependent correction to the 
    wavelength of the SrB4O7:Sm2+ fluorescence line, based on Datchi et al. (2007), Eq. 9.

    Args:
        T_K (float or np.ndarray): Temperature in Kelvin.

    Returns:
        float or np.ndarray: Temperature correction in nanometers.

    References:
        Datchi, F., et al. (2007). Optical pressure sensors for high-pressure–high-temperature 
        studies in a diamond anvil cell. High Pressure Research, 27, 447–463.
    """
    d_T = T_K - 296.0
    return -8.7e-5 * d_T + 4.62e-6 * d_T**2 - 2.38e-9 * d_T**3

def pressure_sensor(
    lbd_nm: float | np.ndarray,
    T_K: float | np.ndarray,
    scale: str = "ruby",
    lbd_0_ruby: float = 694.22,
    lbd_0_borate: float = 685.44,
    borate_scale: str = "datchi2007"
) -> float | np.ndarray:
    """
    Converts the wavelength of fluorescence lines into pressure in GPa.

    This function calculates pressure using either the Ruby R1 fluorescence line 
    or the SrB4O7:Sm2+ fluorescence line, with optional temperature corrections.

    Args:
        lbd_nm (float or np.ndarray): Measured wavelength of the fluorescence line in nanometers.
        T_K (float or np.ndarray): Temperature in Kelvin.
        scale (str, optional): Pressure scale to use ('ruby' or 'borate'). Default is 'ruby'.
        lbd_0_ruby (float, optional): Reference wavelength of the R1 Ruby line at ambient conditions. Default is 694.22 nm.
        lbd_0_borate (float, optional): Reference wavelength of the SrB4O7:Sm2+ line at ambient conditions. Default is 685.44 nm.
        borate_scale (str, optional): Borate pressure scale to use ('datchi2007' or 'rashchenko2015'). Default is 'datchi2007'.

    Returns:
        float or np.ndarray: Pressure in GPa.

    Raises:
        ValueError: If an invalid scale or borate_scale is specified.

    Notes:
        - The Ruby pressure scale follows Shen et al. (2020).
        - The Borate pressure scales are based on Datchi et al. (2007) and Rashchenko et al. (2015).
        - Temperature corrections are applied using Datchi et al. (2007). To disable temperature corrections, set `T_K=296.0`.

    References:
        - Datchi, F., et al. (2007). High Pressure Research, 27, 447–463.
        - Rashchenko, S.V., et al. (2015). Journal of Applied Physics, 117, 145902.
        - Shen, G., et al. (2020). High Pressure Research, 40, 299–314.
    """
    if scale == "ruby":
        lbd_T_corr = lbd_nm - ruby_T_corr(T_K)
        lbd_ratio = (lbd_T_corr-lbd_0_ruby)/lbd_0_ruby
        P_ = 1.87e3 * lbd_ratio * (1.0 + 5.63*lbd_ratio)
        return P_
    elif scale == "borate":
        lbd_T_corr = lbd_nm - borate_T_corr(T_K)
        dlbd = lbd_T_corr-lbd_0_borate
        if borate_scale == "datchi2007":
            A = 3.989
            B = 0.006915
            C = 0.0166
            P_ = A * dlbd * ((1.0 + B*dlbd)/(1.0 + C*dlbd))
        elif borate_scale == "rashchenko2015":
            part1 = (1.0 + dlbd / 685.51)**(-14.3)
            P_ = -2836.0 / 14.3 * (part1 - 1.0)
        else:
            raise ValueError("Choose the scale: datchi2007 or rashchenko2015")
        return P_
    else:
        raise ValueError("Choose the scale: ruby or borate")