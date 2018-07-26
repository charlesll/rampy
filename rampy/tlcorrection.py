# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import numpy as np
from scipy.constants import c, h, k

def tlcorrection(x,y,temp,wave, **kwargs):
    """
        tlcorrection(x,y,temp,wave,**kwargs)

    Parameters
    ----------
    x
        Raman shifts in cm-1
    y
        Intensity values as counts
    temp
        Temperature in °C
    wave
        wavenumber of the laser that excited the sample, in nm

    kwargs
    ------
    correction
        Equation used for the correction. Choose between 'long', 'galeener', or 'hehlen'. Default = 'long'.
    normalisation
        Data normalisation procedure. Choose between 'intensity', 'area', or 'no'. Default = 'area'.
    density
        The density of the studied material in kg m-3, to be used with the 'hehlen' equation. Default = 2210.0 (density of silica).

    Returns
    -------
    x
        The Raman shifts values.
    long
        The corrected intensities.
    eselong
        The errors calculated as sqrt(y) on raw intensities and propagated after the correction.

    Remarks
    -------

    This correction uses the formula reported in Galeener and Sen (1978), Mysen et al. (1982), Brooker et al. (1988) and Hehlen et al. (2010).

    The 'galeener' equation is the exact one reported in Galeener and Sen (1978), which is a modification from Shuker and Gammon (1970) for accounting of (vo - v)^4 dependence of the Raman intensity. See also Brooker et al. (1988) for further discussion.

    The 'long' equation is that of Galeener and Sen (1978) corrected by a vo^3 coefficient for removing the cubic meter dimension of the equation of 'galeener'. This equation has been used in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012).

    The 'hehlen' equation is that reported in Hehlen et al. (2010). It actually originates before this publication (Brooker et al. 1988). It uses a different correction that avoid crushing the signal below 500 cm-1. THerefore, it has the advantage of keeping intact the Boson peak signal in glasses.
    """

    if x[-1] < x[0]: # to raise an error if decreasing x values are provided
        raise ValueError('x values should be increasing.')

    # getting the kwargs
    correction = kwargs.get('correction','long')
    normalisation = kwargs.get('normalisation','area')
    density = kwargs.get('density',2210.0)

    #h = 6.626070040e-34   # J S    Plank constant from NIST
    #hb = 1.054571800e-34 # J S    Reduced Plank constant from NIST
    #k = 1.38064852e-23      # J K-1    Boltzman constant from NIST
    #c = 299792458.         # M S-1    Speed of light from NIST
    nu0 = 1.0/wave*1e9     # nu0 laser is in M-1 (wave is in nm)
    T = temp + 273.15    # K temperature
    # density is in KG M-3

    # Calculate the relative error on data as sqrt(y). If y <= 0, then error = abs(y).
    ese = np.sqrt(np.absolute(y))/np.absolute(y) # relative errors

    # get the Raman shift in m-1
    nu = 100.0*x # cm-1 -> m-1 Raman shift


    # then we proceed to the correction
    try:
        if correction == 'long':
            # Formula used in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012) (corrected for using the Planck constant in the last reference)
            # It is that reported in Brooker et al. (1988) with the addition of a scaling nu0^3 coefficient for adimentionality
            frequency = nu0**3*nu/((nu0-nu)**4) # frequency correction; dimensionless
            boltzman = 1.0 - np.exp(-h*c*nu/(k*T)) # temperature correction with Boltzman distribution; dimensionless
            ycorr = y*frequency*boltzman # correction

        elif correction == 'galeener':
            # This uses the formula reported in Galeener and Sen (1978) and Brooker et al. (1988); it uses the Bose-Einstein / Boltzman distribution
            # Formula from  without the scaling vo^3 coefficient reported in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012)
            frequency = nu/((nu0-nu)**4) # frequency correction; M^3
            boltzman = 1.0 - np.exp(-h*c*nu/(k*T)) # temperature correction with Boltzman distribution; dimensionless
            ycorr = y*frequency*boltzman # correction

        elif correction =='hehlen':
            # this uses the formula reported in Hehlen et al. 2010
            frequency = 1/(nu0**3*density) # frequency + density correction; M/KG
            boltzman = 1.0 - np.exp(-h*c*nu/(k*T)) # dimensionless
            ycorr = nu*y*frequency*boltzman # correction
    except:
        print("Not implemented, choose correction = long, galeener or hehlen.")

    # we take care of the normalisation
    try:
        if normalisation == 'area':
            ycorr = ycorr/np.trapz(ycorr,x) # area normalisation
        elif normalisation == 'intensity':
            ycorr = ycorr/np.max(ycorr) # max. intensity normalisation
        elif normalisation == 'no':
            print("No normalisation...")
    except:
        print("Set the optional normalisation parameter to area, intensity or no.")

    esecorr = ese*ycorr # error calculation

    return x, ycorr, esecorr
