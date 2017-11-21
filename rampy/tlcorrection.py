#############################################################################
#Copyright (c) 2017 Charles Le Losq
#
#The MIT License (MIT)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the #Software without restriction, including without limitation the rights to use, copy, #modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, #and to permit persons to whom the Software is furnished to do so, subject to the #following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, #INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Temperature and laser wavelength effects correction for Raman spectra
# Charles Le Losq
# RSES-ANU 2016
# Long's correction of Raman spectra and normalisation
# last rev. Oct 2016 for convertion in Julia
# ensures strictly increasing values of wavenumber
# calc. e.s.e. as Long cor. norm. sqrt(n_raw) 3d output col.
# exp. format to avoid null ese.
#
# THOSE FORMULA ARE WRITTEN FOR STOKES RAMAN. May be easily adapted for anti-Stokes upon request.
#
# References
# Shuker, Reuben, and Robert Gammon. 1970. Raman-Scattering Selection-Rule Breaking and the Density of States in Amorphous Materials. Physical Review Letters 25 (4): 222–25.
# Galeener, F. L., and Sen, P. N. 1978. Theory of the First-Order Vibrational Spectra of Disordered Solids. Physical Review B 17 (4): 1928–33.
# Mysen, B. O., L. W. Finger, D. Virgo, and F. A. Seifert. 1982. Curve-Fitting of Raman Spectra of Silicate Glasses. American Mineralogist 67: 686–95.
# Brooker et al. 1988 Assessment of correction procedures for reduction of Raman spectra. Journal of Raman Spectroscopy 19(2), 71-78.
# Neuville, D. R., and B. O. Mysen. 1996. Role of Aluminium in the Silicate Network: In Situ, High-Temperature Study of Glasses and Melts on the Join SiO₂-NaAl0₂. Geochimica et Cosmochimica Acta 60: 1727–37.
# Le Losq, C., D. R. Neuville, R. Moretti, and J. Roux. 2012. Determination of Water Content in Silicate Glasses Using Raman Spectrometry: Implications for the Study of Explosive Volcanism. American Mineralogist 97 (5-6): 779–90. doi:10.2138/am.2012.3831.
# Hehlen, B. 2010. Inter-Tetrahedra Bond Angle of Permanently Densified Silicas Extracted from Their Raman Spectra. Journal of Physics: Condensed Matter 22 (2): 025401.
#############################################################################

import numpy as np

def tlcorrection(x,y,temp,wave, **kwargs):
    """
        tlcorrection(x,y,temp,wave,**kwargs)
    
    INPUTS:
    
        data: input spectrum with x and y in first and second columns respectively
    
        temp: Float64, the temperature in °C

        wave: Float64, the wavenumber at which the spectrum was acquirred in nm

    OPTIONS:

        'correction': Equation used for the correction. Choose between 'long', 'galeener', or 'hehlen'. Default = 'long'

        'normalisation': data normalisation procedure. Choose between 'intensity', 'area', or 'no'. Default = 'area'

        'density': the density of the studied material in kg m-3, to be used with the 'hehlen' equation. Default = 2210.0 (density of silica)

    OUTPUTS:

    (are combined in one array if only one output name is given)

        x: the x values

        long: the corrected y values

        eselong: the errors calculated as sqrt(y) on raw data and propagated after the correction
	
    NOTES:

    This correction uses the formula reported in Galeener and Sen (1978), Mysen et al. (1982), Brooker et al. (1988) and Hehlen et al. (2010).

    The 'galeener' equation is the exact one reported in Galeener and Sen (1978), which is a modification from Shuker and Gammon (1970) for accounting of (vo - v)^4 dependence of the Raman intensity. See also Brooker et al. (1988) for further discussion.

    The 'long' equation is that of Galeener and Sen (1978) corrected by a vo^3 coefficient for removing the cubic meter dimension of the equation of 'galeener'. This equation has been used in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012).

    The 'hehlen' equation is that reported in Hehlen et al. (2010). It actually originates before this publication (Brooker et al. 1988). It uses a different correction that avoid crushing the signal below 500 cm-1. THerefore, it has the advantage of keeping intact the Boson peak signal in glasses.
    """


    # getting the kwargs
    correction = kwargs.get('correction','long')
    normalisation = kwargs.get('normalisation','area') 
    density = kwargs.get('density',2210.0)

    h = 6.626070040e-34   # J S    Plank constant from NIST
    hb = 1.054571800e-34 # J S    Reduced Plank constant from NIST
    k = 1.38064852e-23      # J K-1    Boltzman constant from NIST
    c = 299792458.         # M S-1    Speed of light from NIST
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
            ycorr = y*frequency*boltzman # correction
    except:
        print "Not implemented, choose correction = long, galeener or hehlen."
    
    # we take care of the normalisation
    try:
        if normalisation == 'area':
            ycorr = ycorr/np.trapz(ycorr,x) # area normalisation
        elif normalisation == 'intensity':
            ycorr = ycorr/np.max(ycorr) # max. intensity normalisation
        elif normalisation == 'no':
            print "No normalisation..."
    except: 
        print "Set the optional normalisation parameter to area, intensity or no."
    
    esecorr = ese*ycorr # error calculation
    
    return x, ycorr, esecorr