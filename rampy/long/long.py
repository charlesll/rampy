import numpy as np

def longcorr(data,temp,wave): # input are a two column matrix of data, temperature in C, wave as the wavelength in cm-1
    """
    # Long Correction
    # Charles Le Losq
    # CIW Washington 2014
    
    # Long's correction of Raman spectra and normalisation
    # last rev. Oct 2010, converted to Matlab and then to Python
    # ensures strictly increasing values of wavenumber
    # calc. e.s.e. as Long cor. norm. sqrt(n_raw) 3d output col.
    # exp. format to avoid null ese.
    # program long3;
    
    #Original program in pascal from J. Roux, modified for Python C. Le Losq.
    
    # See Shucker and Gammon, Phys. Rev. Lett. 1970; Galeener and Sen, Phys. Rev. B 1978; Neuville and Mysen, GCA 1996; Le Losq et al. AM 2012 and GCA 2014 for equations and theory
    """
    h = 6.62606896*10**-34   # J.s    Plank constant
    k = 1.38066e-23;     # J/K    Boltzman
    c = 2.9979e8;        # m/s    Speed of light
    v = wave;            # cm-1   Excitating laser line
    nu0 = 1.0/v*10**9;    # conversion of v in m
    T = temp + 273.15;   # C input temperature in K

    x = data[:,0]
    y = data[:,1]
    
    
    # Calculate the error on data as sqrt(y). If y <= 0, then error = sqrt of absolute value of y.
    # for doing that we construct an artificial array y containing only positive values 
    
    ese = np.sqrt(np.absolute(data[:,1])) # we assume that errors on raw data are in sqrt(y)
    
    # For retrieving the good errors after long correction, one simple way is
    # to work with relative errors...
    error = ese/y;
    
    # then we proceed to the correction (Neuville and Mysen, 1996; Le Losq et
    # al., 2012)
    
    nu = 100.0*x; # cm-1 -> m-1 Raman shift
    rnu = nu0-nu; # nu0 is in m-1
    t0 = nu0*nu0*nu0*nu/rnu/rnu/rnu/rnu;
    t1 = -h*c*nu/k/T; # c in m/s  : t1 dimensionless
    t2 = 1 - np.exp(t1);
    longsp = y*t0*t2; # for y values
    #long2 = ese*t0*t2; # for errors, as comment as we use relative errors
    
    
    # normalized to max intensity
    # tried area with np.trapz but their is an issue for now
    norm = np.trapz(longsp,x)
    #norm = np.max(longsp)
    longsp = longsp/norm
    
    eselong = error*longsp
    
    spectreout = np.zeros((len(x),3))
    spectreout[:,0] = x
    spectreout[:,1] = longsp
    spectreout[:,2] = eselong
    
    return spectreout