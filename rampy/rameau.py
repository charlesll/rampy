# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################

import numpy as np
import pandas as pd
import scipy
import rampy as rp

class rameau:
    """treat Raman spectra of glass to retrieve the glass water content

    Parameters
    ==========
    data_liste: Pandas dataframe
        A Pandas dataframe containing the data and various meta information.

    Attributes
    ==========
    x: ndarray
        a 1D array (Nx1) containing the common x axis (wavelength) of the spectra.
    y: ndarray
        a NxM array (with M the number of spectra) containing the raw intensities of the spectra.
    y_corr: ndarray
        a NxM array (with M the number of spectra) containing the corrected intensities of the spectra.
    y_base: ndarray
        a NxM array (with M the number of spectra) containing the backgrounds of the spectra.
    rws: ndarray
        a 1D array (Nx1) containing the ratio between the integrated intensities of the water and silicate signals.
    rw: ndarray
        a 1D array (Nx1) containing the integrated intensities of the water signal.
    rs:
        a 1D array (Nx1) containing the integrated intensities of the silicate signals.
    water:
        the known glass water content provided in data_liste (set to 0 if predicting for unknowns)
    water_predicted:
        the predicted glass water content provided in data_liste (set to 0 if predicting for unknowns)
    p: ndarray
        calibration coefficient(s) of the LL2012 or DG2017 method
    names: pandas dataframe
        filenames indicated in the data_liste input 

    Note
    ====
    Uses either the LL2012 method (Le Losq et al., 2012) or the DG2017 (Di Genova et al., 2017) method. See references.

    In the LL2012 method, a cubic spline is fitted to the regions of interest provided in self.data_liste (see example).
    The spline is smoothed by the spline_coeff of the data_reduction method. The water content is calculated following eq. (3) of LL2012, with the A coefficient either provided or calculated by the method self.calibrate().

    In the DG2017 method, a third-order polynomial is fitted to the spectra following the instructions of Di Genova et al. (2017).
    The water content is calculated as wt% H2O = Rws * (a * [FeO wt%] + b) with a and b the coefficients either provided or calculated by the method self.calibrate().

    References
    ==========
    LL2102: C. Le Losq, D. R. Neuville, R. Moretti, J. Roux, Determination of water content in silicate glasses using Raman spectrometry: Implications for the study of explosive volcanism. American Mineralogist. 97, 779–790 (2012).
    DG 2017 D. Di Genova et al., Effect of iron and nanolites on Raman spectra of volcanic glasses: A reassessment of existing strategies to estimate the water content. Chemical Geology. 475, 76–86 (2017).


    """
    # class variables
    x = []
    y = []
    y_corr = []
    y_base = []
    rws = []
    rw = []
    rs = []
    water = []
    water_predicted = []
    p = []
    names = []

    def __init__(self,data_liste):
        self.data_liste = data_liste
        self.water = np.asarray(data_liste["Water, wt%"])
        self.names = data_liste["Name"]

        try: # we test if gcvspline is installed
            import gcvspline
        except ImportError:
            print('gcvspline module not found. Please install it and check it works (requires a FORTRAN compiler).')

    def data_reduction(self,method="LL2012",delim='\t',path_in='./raw/',laser=514.532,spline_coeff=0.001,poly_coeff=3):
        """process Raman spectra of glass to calculate the Rws ratio

        Parameters
        ==========
        self: object
            a rameau object that has been initiated.
        method: string
            The used method. LL2012: Le Losq et al. (2012); DG2017: Di Genova et al. (2017). See references. Default = "LL2012".
        delim: string
            File delimiter. Use '\t' for tabulated text or ',' for comma separated text. Default = '\t'.
        path_in: string
            Path for the spectra. Default = './raw/'
        laser: float
            Laser line wavelength in nm. Default = 514.532.
        spline_coeff: float
            Smoothing coefficient for the spline baseline. An array of size len(data_liste) can be provided. Default = 0.001.
        poly_coeff: int
            Polynomial coefficient for the polynomial baseline function. Default = 3 (DG2017 method; set to 2 for Behrens et al. (2006) method).

        Returns
        =======
        self.x: ndarray
            Common x axis.
        self.y_all: ndarray
            All raw spectra from data_liste in an array of length len(x) and with as many column as spectra.
        self.y_all_corr: ndarray
            All corrected spectra from data_liste in an array of length len(x) and with as many column as spectra.
        self.y_all_base: ndarray
            All baselines for spectra from data_liste in an array of length len(x) and with as many column as spectra.
        self.rws: ndarray
            The ratio of the water integrated intensity over that of silicate signals.
        self.rw: ndarray
            The integrated intensity of water signal.
        self.rs: ndarray
            The integrated intensity of silicate signals.
        """
        self.x, self.y, self.y_corr, self.y_base, self.rws, self.rw, self.rs = fit_spectra(self.data_liste,method=method,delim=delim,path_in=path_in,spline_coeff=spline_coeff,poly_coeff=poly_coeff)

    def calibrate(self,method="LL2012"):
        """Fit a calibration by optimizing the K coefficient(s)

        Parameters
        ----------
        self: object
            rameau object with treated spectra (see data_reduction method)
        method: string
            the method used; choose between "LL2012" or "DG2017", default = "LL2012".

        Returns
        -------
        popt: ndarray or float
            The optimized parameter(s);
            if method = "DG2017", popt=np.array([a,b]), parameters of the equation K = a * [FeO wt%] + b.
            if method = "LL2017", popt = A (float), with A parameter in the equation (3) of Le Losq et al. (2012).
        """
        dictio = {"water": self.water,
         "feo": np.asarray(self.data_liste["FeO"]),
         "rws": self.rws}
        try:
            if method == "LL2012":
                self.p = LL2012_calibrate(dictio)
            elif method == "DG2017":
                self.p = DG2017_calibrate(dictio)
        except:
            raise TypeError('Set method: Choose between LL2012 and DG2017 methods.')

    def predict(self,method="LL2012"):
        """predict the water content from the Rws

        Parameters
        ----------
        self: object
            rameau object with treated spectra (see data_reduction method).
        method: string
            the method used; choose between "LL2012" or "DG2017", default = "LL2012".

        Returns
        -------
        H2O
            The glass water contents in wt%
        """
        dictio = {"feo": np.asarray(self.data_liste["FeO"]),
         "rws": self.rws}

        if method == "LL2012":
            try:
                if self.p.size > 0:
                    print("Using adjusted A coefficient: %f" % self.p)
                    self.water_predicted = LL2012_predict(dictio,A=self.p)
                else:
                    print("Using the A coefficient = 0.007609 from Le Losq et al. 2012...")
                    self.water_predicted = LL2012_predict(dictio)
            except:
                raise TypeError("Bad p coefficient. Did you try predicting values with a different method than the calibration?")
        elif method == "DG2017":
            try:
                if self.p.size > 0:
                    print("Using adjusted p coefficients: %f and %f" % (self.p[0],self.p[1]))
                    self.water_predicted = DG2017_predict(dictio,a=self.p[0], b=self.p[1])
                else:
                    print("Using the p coefficients 0.096 and 0.663 from Di Genova et al. 2017...")
                    self.water_predicted = DG2017_predict(dictio)
            except:
                raise TypeError("Bad p coefficients. Did you try predicting values with a different method than the calibration?")

def fit_spectra(data_liste,method="LL2012",delim='\t',path_in='./raw/',laser=514.532,spline_coeff=0.001, poly_coeff=3):
    """Calculate the ratios of water and silicate signals from Raman spectra

    Parameters
    ----------
    data_liste: Pandas DataFrame
        Contains the list of spectra, see provided file as an example
    method: string
        The used method. LL2012: Le Losq et al. (2012); DG2017: Di Genova et al. (2017). See references.
    delim: string
        File delimiter. Use '\t' for tabulated text or ',' for comma separated text.
    path_in: string
        Path for the spectra
    laser: float
        Laser line wavelength in nm
    spline_coeff: float
        Smoothing coefficient for the spline baseline. An array of size len(data_liste) can be provided. Default = 0.001.
    poly_coeff: int
        Polynomial coefficient for the polynomial baseline function. Default = 3 (DG2017 method). Set to 2 for Behrens et al. (2006) method.

    Returns
    -------
    x: ndarray
        Common x axis.
    y_all: ndarray
        All raw spectra from data_liste in an array of length len(x) and with as many column as spectra.
    y_all_corr: ndarray
        All corrected spectra from data_liste in an array of length len(x) and with as many column as spectra.
    y_all_base: ndarray
        All baselines for spectra from data_liste in an array of length len(x) and with as many column as spectra.
    rws: ndarray
        The ratio of the water integrated intensity over that of silicate signals.
    rw: ndarray
        The integrated intensity of water signal.
    rs: ndarray
        The integrated intensity of silicate signals.

    Raises
    ------
    IOError
        If method is not set to LL2012 or DG2017.

    References
    ----------
    C. Le Losq, D. R. Neuville, R. Moretti, J. Roux, Determination of water content in silicate glasses using Raman spectrometry: Implications for the study of explosive volcanism. American Mineralogist. 97, 779–790 (2012).
    D. Di Genova et al., Effect of iron and nanolites on Raman spectra of volcanic glasses: A reassessment of existing strategies to estimate the water content. Chemical Geology. 475, 76–86 (2017).
    """

    import gcvspline
    
    x_all_lf = np.arange(50,1400,1.0)
    x_all_hf = np.arange(2800,3800,1.0)
    x = np.hstack((x_all_lf,x_all_hf))
    y_all = np.zeros((len(x),len(data_liste)))
    y_all_base = np.copy(y_all)
    y_all_corr = np.copy(y_all)

    rws = np.ones(len(data_liste))
    rw = np.ones(len(data_liste))
    rs = np.ones(len(data_liste))

    record_std = np.zeros((len(data_liste),2))

    rois = data_liste.loc[:,"ROI1 lb":"ROI6 hb"]

    for i in range(len(data_liste)):

        # importing the spectra
        sp = np.genfromtxt("./raw/"+data_liste["Name"][i],delimiter=delim,skip_header=1)

        # constructing an interpolator: this will allow an output of all data with the same X axis
        f = scipy.interpolate.interp1d(sp[:,0], sp[:,1],fill_value="extrapolate")

        # temperature and excitation line correction (see Rameau help)
        x, y_all[:,i], sdf = rp.tlcorrection(x,f(x),23.0,laser,normalisation='intensity')

        # getting the roi
        roi = np.array(rois.loc[i]).reshape(int(len(rois.loc[i])/2),2)

        # calculating baseline
        if method == "LL2012": # spline
            
            try:
                c_hf, b_hf = rp.baseline(x,y_all[:,i],roi,"gcvspline",s=spline_coeff)
            except:
                break

            y_all_corr[:,i]=c_hf[:,0]
            y_all_base[:,i]=b_hf[:,0]

        elif method == "DG2017": # polynomial 3 following DG2017 method

            # getting portion of interrest
            x_lf = x[np.where(x<2000.)].reshape(-1)
            x_hf = x[np.where(x>2000.)].reshape(-1)

            y_lf = y_all[np.where(x<2000.),i].reshape(-1)
            y_hf = y_all[np.where(x>2000.),i].reshape(-1)

            c_lf, b_lf = rp.baseline(x_lf,y_lf,np.array([[0,200],[1240,1500]]),"poly",polynomial_order = poly_coeff)
            c_hf, b_hf = rp.baseline(x_hf,y_hf,np.array([[2500,3100],[3750,3900]]),"poly",polynomial_order = poly_coeff)

            y_all_corr[:,i] = np.hstack((c_lf.reshape(-1),c_hf.reshape(-1)))
            y_all_base[:,i] = np.hstack((b_lf.reshape(-1),b_hf.reshape(-1)))

        else:
            raise TypeError('method should be set to LL2012 or DG2017')

        # Area / Integrated Intensity calculation
        S = np.trapz(y_all_corr[np.where((x>150)&(x<1250)),i],x[np.where((x>150)&(x<1250))])
        W = np.trapz(y_all_corr[np.where((x>3100)&(x<3750)),i],x[np.where((x>3100)&(x<3750))])

        # updating the Pandas dataframe rws
        rs[i] = S[0]
        rw[i] = W[0]
        rws[i] = W[0]/S[0]

    return x, y_all, y_all_corr, y_all_base, rws, rw, rs

def DG2017_predict(dictio,a=0.096,b=0.663):
    """Calculate the K coefficient for the DG2017 method.

    Parameters
    ----------
    dictio: dict
        a dictionary with ndarrays named "feo" and "rws"
    a and b: float
        factors in the equation: K = a * [FeO wt%] + b; default values from Di Genova et al. (2017)

    Returns
    -------
    H2O (wt %): ndarray
        The water content of the glasses calculated as Rws * (a * [FeO wt%] + b)

    """
    return dictio["rws"]*(dictio["feo"]*a + b)

def DG2017_calibrate(dictio):
    """Fit a calibration by optimizing the K coefficient in the DG2017 method

    Parameters
    ----------
    dictio: dictionary
        dictionary with arrays named "feo", "rws" and "water".

    Returns
    -------
    popt: ndarray
        The optimize a and b parameters of the equation K = a * [FeO wt%] + b.
    """
    popt, pcov = scipy.optimize.curve_fit(DG2017_predict, dictio, dictio["water"])
    return popt

def LL2012_predict(dictio,A=0.007609):
    """Predict the water content using the equation (3) from Le Losq et al. (2012)

    Parameters
    ----------
    dictio: dict
        a dictionary with ndarray named "rws"

    Returns
    -------
    H2O
        The glass water contents in wt%
    """
    return 100*A*dictio["rws"]/(1 + A*dictio["rws"])

def LL2012_calibrate(dictio):
    """Fit a calibration line following equations (2) and (3) from Le Losq et al. (2012)

    Parameters
    ----------
    dictio
        dictionary with arrays named "feo", "rws" and "water".

    Returns
    -------
    A: float
        The parameter in the equation (3) of Le Losq et al. (2012).

    """
    A, pcov = scipy.optimize.curve_fit(LL2012_predict,dictio,dictio["water"])
    return A
