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

def rameau(data_liste,method="LL2012",delim='\t',spline_coeff=0.001,poly_coeff=3):
    """Calculate the ratios of water and silicate signals from Raman spectra

    Parameters
    ----------
    data_liste
        The first parameter.
    method
        The used method. LL2012: Le Losq et al. (2012); DG2017: Di Genova et al. (2017). See references.
    delim
        File delimiter. Use '\t' for tabulated text or ',' for comma separated text.
    spline_coeff
        Smoothing coefficient for the spline baseline. An array of size len(data_liste) can be provided. Default = 0.001.
    poly_coeff
        Polynomial coefficient for the polynomial baseline function. Default = 3 (DG2017 method). Set to 2 for Behrens et al. (2006) method.

    Returns
    -------
    x
        Common x axis.
    y_all
        All raw spectra from data_liste in an array of length len(x) and with as many column as spectra.
    y_all_corr
        All corrected spectra from data_liste in an array of length len(x) and with as many column as spectra.
    y_all_base
        All baselines for spectra from data_liste in an array of length len(x) and with as many column as spectra.
    rws
        The ratio of the water integrated intensity over that of silicate signals.

    Raises
    ------
    IOError
        If method is not set to LL2012 or DG2017.

    References
    ----------
    C. Le Losq, D. R. Neuville, R. Moretti, J. Roux, Determination of water content in silicate glasses using Raman spectrometry: Implications for the study of explosive volcanism. American Mineralogist. 97, 779–790 (2012).
    D. Di Genova et al., Effect of iron and nanolites on Raman spectra of volcanic glasses: A reassessment of existing strategies to estimate the water content. Chemical Geology. 475, 76–86 (2017).
    """
    x_all_lf = np.arange(50,1400,1.0)
    x_all_hf = np.arange(2800,3800,1.0)
    x = np.hstack((x_all_lf,x_all_hf))
    y_all = np.zeros((len(x),len(data_liste)))
    y_all_base = np.copy(y_all)
    y_all_corr = np.copy(y_all)

    rws = pd.DataFrame(columns=['S', 'W', 'Rws'])

    record_std = np.zeros((len(data_liste),2))

    rois = data_liste.loc[:,"ROI1 lb":"ROI6 hb"]

    for i in range(len(data_liste)):

        # importing the spectra
        sp = np.genfromtxt("./raw/"+data_liste["Name"][i],delimiter=delim,skip_header=1)

        # constructing an interpolator: this will allow an output of all data with the same X axis
        f = scipy.interpolate.interp1d(sp[:,0], sp[:,1],fill_value="extrapolate")

        for_tl = f(x) - f(x).min() +0.0001

        # temperature and excitation line correction (see Rameau help)
        x, y_all[:,i], sdf = rp.tlcorrection(x,f(x),23.0,514.532,normalisation='intensity')

        # getting the roi
        roi = np.array(rois.loc[i]).reshape(len(rois.loc[i])/2,2)

        # calculating baseline
        if method == "LL2012": # spline

            c_hf, b_hf = rp.baseline(x,y_all[:,i],roi,"gcvspline",s=spline_coeff)

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

        else: IOError('method should be set to LL2012 or DG2017')

        # Area / Integrated Intensity calculation
        S = np.trapz(y_all_corr[np.where((x>150)&(x<1250)),i],x[np.where((x>150)&(x<1250))])
        W = np.trapz(y_all_corr[np.where((x>3100)&(x<3750)),i],x[np.where((x>3100)&(x<3750))])

        # updating the Pandas dataframe rws
        rws['S'].loc[i] = S[0]
        rws['W'].loc[i] = W[0]
        rws['Rws'].loc[i] = W[0]/S[0]

    return x, y_all, y_all_corr, y_all_base, rws

def DGG2017_direct(x,a=0.096,b=0.663):
    """Calculate the K coefficient for the DG2017 method.

    Parameters
    ----------
    x
        a dictionary with arrays named "feo" and "rws"
    a and b
        factors in the equation: K = a * [FeO wt%] + b; default values from Di Genova et al. (2017)

    Returns
    -------
    H2O (wt %)
        The water content of the glasses calculated as Rws * (a * [FeO wt%] + b)

    """
    return x["rws"]*(x["feo"]*a + b)

def DGG2017_optimize(dataset):
    """Optimize the K coefficient in the DG2017 method

    Parameters
    ----------
    dataset
        a dictionary with arrays named "feo", "rws" and "water"

    Returns
    -------
    popt
        The optimize a and b parameters of the equation K = a * [FeO wt%] + b
    H2O (wt %)
        the glass water content in wt %
    """
    popt, pcov = scipy.optimize.curve_fit(DGG2017_direct, dataset, dataset["water"])
    return popt, DGG2017_direct(dataset,popt[0],popt[1])

def eq3_LL2012(rws,A=0.007609):
    """Predict the water content using the equation (3) from Le Losq et al. (2012)

    Parameters
    ----------
    rws
        The ratio of the water integrated intensity over that of silicate signals

    Returns
    -------
    H2O
        The glass water contents in wt%
    """
    return 100*A*rws/(1 + A*rws)

def eq3_LL2012_optimized(dictio):
    """Fit a calibration line following equations (2) and (3) from Le Losq et al. (2012)

    Parameters
    ----------
    dictio
        Dictionnary with "water" and "rws" keys for the sample water content and rws.

    Returns
    -------
    H2O
        The water content in wt% of the glasses calculated from eq3_LL2012
    A
        The parameter in the equation (3) of Le Losq et al. (2012)

    """
    A, pcov = scipy.optimize.curve_fit(eq3_LL2012,dictio["rws"],dictio["water"])
    pred = eq3_LL2012(dictio["rws"],A)
    return pred, A
