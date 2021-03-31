#!/usr/bin/env python
#-*- coding: utf-8 -*-
import numpy as np
import rampy as rp

def peakarea(shape,**options):
    """returns the area of a peak

    (!experimental!)

    gaussian peak area is calculated analytically; areas for other peak shapes
    are calculated using trapezoidal integration.

    Inputs
    ------
    shape : string
        gaussian, lorentzian, pseudovoigt or pearson7

    Options
    -------
    amp : float or ndarray
        amplitude of the peak
    pos : float or ndarray
        peak position
    HWHM : float or ndarray
        half-width at half-maximum
    a3 : float or ndarray
        a3 parameters for pearson7
    eseAmplitude : float or ndarray
        standard deviation on amp; Default = None
    eseHWHM : float or ndarray
        standard deviation on HWHM; Default = None

    Returns
    -------
    area : float or ndarray
        peak area
    esearea : float or ndarray
        error on peak area; will be None if no errors on amp and HWHM were provided.
    """
    if shape == "gaussian":
        try:
            amp=options.get("amp")
            HWHM=options.get("HWHM")
            area = np.sqrt(np.pi/np.log(2))*amp*HWHM
            if options.get("eseAmplitude") != None:
                eseAmplitude = options.get("eseAmplitude")
                if options.get("eseHWHM") != None:
                    eseHWHM = options.get("eseHWHM")
                    esearea = np.sqrt((np.pi/np.log(2)*HWHM)**2 * eseAmplitude**2 + (np.pi/np.log(2)*amp)**2 * eseHWHM**2)
                return area, esearea
            else:
                return area

        except ValueError:
            print("amp and HWHM should be provided")

    elif (shape == "lorentzian") or (shape == "pseudovoigt") or (shape == "pearson7"): # other shapes: trapezoidal integration
        try:
            amp=options.get("amp")
            pos=options.get("pos")
            HWHM=options.get("HWHM")
        except ValueError:
            print("amp,pos,HWHM should be provided")

        # we construct a fake X axis with lots of samples
        x_int = np.linspace(pos-(10*HWHM),pos+10*HWHM, 10000)

        # now calculate y, also ask for any additional parameters
        if shape == "pearson7":
            try:
                a3=options.get("a3")
            except ValueError:
                print("a3 should be provided for pearson7 peak")
            y = rp.pearson7(x_int,amp,pos,HWHM,a3)
        elif shape == "pseudovoigt":
            try:
                LGratio=options.get("LGratio")
            except ValueError:
                print("LGratio should be provided for pseudovoigt peak")
            y = rp.pseudovoigt(x_int,amp,pos,HWHM,LGratio)
        else: # for lorentzian peak
            y = rp.lorentzian(x_int,amp,pos,HWHM)

        return np.trapz(y,x_int)
    else:
        raise ValueError("Shape should be gaussian, lorentzian, pseudovoigt or pearson7")


def gaussianarea(amp,HWHM,**options):
    """returns the area of a Gaussian peak

    Inputs
    ------
    amp : float or ndarray
        amplitude of the peak
    HWHM : float or ndarray
        half-width at half-maximum

    Options
    -------
    eseAmplitude : float or ndarray
        standard deviation on amp; Default = None
    eseHWHM : float or ndarray
        standard deviation on HWHM; Default = None

    Returns
    -------
    area : float or ndarray
        peak area
    esearea : float or ndarray
        error on peak area; will be None if no errors on amp and HWHM were provided.
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
