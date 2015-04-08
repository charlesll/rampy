# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 14:34:21 2014

Heavily modified on March 3, 2015

@author: charleslelosq
Carnegie Institution for Science

This script calculates pressure (MPa) in DAC experiments from Raman shift of 13C diamonds
Put it anywhere, you need however the lib-charles library as well as numpy, scipy, matplotlib, and Tkinter
that usually come with any python distribution
"""
import sys, os
sys.path.append("/Users/closq/Google Drive/Rampy/lib-charles/")

import numpy as np
import scipy
import matplotlib
import matplotlib.gridspec as gridspec
from pylab import *
from StringIO import StringIO
from scipy import interpolate
from scipy.optimize import curve_fit

# to fit spectra we use the lmfit software of Matt Newville, CARS, university of Chicago, available on the web
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, fit_report

from spectratools import *

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename

### residual function
def residual(pars, x, data=None, eps=None):
    # unpack parameters:
    #  extract .value attribute for each parameter
    a1 = pars['a1'].value
    a2 = pars['a2'].value
    a3 = pars['a3'].value
    
    f1 = pars['f1'].value
    f2 = pars['f2'].value
    f3 = pars['f3'].value
    
    l1 = pars['l1'].value
    l2 = pars['l2'].value
    l3 = pars['l3'].value
    
    e1 = pars['e1'].value
    e2 = pars['e2'].value
    e3 = pars['e3'].value
    
    # Pseuydovoigts peaks
    peak1 = pseudovoigt(x,a1,f1,l1,e1)
    peak2 = pseudovoigt(x,a2,f2,l2,e2)
    peak3 = pseudovoigt(x,a3,f3,l3,e3)    
    
    model = peak1 + peak2 + peak3
    
    if data is None:
        return model, peak1, peak2, peak3
    elif eps is None:
        return (model - data)
    else:
        return (model - data)/eps


#### DATA PATHS AND INPUT
tkMessageBox.showinfo(
            "Open ",
            "Please open the list of spectra")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

datalist = np.genfromtxt(filename,delimiter = '\t', dtype=None,skip_header=0)
RTdiamond = datalist['f0']
HTdiamond = datalist['f1']
temperature = datalist['f2']

# for the results
outputs = np.zeros((len(RTdiamond),4)) 

for lg in range(len(RTdiamond)):
    inputRT = np.genfromtxt(RTdiamond[lg],skip_header=20, skip_footer=43) # get the sample to deconvolute
    inputHT = np.genfromtxt(HTdiamond[lg],skip_header=20, skip_footer=43)

    # quick normalisation of the maximal intensity to 100
    # we will also take care that no intensity goes beyond 0

    inputRT[:,1] = inputRT[:,1]/np.max(inputRT[:,1])*100
    inputHT[:,1] = inputHT[:,1]/np.max(inputHT[:,1])*100

    figure(lg)   
    gs = gridspec.GridSpec(1, 2)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    ax1.set_title('RT')
    ax2.set_title(temperature[lg])
    
    ax1.plot(inputRT[:,0], inputRT[:,1],'k-')
    ax2.plot(inputHT[:,0],inputHT[:,1],'r-')
    params = Parameters()
    #               (Name,  Value,  Vary,   Min,  Max,  Expr)
    params.add_many(('a1',   10,   True,  0,      None,  None),
                    ('f1',   1261,   True, 1260,    1300,  None),
                    ('l1',   26,   True,  0,      50,  None),
                    ('e1',   0.5,   True,  0,      1,  None),
                    ('a2',   100,   True,  0,      None,  None),
                    ('f2',   1332,   True, 1300,    1330,  None),
                    ('l2',   26,   True,  0,      50,  None),
                    ('e2',   0.5,   True,  0,      1,  None),
                    ('a3',   2,   True,  0,      None,  None),
                    ('f3',   1882,   True, 1800,    1920,  None),
                    ('l3',   4,   True,  0,      50,  None),
                    ('e3',   0.5,   True,  0,      1,  None))
                    
    paramsHT = Parameters()
    #               (Name,  Value,  Vary,   Min,  Max,  Expr)
    paramsHT.add_many(('a1',   10,   True,  0,      None,  None),
                    ('f1',   1261,   True, 1260,    1280,  None),
                    ('l1',   26,   True,  0,      50,  None),
                    ('e1',   0.5,   True,  0,      1,  None),
                    ('a2',   100,   True,  0,      None,  None),
                    ('f2',   1315,   True, 1300,    1330,  None),
                    ('l2',   26,   True,  0,      50,  None),
                    ('e2',   0.5,   True,  0,      1,  None),
                    ('a3',   2,   True,  0,      None,  None),
                    ('f3',   1882,   True, 1800,    1920,  None),
                    ('l3',   4,   True,  0,      50,  None),
                    ('e3',   0.5,   True,  0,      1,  None))
                        
     
    # Fitting the RT spectrum
    bir = np.array([(1225,1240),(1370,1500),(1700,1820),(1920,1950)]) # BIR diamond        
    ycorrRT, baselineRT, coeffsDRT = linbaseline(inputRT[:,0:2],bir,'gcvspline',0.15)
    result = minimize(residual, params,method = "leastsq", args=(ycorrRT[:,0], ycorrRT[:,1])) # fit data with  model from scipy
    youRT, peak1RT, peak2RT, peak3RT = residual(params, inputRT[:,0]) # the different peaks
    
    # Fitting the HT spectrum
    ycorrHT, baselineHT, coeffsDHT = linbaseline(inputHT[:,0:2],bir,'gcvspline',0.15)
    resultHT = minimize(residual, paramsHT,method = "leastsq", args=(ycorrHT[:,0], ycorrHT[:,1]))
    youHT, peak1HT, peak2HT, peak3HT = residual(paramsHT, inputHT[:,0]) # the different peaks
        
    ax1.plot(inputRT[:,0], youRT+baselineRT[:,1],'b-',inputRT[:,0],peak1RT,'g-',inputRT[:,0],peak2RT,'g-',inputRT[:,0],peak3RT,'g-',inputRT[:,0],baselineRT[:,1],'g-')
    ax2.plot(inputHT[:,0], youHT+baselineHT[:,1],'b-',inputHT[:,0],peak1HT,'g-',inputHT[:,0],peak2HT,'g-',inputHT[:,0],peak3HT,'g-',inputHT[:,0],baselineHT[:,1],'g-')
    
    Freq13C_RT = params['f1'].value - (params['f3'].value - 1878.28) #f1 and f3 are the frequency of 13C and neon line
    Freq13C_HT = paramsHT['f1'].value - (paramsHT['f3'].value - 1878.28) #f1 and f3 are the frequency of 13C and neon line

    # WE CALCULATE PRESSURE WITH EQ. FROM MYSEN YAMASHITA 2010
    outputs[lg,0] = temperature[lg]
    outputs[lg,1] = Freq13C_RT
    outputs[lg,2] = Freq13C_HT
    outputs[lg,3] = (Freq13C_HT-Freq13C_RT+1.065*10**-2*temperature[lg] + 1.769*10**-5*(temperature[lg]**2)) / 0.002707 #MPa
    
    


