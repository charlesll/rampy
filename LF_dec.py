# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 07:54:05 2014

@author: charleslelosq
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 10:19:24 2014

@author: charleslelosq
Carnegie Institution for Science

This script calculates pressure (MPa) in DAC experiments from Raman shift of 13C diamonds
Put it anywhere, you need however the lib-charles library as well as numpy, scipy, matplotlib, and Tkinter
that usually come with any python distribution
"""
import sys
sys.path.append("/Users/charleslelosq/anaconda/lib/python2.7/lib-charles")

import csv
import numpy as np
import scipy
import matplotlib
import matplotlib.gridspec as gridspec
from pylab import *
from StringIO import StringIO
from scipy import interpolate

# to fit spectra we use the lmfit software of Max Newville, CARS, university of Chicago, available on the web
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, fit_report

from spectratools import *  #Charles' libraries and functions

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename

#### We define a set of functions that will be used for fitting data
#### unfortunatly, as we use lmfit (which is convenient because it can fix or release 
#### easily the parameters) we are not able to use arrays for parameters... 
#### so it is a little bit long to write all the things, but in a way quite robust also...
#### gaussian and pseudovoigt functions are available in spectratools
#### if you need a voigt, fix the gaussian-to-lorentzian ratio to 1 in the parameter definition before
#### doing the data fit
def residual_melt(pars, x, data=None, eps=None):
    # unpack parameters:
    #  extract .value attribute for each parameter
    a1 = pars['a1'].value
    a2 = pars['a2'].value
    a3 = pars['a3'].value
    a4 = pars['a4'].value
    a5 = pars['a5'].value
    
    f1 = pars['f1'].value
    f2 = pars['f2'].value
    f3 = pars['f3'].value
    f4 = pars['f4'].value
    f5 = pars['f5'].value    
    
    l1 = pars['l1'].value
    l2 = pars['l2'].value
    l3 = pars['l3'].value
    l4 = pars['l4'].value
    l5 = pars['l5'].value
    
    # Gaussian model
    model = gaussian(x,a1,f1,l1) + \
            gaussian(x,a2,f2,l2) + \
            gaussian(x,a3,f3,l3) + \
            gaussian(x,a4,f4,l4) + \
            gaussian(x,a5,f5,l5)
    
    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model - data)/eps
    
def residual_fluid(pars, x, data=None, eps=None):
    # unpack parameters:
    #  extract .value attribute for each parameter
    a1 = pars['a1'].value
    a2 = pars['a2'].value
    a3 = pars['a3'].value
    a4 = pars['a4'].value
    a5 = pars['a5'].value
    a6 = pars['a6'].value
    
    f1 = pars['f1'].value
    f2 = pars['f2'].value
    f3 = pars['f3'].value
    f4 = pars['f4'].value
    f5 = pars['f5'].value  
    f6 = pars['f6'].value
    
    l1 = pars['l1'].value
    l2 = pars['l2'].value
    l3 = pars['l3'].value
    l4 = pars['l4'].value
    l5 = pars['l5'].value
    l6 = pars['l6'].value
    
    phy = pars['phy'].value
    
    # Gaussian and pseudovoigt model
    model = gaussian(x,a1,f1,l1) + \
            gaussian(x,a2,f2,l2) + \
            pseudovoigt(x,a3,f3,l3,phy) + \
            gaussian(x,a4,f4,l4) + \
            gaussian(x,a5,f5,l5) + \
            gaussian(x,a6,f6,l6)
    
    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model - data)/eps

def peaks_melt(pars, x):
     # unpack parameters:
    #  extract .value attribute for each parameter
    a1 = pars['a1'].value
    a2 = pars['a2'].value
    a3 = pars['a3'].value
    a4 = pars['a4'].value
    a5 = pars['a5'].value
    
    f1 = pars['f1'].value
    f2 = pars['f2'].value
    f3 = pars['f3'].value
    f4 = pars['f4'].value
    f5 = pars['f5'].value    
    
    l1 = pars['l1'].value
    l2 = pars['l2'].value
    l3 = pars['l3'].value
    l4 = pars['l4'].value
    l5 = pars['l5'].value
    
    peak1 = gaussian(x,a1,f1,l1)
    peak2 = gaussian(x,a2,f2,l2)
    peak3 = gaussian(x,a3,f3,l3)
    peak4 = gaussian(x,a4,f4,l4)
    peak5 = gaussian(x,a5,f5,l5)
    
    return peak1,peak2,peak3,peak4,peak5
    
def peaks_fluid(pars, x):
    # unpack parameters:
    #  extract .value attribute for each parameter
    a1 = pars['a1'].value
    a2 = pars['a2'].value
    a3 = pars['a3'].value
    a4 = pars['a4'].value
    a5 = pars['a5'].value
    a6 = pars['a6'].value
    
    f1 = pars['f1'].value
    f2 = pars['f2'].value
    f3 = pars['f3'].value
    f4 = pars['f4'].value
    f5 = pars['f5'].value  
    f6 = pars['f6'].value
    
    l1 = pars['l1'].value
    l2 = pars['l2'].value
    l3 = pars['l3'].value
    l4 = pars['l4'].value
    l5 = pars['l5'].value
    l6 = pars['l6'].value
    
    phy = pars['phy'].value
    
    # Gaussian and pseudovoigt model
    peak1 = gaussian(x,a1,f1,l1)
    peak2 = gaussian(x,a2,f2,l2)
    peak3 = pseudovoigt(x,a3,f3,l3,phy)
    peak4 = gaussian(x,a4,f4,l4)
    peak5 = gaussian(x,a5,f5,l5)
    peak6 = gaussian(x,a6,f6,l6)
    
    return peak1,peak2,peak3,peak4,peak5,peak6


##### CORE OF THE CALCULATION BELOW

#### CALLING THE DATA NAMES
tkMessageBox.showinfo(
            "Open file",
            "Please open the list of spectra")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
with open(filename) as inputfile:
    results = list(csv.reader(inputfile)) # we read the data list

#### LOOP FOR BEING ABLE TO TREAT MULTIPLE DATA
#### WARNING: OUTPUT ARE AUTOMATICALLY GENERATED IN A DIRECTORY CALLED "DECONV"
#### (see end) THAT SHOULD BE PRESENT !!!!!!!!!!
for lg in range(len(results)):
    name = str(results[lg]).strip('[]')
    name = name[1:-1] # to remove unwanted ""
    sample = np.genfromtxt(name) # get the sample to deconvolute

    # we set here the lower and higher bonds for the interest region
    lb = 450
    hb = 1250
    
    interestspectra = sample[np.where((sample[:,0] > lb)&(sample[:,0] < hb))]
    ese0 = interestspectra[:,2]/abs(interestspectra[:,1]) #take ese  as a percentage, we assume that the treatment was made correctly for error determination... if not, please put  sigma = None
    interestspectra[:,1] = interestspectra[:,1]/np.amax(interestspectra[:,1])*100 # normalise spectra to maximum, easier to handle after 
    sigma = abs(ese0*interestspectra[:,1]) #calculate good ese
    #sigma = None # you can activate that if you are not sure about the errors

    xfit = interestspectra[:,0] # region to be fitted
    data = interestspectra[:,1] # region to be fitted

    params = Parameters()
    ####################### FOR MELT:
    ####################### COMMENT IF NOT WANTED
    #           (Name,  Value,  Vary,   Min,  Max,  Expr)
#    params.add_many(('a1',   32,  True, 0,      None,  None),
#               ('f1',   777,  True, 750,    None,  None),
#               ('l1',   016,  True, 0,      None,  None),
#               ('a2',   017,  True, 0,      None,  None),
#               ('f2',   820,  True, None,   None,  None),
#               ('l2',   032,  True, None,   None,  None),  
#               ('a3',   40.51,True, 0,      None,  None),
#               ('f3',   880,  True, None,   None,  None),
#               ('l3',   044,  True, None,   None,  None),  
#               ('a4',   041,  True, 0,      None,  None),
#               ('f4',   1020, True, 900,    None,  None),
#               ('l4',   047,  True, None,   None,  None),  
#               ('a5',   077,  True, 0,      None,  None),
#               ('f5',   1060, True, 1000,   None,  None),
#               ('l5',   037,  True, None,   None,  None))  
    
    ####################### FOR FLUID:
    ####################### COMMENT IF NOT WANTED
    #           (Name,  Value,  Vary,   Min,  Max,  Expr)
    params.add_many(('a1',   032,  True, 0,      None,  None),
                    ('f1',   600,  True, None,   None,  None),
                    ('l1',   020,  True, None,   None,  None),
                    ('a2',   030,  True, 0,      None,  None),
                    ('f2',   700,  True, None,   None,  None),
                    ('l2',   031,  True, None,   None,  None),  
                    ('a3',   80,  True, 0,      None,  None),
                    ('f3',   763,  True, None,   None,  None),
                    ('l3',   010,  True, None,   None,  None),  
                    ('phy',  0.5,  True, 0,       1,    None),
                    ('a4',   041,  True, 0,      None,  None),
                    ('f4',   810,  True, 790,    None,  None),
                    ('l4',   030,  True, None,   None,  None),  
                    ('a5',   020,  True, 0,      None,  None),
                    ('f5',   880,  True, 850,    None,  None),
                    ('l5',   040,  True, None,   None,  None),  
                    ('a6',   030,  True, 0,      None,  None),
                    ('f6',   1060, True, 1000,   None,  None),
                    ('l6',   030,  True, None,   None,  None))  
    
    result = minimize(residual_fluid, params, args=(xfit, data)) # fit data with leastsq model from scipy
    yout = data + result.residual # the model line
    model = fit_report(params) # the report
    peak1,peak2,peak3,peak4,peak5,peak6 = peaks_fluid(params,xfit) # the different peaks
    
    ##### WE DO A NICE FIGURE THAT CAN BE IMPROVED FOR PUBLICATION
    fig = figure()
    plot(xfit,data,'k-')
    plot(xfit,yout,'r-')
    plot(xfit,peak1,'b-')
    plot(xfit,peak2,'b-')
    plot(xfit,peak3,'b-')
    plot(xfit,peak4,'b-')
    plot(xfit,peak5,'b-')
    plot(xfit,peak6,'b-')
    
    xlim(lb,hb)
    ylim(-5,105)
    xlabel("Raman shift, cm$^{-1}$", fontsize = 18, fontweight = "bold")
    ylabel("Normalized intensity, a. u.", fontsize = 18, fontweight = "bold")
    text(650,110,name,fontsize = 18,color='r')     

    ##### output of data, fitted peaks, parameters, and the figure in pdf
    ##### all goes into the ./deconv/ folder
    name.rfind('/')
    nameout = name[name.rfind('/')+1::]
    namesample = nameout[0:nameout.find('.')]
    pathbeg = filename[0:filename.rfind('/')]
    pathint = str('/deconv/') # the output folder
    ext1 = '_ydec.txt'
    ext2 = '_params.txt'
    ext3 = '.pdf'
    pathout1 = pathbeg+pathint+namesample+ext1
    pathout2 = pathbeg+pathint+namesample+ext2
    pathout3 = pathbeg+pathint+namesample+ext3
    np.savetxt(pathout1,np.concatenate((xfit,yout,peak1,peak2,peak3,peak4,peak5),1))
    f = open(pathout2,'w')
    f.write(model)
    f.close()    
    #np.savetxt(pathout2,model)
    savefig(pathout3)
 