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


"""
import sys
sys.path.append("/Users/charleslelosq/Documents/RamPy/lib-charles/")

import csv
import numpy as np
import scipy
import matplotlib
import matplotlib.gridspec as gridspec
from pylab import *
from StringIO import StringIO
from scipy import interpolate

# to fit spectra we use the lmfit software of Matt Newville, CARS, university of Chicago, available on the web
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
def residual(pars, x, data=None, eps=None):
    # unpack parameters:
    #  extract .value attribute for each parameter
    a1 = pars['a1'].value
    a2 = pars['a2'].value
    
    f1 = pars['f1'].value
    f2 = pars['f2'].value    
    
    l1 = pars['l1'].value
    l2 = pars['l2'].value
    
    # Gaussian model
    
    peak1 = gaussian(x,a1,f1,l1)
    peak2 = gaussian(x,a2,f2,l2)    
    
    model = peak1 + peak2
    
    if data is None:
        return model, peak1, peak2
    if eps is None:
        return (model - data)
    return (model - data)/eps
    
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
    lb = 4700 ### MAY NEED TO AJUST THAT
    hb = 6000
    
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
    #               (Name,  Value,  Vary,   Min,  Max,  Expr)
    params.add_many(('a1',   1,   True, 0,      None,  None),
                    ('f1',   5200,  True, 750,    None,  None),
                    ('l1',   1,  True, 0,      None,  None),
                    ('a2',   1,  True, 0,      None,  None),
                    ('f2',   5400,  True, None,   None,  None),
                    ('l2',   1,  True, None,   None,  None))  
                         
    result = minimize(residual_melt, params, args=(xfit, data)) # fit data with leastsq model from scipy
    model = fit_report(params) # the report
    yout, peak1,peak2,= residual_melt(params,xfit) # the different peaks
    
    #### We just calculate the different areas up to 4700 cmm-1 and those of the gaussians
    # Select interest areas for calculating the areas of OH and H2Omol peaks
    intarea45 = sample[np.where((sample[:,0]> 4100) & (sample[:,0]<4700))]
    area4500 = np.trapz(intarea45[:,1],intarea45[:,0])
    esearea4500 = 1/sqrt(area4500) # We assume that RELATIVE errors on areas are globally equal to 1/sqrt(Area)
      
    # now for the gaussians
    # unpack parameters:
    #  extract .value attribute for each parameter
    a1 = pars['a1'].value
    a2 = pars['a2'].value
    
    l1 = pars['l1'].value
    l2 = pars['l2'].value
    
    AireG1 = gaussianarea(a1,l1)
    AireG2 = gaussianarea(a2,l2)
    
    ##### WE DO A NICE FIGURE THAT CAN BE IMPROVED FOR PUBLICATION
    fig = figure()
    plot(sample[:,0],sample[:,1],'k-')
    plot(xfit,yout,'r-')
    plot(xfit,peak1,'b-')
    plot(xfit,peak2,'b-')
    
    xlim(lb,hb)
    ylim(0,np.max(sample[:,1]))
    xlabel("Wavenumber, cm$^{-1}$", fontsize = 18, fontweight = "bold")
    ylabel("Absorption, a. u.", fontsize = 18, fontweight = "bold")
    
    text(4000,np.max(intarea45[:,1])+0.03*np.max(intarea45[:,1]),('Area OH: \n'+'%.1f' % area4500),color='b',fontsize = 16)
    text(4650,a1 + 0.05*a1,('Area pic 1$: \n'+ '%.1f' % AireG1),color='b',fontsize = 16)
    text(5000,a2 + 0.05*a2,('OH/H$_2$O$_{mol}$: \n'+'%.3f' % ratioOH_H2O+'\n+/-'+'%.3f' % eseratioOH_H2O),color='r',fontsize = 16)
    
    ##### output of data, fitted peaks, parameters, and the figure in pdf
    ##### all goes into the ./deconv/ folder
    name.rfind('/')
    nameout = name[name.rfind('/')+1::]
    namesample = nameout[0:nameout.find('.')]
    pathint = str('/deconv/') # the output folder
    ext1 = '_ydec.txt'
    ext2 = '_params.txt'
    ext3 = '.pdf'
    pathout1 = pathbeg+pathint+namesample+ext1
    pathout2 = pathbeg+pathint+namesample+ext2
    pathout3 = pathbeg+pathint+namesample+ext3
    matout = np.vstack((xfit,data,yout,peak1,peak2))
    matout = np.transpose(matout)
    np.savetxt(pathout1,matout) # saving the arrays of spectra
    fd = os.open( pathout2, os.O_RDWR|os.O_CREAT ) # Open a file and create it if it do not exist
    fo = os.fdopen(fd, "w+") # Now get a file object for the above file.
    fo.write(model) # write the parameters in it
    fo.close()
    savefig(pathout3) # save the figure
 
