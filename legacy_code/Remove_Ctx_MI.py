# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 17:52 2014
Modified 16 nov 2014 for "auto" treatment

@author: charleslelosq
Carnegie Institution for Science

This script is used to subtract the second order of diamond
in the 2000-4000 cm-1 frequency range of Raman spectra from Diamond Anvil Cell
experiments.

Put it anywhere, you need however to properly set up the path of /lib-charles 
and /lib-charles/gcvspl/ libraries (as well as numpy, scipy, matplotlib, 
and Tkinter that usually come with any python distribution)
"""

import sys
sys.path.append("/Users/charles/Tresors/Charles/Software/RamPy/lib-charles/")
sys.path.append("/Users/charles/Tresors/Charles/Software/RamPy/lib-charles/gcvspl/")

import numpy as np
import scipy
import matplotlib

import matplotlib.gridspec as gridspec
from pylab import *
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
# to fit spectra we use the lmfit software of Matt Newville, CARS, university of Chicago, available on the web
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, fit_report
from ast import literal_eval

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename
     
# Home made modules
from spectratools import *

###############################################################################
######## HERE ARE DEFINED THE BONDARIES / BIRS USED IN THE CODE ###############
###############################################################################

#### DATA PATHS AND INPUT
tkMessageBox.showinfo(
            "Open ",
            "Please open the list of spectra")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
samplename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

pathbeg = samplename[0:samplename.rfind('/')] # The path where the raw and treated folders are 

# we import the information in an array, skipping the first line
dataliste = np.genfromtxt(samplename,dtype = 'string',delimiter = '\t', skip_header=0,skip_footer=0)
pathglass = (dataliste[:,0])
smo_glass = dataliste[:,1]
BIRS_glass = dataliste[:,2]
pathctx = (dataliste[:,3])
smo_ctx = dataliste[:,4]
BIRS_ctx = dataliste[:,5]

output = np.zeros((len(pathglass),2))

long_corr = 0
linear_correction = 1

# SAME FOR ALL SPECTRA
# DEFINED HERE TO SIMPLIFY THE USE OF THIS CODE

x = np.arange(80,3950,0.4) # X scale for correction of X deviations, take care to put correct bondaries

ctx_fraction = np.zeros((len(pathglass),1))

for i in range(len(pathglass)): # We loop over in dataliste

    rawdiamond = np.genfromtxt(pathbeg+'/raw/'+pathctx[i],skip_header=1) # Skipping lines from JASCO files
    rawsample = np.genfromtxt(pathbeg+'/raw/'+pathglass[i],skip_header=1)
    
    rawdiamond.view('f8,f8').sort(order=['f0'], axis=0) 
    rawsample.view('f8,f8').sort(order=['f0'], axis=0) 

    birSample = np.array(([1260,2000],[2000,2750],[3800,4000]))
    birDiamond = np.array(([1260,2000],[2000,4000]))
    
    rg1 = np.where((rawsample[:,0]>birSample[0,0]) & (rawsample[:,0]<=birSample[0,1]))
    baseline_sample_1 = np.zeros((rg1[-1][-1]+1,2))
    p1 = np.poly1d(np.polyfit(np.ravel(rawsample[rg1,0]),np.ravel(rawsample[rg1,1]),1)) 
    baseline_sample_1[:,0] = rawsample[(rawsample[:,0]<=birSample[0,1]),0]
    baseline_sample_1[:,1] = np.polyval(p1,rawsample[(rawsample[:,0]<=birSample[0,1]),0])        
    
    corr_trash, baseline_sample_2, coeffs = linbaseline(rawsample[(rawsample[:,0]>birSample[0,1])],birSample[1::,:],'gcvspline',0.1)
    
    rg2 = np.where((rawdiamond[:,0]>birDiamond[0,0]) & (rawdiamond[:,0]<=birDiamond[0,1]))
    baseline_diamond_1 = np.zeros((rg2[-1][-1]+1,2))
    p2 = np.poly1d(np.polyfit(np.ravel(rawdiamond[rg2,0]),np.ravel(rawdiamond[rg2,1]),1))
    baseline_diamond_1[:,0] = rawdiamond[(rawdiamond[:,0]<=birDiamond[0,1]),0]        
    baseline_diamond_1[:,1] = np.polyval(p2,rawdiamond[(rawdiamond[:,0]<=birDiamond[0,1]),0])  
    
    corr_trash2, baseline_diamond_2, coeffs2 = linbaseline(rawdiamond[(rawdiamond[:,0]>birDiamond[0,1])],birDiamond[1::,:],'gcvspline',0.1)
    
    baselineS = np.concatenate((baseline_sample_1,baseline_sample_2))
    baselineD = np.concatenate((baseline_diamond_1,baseline_diamond_2))
    
    corrdiamond = np.zeros(rawdiamond.shape)
    corrsample = np.zeros(rawsample.shape)
    
    corrdiamond[:,0] = rawdiamond[:,0]
    corrsample[:,0] = rawsample[:,0]        
    
    corrdiamond[:,1] = rawdiamond[:,1] - baselineD[:,1]
    corrsample[:,1] = rawsample[:,1] - baselineS[:,1]
    
    # variante for ctx signal:
    # reading the birs
    birDiamond_tuple = literal_eval(BIRS_ctx[i]) # the BIRs
    birDiamond = np.zeros((len(birDiamond_tuple),2))   
    birDiamond[:,:] = birDiamond_tuple[:][:]
    
    # FOR GCVSPL: errors = sqrt(y), directly calculated in spectratools.linbaseline
    f_ctx = np.fromstring(smo_ctx[i], dtype = float, sep = ' ')
    f_gls = np.fromstring(smo_glass[i], dtype = float, sep = ' ')
    corrdiamond, baselineD, coeffsD = linbaseline(rawdiamond,birDiamond,'gcvspline',f_ctx) # SUbtract a  baseline below Diamond spectra  
    
    ###########################################################################
    ########################## CORRECTION OF X SHIFTS #########################
    ###########################################################################
        
    ###########################################################################
    # Here we put some objective function we need
    def residual(pars, spforcorr, sptarget=None):
        # unpack parameters:
        #  extract .value attribute for each parameter
        xs = pars['xshift'].value
        yx = pars['yshift'].value
        yx2 = pars['yshiftcarr'].value
        spcorr = np.zeros((shape(spforcorr)))
    
        # we need to resample the spectra to compare them
        tck = interpolate.splrep(spforcorr[:,0]-xs,spforcorr[:,1]*yx+spforcorr[:,1]**2*yx2,s=0)
        spcorr[:,0] = spforcorr[:,0]    
        spcorr[:,1] = interpolate.splev(spforcorr[:,0],tck,der=0)
        
        if sptarget is None: #in such case we return the corrected spectrum
            return spcorr
        return (spcorr[:,1] - sptarget[:,1])
    ###########################################################################        
        
    # Now we choose the portion of spectra to fit
    DiamondtoFit = corrdiamond[np.where((corrdiamond[:,0]> 655) & (corrdiamond[:,0] < 670))]
    SampletoFit =  corrsample[np.where((corrsample[:,0]> 655) & (corrsample[:,0] < 670))]
   
    # Now we enter the model parameters
    params = Parameters()
    params.add_many(('xshift',   1,   True,  -15,      15,  None),
                    ('yshift',   0,   True, None,    None,  None),
                    ('yshiftcarr',   1e-2,   True, None,    None,  None))
        
    # Now we chose the algorithm and run the optimization
    algorithm = "leastsq"
    result = minimize(residual, params,method = algorithm, args=(DiamondtoFit, SampletoFit))

    cds = residual(result.params,corrdiamond)    
    
    # To apply the correction for x shift, we need to interpolate to create new datasets
    tck = interpolate.splrep(corrsample[:,0],corrsample[:,1],s=0)
    tck2 = interpolate.splrep(cds[:,0],cds[:,1],s=0)
        
    # The following rows contain the spectra corrected from x and y shifts
    diamondfinal = np.zeros((len(x),2))
    samplefinal = np.zeros((len(x),2))
    diamondfinal[:,0] = x
    samplefinal[:,0] = x
    diamondfinal[:,1] = interpolate.splev(x,tck2,der=0)
    samplefinal[:,1] = interpolate.splev(x,tck,der=0)
        
    ###########################################################################
    ############## ADJUSTEMENT OF THE PROPORTION OF CRISTAL ###################
    ###########################################################################    
    
    ###########################################################################
    # Here we put some objective function we need
    def residual2(pars, spforcorr, spctx, sptarget=None):
        ctx_fr = pars['ctx_fr'].value # unpack parameters
        
        measured_calc = np.zeros(spforcorr.shape)
        measured_calc[:,0] = spforcorr[:,0]    
        measured_calc[:,1] = spforcorr[:,1] - ctx_fr * spctx[:,1]

        if sptarget is None: #in such case we return the corrected spectrum
            return measured_calc
        return (measured_calc[:,1] - sptarget[:,1])
    ###########################################################################     
        
    # Portion of spectra without ctx signal
    bir_final = np.array(([52,98],[160,172],[446,474],[577,594],[736,813],[883,928],[1130,4000]))
          
    #First estimation of ctx fraction
    diamond660 = np.max(diamondfinal[(diamondfinal[:,0]> 655) & (diamondfinal[:,0] < 670)])
    sample660 = np.max(samplefinal[(samplefinal[:,0]> 655) & (samplefinal[:,0] < 670)])
    init_ctx = sample660/diamond660
    
    smo_sample, baselineSc, coeffsSc = linbaseline(samplefinal,bir_final,'gcvspline',0.15)
    params2 = Parameters()
    params2.add_many(('ctx_fr',   0.2,   True,  0,      2,  None))
        
    # Now we chose the algorithm and run the optimization
    algorithm = "leastsq"
    result2 = minimize(residual2, params2,method = algorithm, args=(samplefinal,diamondfinal,baselineSc))

    sampleultimateINT = residual2(result2.params,samplefinal,diamondfinal)        
    
    ctx_fraction[i,0] = result2.params['ctx_fr'].value
            
    figure(figsize=(10,6))
    gs = matplotlib.gridspec.GridSpec(1, 3)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    ax1.plot(rawsample[:,0],rawsample[:,1],'k-')
    ax1.plot(rawdiamond[:,0],rawdiamond[:,1],'r-')
    ax1.plot(baselineD[:,0],baselineD[:,1],'b--')
    ax1.plot(baselineS[:,0],baselineS[:,1],'b--')
    
    ax2.plot(corrsample[:,0],corrsample[:,1],'k-')
    ax2.plot(baselineSc[:,0],baselineSc[:,1],'b-')
    ax2.plot(corrdiamond[:,0],corrdiamond[:,1],'r-')
    ax2.plot(diamondfinal[:,0],diamondfinal[:,1],'g-')
    
    # Intensity is normalized for representation
    ax3.plot(sampleultimateINT[:,0],sampleultimateINT[:,1] ,'k-')
    
    # Limits
    ax1.set_xlim(50,3850)
    ax2.set_xlim(50,3850)
    ax3.set_xlim(50,3850)
    
    # Labels:
    ax1.set_ylabel("Intensity, a. u.", fontsize = 18, fontweight = "bold")
    ax2.set_xlabel("Raman shift, cm$^{-1}$",fontsize = 18,fontweight = "bold")
    
    plt.tight_layout()


    
    #### DATA and FIGURE OUTPUT
    np.savetxt(pathbeg+'/python_treated/'+pathglass[i],sampleultimateINT)
    #name = pathsave[i]
    #namefig = name[0:name.rfind('.')]+'.pdf'
    #savefig(pathbeg+namefig) # save the figure
     
