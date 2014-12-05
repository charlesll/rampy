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
sys.path.append("/Users/charleslelosq/Documents/RamPy/lib-charles/")
sys.path.append("/Users/charleslelosq/Documents/RamPy/lib-charles/gcvspl/")

import numpy as np
import scipy
import matplotlib

import matplotlib.gridspec as gridspec
from pylab import *
from scipy import interpolate
from scipy.optimize import curve_fit

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename
     
# Home made modules
from spectratools import *

#### DATA PATHS AND INPUT
tkMessageBox.showinfo(
            "Open ",
            "Please open the list of spectra")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
samplename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

# we import the information in an array, skipping the first line
dataliste = np.genfromtxt(samplename,dtype = 'string', skip_header=0,skip_footer=0)
pathdiamond = (dataliste[:,0])
pathsample = (dataliste[:,1])
pathsave = (dataliste[:,2])

output = np.zeros((len(pathdiamond),2))

###############################################################################
######## HERE ARE DEFINED THE BONDARIES / BIRS USED IN THE CODE ###############
###############################################################################
# SAME FOR ALL SPECTRA
# DEFINED HERE TO SIMPLIFY THE USE OF THIS CODE

x = np.arange(2035,3800,0.2) # X scale for correction of X deviations, take care to put correct bondaries

# Background Interpolation Regions for baseline subtraction
birDiamond = np.array([(2035,2055),(2810,3849)]) # BIR diamond
smofd = 0.03 # smoothing factor
birSample = np.array([(2035,2055),(2840,2860),(3010,3030),(3790,3849)]) # BIR sample: typical 2860-2950 for melt; For fluid, may add (3000,3040)
smofs = 0.05 # smoothing factor

#### DO YOU NEED A FILTER BEFORE ADJUSTING DIAMOND/SAMPLE SIGNALS???
filterSwitch1 = 0 #0 to switch off, 1 to turn on
cutfq = np.array([0.05]) # cutoff Frequency
filter1type = 'low' # 'low', 'high', 'bandstop', 'bandpass', change cutfq as a function of the filter type (see doc)

bmaxd = np.array([(2400,2500)]) # 2400 2475 Here is were the program have to search for the frequency of the peak that is used for x calibration
dconvfit = np.array([(2150,2360)]) # Here is were we calculate the conversion factor between the diamond and sample spectra. Unless you know what you are doing, DO NOT MODIFY!

#### DO YOU NEED A SECOND FILTER BEFORE AREA CALCULATION, IN CASE OF "INTERFRINGEANCE" BANDS FOR INSTANCE?
filterSwitch2 = 0 #0 to switch off, 1 to turn on
cutfq2 = np.array([0.013]) # cutoff Frequency
filter2type = 'low' # 'low', 'high', 'bandstop', 'bandpass', change cutfq2 as a function of the filter type (see doc)

# Here are the bondaries for the calculation of the areas. Please CHANGE THEM AS A FUNCTION OF YOUR DATA!!!!!!
lb = 2200 # Lower Bondary
hb = 3800 # Upper Bondary
mb1 = 2850 # Intermediate 1 (end of OD peak) good for fluid at 2860 (D2 doublet near 2900 to avoid), for melt at 2800, depends on temperature
mb2 = 3030 # Intermediate 2 (beginning of OH peak) set as mb1 for melt, but for fluid at 3000 it allows avoiding the D2 doublet near 2900 cm-1

# In case some of the signal is in the negative portion at the end, you want to activate that:
birFinalSwitch = 0 # If 0, deactivated, if 1, activated
birFinal = np.array([(2036,2150),(2840,2850),(3777,3779)])
basetype = 'poly'
smofu = 2 #final smo/poly factor


for i in range(len(pathdiamond)): # We loop over in dataliste

    rawsample = np.genfromtxt(pathsample[i],skip_header=20, skip_footer=43) # Skipping lines from JASCO files
    diamond = np.genfromtxt(pathdiamond[i],skip_header=20, skip_footer=43)

    #### SUBTRACTION OF THE DIAMOND SIGNAL FROM SAMPLE SPECTRA
    # WE INTERPOLATE SPECTRA WITH SAME X in case its is not the case
    tck = interpolate.splrep(diamond[:,0],diamond[:,1],s=0)
    y = interpolate.splev(rawsample[:,0],tck,der=0) # interp with x from rawsample
    rawdiamond = np.zeros(shape(rawsample))
    rawdiamond[:,1] = y
    rawdiamond[:,0] = rawsample[:,0]

    # If needed to suppress high-frequency oscillations, here a filter:
    if filterSwitch1 == 1:
        rawsample = spectrafilter(rawsample,filter1type,cutfq,1,np.array([1]))
        rawdiamond = spectrafilter(rawdiamond,filter1type,cutfq,1,np.array([1]))

    # FOR GCVSPL: errors = sqrt(y), directly calculated in spectratools.linbaseline
    corrdiamond, baselineD, coeffsD = linbaseline(rawdiamond,birDiamond,'gcvspline',smofd) # SUbtract a  baseline below Diamond spectra  
    corrsample, baselineS, coeffsS = linbaseline(rawsample,birSample,'gcvspline',smofs)   # SUbtract a  baseline below Sample spectra
    
    #### CORRECTION OF ANY X SHIFTS    

    # We look at the Raman shift of the maximum intensity second order diamond peak near 2444-2475 cm-1
    # If needed please use the 2160 cm-1 peak and change values below!!!!! (see intro)
    # We search the maximum values in the spectra
    maxDiamond = corrdiamond[np.where((corrdiamond[:,0] > bmaxd[0,0]) & (corrdiamond[:,0] < bmaxd[0,1]))] 
    maxD = np.max(maxDiamond[:,1])
    maxSample = corrsample[np.where((corrsample[:,0] > bmaxd[0,0]) & (corrsample[:,0] < bmaxd[0,1]))]
    maxS = np.max(maxSample[:,1])
    # We retrieve the corresponding x value
    bebe = corrsample[np.where((corrsample[:,1]==maxS))]
    cece = corrdiamond[np.where((corrdiamond[:,1]==maxD))]
    # we calculate the correction factor
    corrx = bebe[:,0] - cece[:,0]
    # To apply the correction for x shift, we need to interpolate to create new datasets
    tck = interpolate.splrep(corrsample[:,0],corrsample[:,1],s=0)
    y = interpolate.splev(x,tck,der=0) # sample spectrum
    tck = interpolate.splrep(corrdiamond[:,0]+corrx,corrdiamond[:,1],s=0)
    y2 = interpolate.splev(x,tck,der=0) # diamond spectrum

    # The following rows contain the spectra corrected from x shifts
    diamondfinal = np.zeros((len(x),2))
    samplefinal = np.zeros((len(x),2))
    diamondfinal[:,0] = x
    samplefinal[:,0] = x
    diamondfinal[:,1] = y2
    samplefinal[:,1] = y

    ##### CALCULATION OF THE INTENSITY CONVERSION FACTOR DIAMOND-SAMPLE

    # We use least-square fit to find this conversion factor, between 2100 and 2360 cm-1
    # Their is few to no expected D2O-OD signals in this region, and two nice peak from diamond
    DiamondtoFit = diamondfinal[np.where((diamondfinal[:,0]> dconvfit[0,0]) & (diamondfinal[:,0] < dconvfit[0,1]))]
    SampletoFit =  samplefinal[np.where((samplefinal[:,0]> dconvfit[0,0]) & (samplefinal[:,0] < dconvfit[0,1]))]
    corrY, cov_out = curve_fit(fun, DiamondtoFit[:,1], SampletoFit[:,1],np.array([(1)])) #Fitting the peak

    # we record anorther spectra for diamond (normalized) and also for sample (subtracted from diamond and
    # area normalized)
    diamondfinal[:,1] = fun(diamondfinal[:,1],corrY) # We correct the diamond spectra
    sampleultimateINT = np.zeros(shape(samplefinal)) # we create the output array containing the good sample spectrum
    sampleultimateINT[:,0] = samplefinal[:,0] 
    sampleultimateINT[:,1] = samplefinal[:,1] - diamondfinal[:,1] #We subtract the second order diamond from sample

    # Last correction, if their is a portion of the "ultimate" signal in the negative domain:
    if birFinalSwitch == 0:
        sampleultimate = sampleultimateINT
        if filterSwitch2 == 1:
            sampleultimate = spectrafilter(sampleultimate,filter2type,cutfq2,1,np.array([1]))
    else:
        if filterSwitch2 == 1:
            sampleultimateINT2 = spectrafilter(sampleultimateINT,filter2type,cutfq2,1,np.array([1]))
            sampleultimate, baselineU, coeffsU = linbaseline(sampleultimateINT2,birFinal,basetype,smofu)
        else:
            sampleultimate, baselineU, coeffsU = linbaseline(sampleultimateINT,birFinal,basetype,smofu)
    
    # Errors as sqrt(corrected y)
    # would be relatively higher than initial errors and may take a part of the correction process into account
    ese = np.zeros(np.shape(sampleultimate))
    ese[:,0] = sampleultimate[:,0]
    ese[:,1] = np.sqrt(abs(samplefinal[:,1]))/samplefinal[:,1]*sampleultimate[:,1]

    # Now we can do whatever we want
    # In the following we calculate the ratio between the D and H peaks
    # Values have to be adjusted manally in the following, because the value that separate OH and OD signals
    # varies! SEE INTRODUCTION... No automatic process also to find it (as for instance the minimum), because the D2 signal
    # mess up the things here for such a process... 
    peakOD = sampleultimate[np.where((sampleultimate[:,0]> lb) & (sampleultimate[:,0] < mb1))]
    peakOH = sampleultimate[np.where((sampleultimate[:,0]> mb2) & (sampleultimate[:,0] < hb))]
    esepeakOD = ese[np.where((ese[:,0]> lb) & (ese[:,0] < mb1))]
    esepeakOH = ese[np.where((ese[:,0]> mb2) & (ese[:,0] < hb))]
    AireOD = np.trapz(peakOD[:,1],peakOD[:,0])
    AireOH = np.trapz(peakOH[:,1],peakOH[:,0])
    eseAireOD = np.trapz(esepeakOD[:,1],esepeakOD[:,0])
    eseAireOH = np.trapz(esepeakOH[:,1],esepeakOH[:,0])
    ratioDH = AireOD/AireOH
    eseratioDH = np.sqrt((1/AireOH)**2*eseAireOD**2+((AireOD-AireOH)/(AireOH**2))**2*eseAireOH**2)
    output[i,0] = ratioDH
    output[i,1] = eseratioDH
    
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
    ax2.plot(corrdiamond[:,0],corrdiamond[:,1],'r-')
    ax2.plot(diamondfinal[:,0],diamondfinal[:,1],'g-')
    
    # Intensity is normalized for representation
    ax3.plot(sampleultimateINT[:,0],sampleultimateINT[:,1]/amax(peakOD[:,1])*100 ,'k-')
    ax3.plot(sampleultimate[:,0],sampleultimate[:,1]/amax(peakOD[:,1])*100 ,'r-')
    if birFinalSwitch == 1:
        ax3.plot(baselineU[:,0],baselineU[:,1]/amax(peakOD[:,1])*100 ,'g-')
    
    # Limits
    ax1.set_xlim(2000,3850)
    ax2.set_xlim(2000,3850)
    ax3.set_xlim(2000,3850)
    
    # we search the lower limit for ax2 and ax3 but the higher free.
    ax2.set_ylim(np.amin(corrdiamond[:,1])-5/100*np.amin(corrdiamond[:,1]),)#np.amax(corrdiamond[:,1])+10/100*np.amax(corrdiamond[:,1])
    ax3.set_ylim(np.amin(sampleultimate[:,1]/amax(peakOD[:,1])*100)-5/100*np.amin(sampleultimate[:,1]/amax(peakOD[:,1])*100),) #np.amax(sampleultimate[:,1]/amax(peakOD[:,1])*100)+10/100*np.amax(sampleultimate[:,1]/amax(peakOD[:,1])*100)
    
    # Labels:
    ax1.set_ylabel("Intensity, a. u.", fontsize = 18, fontweight = "bold")
    ax2.set_xlabel("Raman shift, cm$^{-1}$",fontsize = 18,fontweight = "bold")
    
    plt.tight_layout()

    #### DATA and FIGURE OUTPUT
    np.savetxt(pathsave[i],sampleultimate)
    name = pathsave[i]
    namefig = name[0:name.rfind('.')]+'.pdf'
    savefig(namefig) # save the figure
     
