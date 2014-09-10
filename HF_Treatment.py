# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 17:52 2014

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

import numpy
import scipy
import matplotlib

import matplotlib.gridspec as gridspec
from pylab import *
from scipy import interpolate
from scipy.optimize import curve_fit

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfilename
     
# Home made modules
from spectratools import *

####### First thing: subtracting the second order of diamond from the experimental spectra

# Collecting the sample spectrum
tkMessageBox.showinfo(
            "Open file",
            "Please open the SAMPLE spectrum")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
rawsample = np.genfromtxt(filename,skip_header=20, skip_footer=43) # Skipping lines from JASCO files

# Same for diamond
tkMessageBox.showinfo(
            "Open file",
            "Please open the DIAMOND spectrum")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename2 = askopenfilename() # show an "Open" dialog box and return the path to the selected file
diamond = np.genfromtxt(filename2,skip_header=20, skip_footer=43)


# WE INTERPOLATE SPECTRA WITH SAME X in case its is not the case
tck = interpolate.splrep(diamond[:,0],diamond[:,1],s=0)
y = interpolate.splev(rawsample[:,0],tck,der=0) # interp with x from rawsample
rawdiamond = np.zeros(shape(rawsample))
rawdiamond[:,1] = y
rawdiamond[:,0] = rawsample[:,0]

# SUbtract a  baseline below Diamond spectra
birDiamond = np.array([(2037,2064),(2730,3850)])
# FOR GCVSPL: errors = sqrt(y), directly calculated in spectratools.linbaseline
corrdiamond, baselineD, coeffsD = linbaseline(rawdiamond,birDiamond,'gcvspline',0.05)

# SUbtract a  baseline below Sample spectra, please change ROI if needed, check that !!!!
birSample = np.array([(2037,2054),(2780,2805),(3750,3850)])
corrsample, baselineS, coeffsS = linbaseline(rawsample,birSample,'gcvspline',0.05)
# If needed to suppress high-frequency oscillations, here a filter:
cutfq = np.array([0.1])
corrsample = spectrafilter(corrsample,'low',cutfq,1,np.array([1]))
corrdiamond = spectrafilter(corrdiamond,'low',cutfq,1,np.array([1]))

# In the following we correct the x axis from any shift
# We look at the Raman shift of the maximum intensity second order diamond peak near 2460 cm-1
# If needed please use the 2160 cm-1 peak and change values below!!!!!
# We search the maximum values in the spectra
maxDiamond = corrdiamond[np.where((corrdiamond[:,0]> 2410) & (corrdiamond[:,0] < 2480))]
maxD = np.max(maxDiamond[:,1])
maxSample = corrsample[np.where((corrsample[:,0]> 2410) & (corrsample[:,0] < 2480))]
maxS = np.max(maxSample[:,1])
# We retrieve the corresponding x value
bebe = corrsample[np.where((corrsample[:,1]==maxS))]
cece = corrdiamond[np.where((corrdiamond[:,1]==maxD))]
# we calculate the correction factor
corrx = bebe[:,0] - cece[:,0]
# To apply the correctin for x shift, we need to interpolate to create new datasets !!!!!!
x = np.arange(2040,3800,0.2) # The choosen x value
tck = interpolate.splrep(corrsample[:,0],corrsample[:,1],s=0)
y = interpolate.splev(x,tck,der=0) # sample spectrum
# On interpole pour corriger des déviations de x
tck = interpolate.splrep(corrdiamond[:,0]+corrx,corrdiamond[:,1],s=0)
y2 = interpolate.splev(x,tck,der=0) # diamond spectrum

# The following rows contain the spectra corrected from x shifts
diamondfinal = np.zeros((len(x),2))
samplefinal = np.zeros((len(x),2))
diamondfinal[:,0] = x
samplefinal[:,0] = x
diamondfinal[:,1] = y2
samplefinal[:,1] = y


# Now we need to calculate the conversion factor between the diamond and sample spectra
# We use least-square fit to find this conversion factor, between 2100 and 2360 cm-1
# Their is no expected D2O-OD signals in this region, and two nice peak from diamond
DiamondtoFit = diamondfinal[np.where((diamondfinal[:,0]> 2140) & (diamondfinal[:,0] < 2340))]
SampletoFit =  samplefinal[np.where((samplefinal[:,0]> 2140) & (samplefinal[:,0] < 2340))]
corrY, cov_out = curve_fit(fun, DiamondtoFit[:,1], SampletoFit[:,1],np.array([(1)])) #Fitting the peak

# we record anorther spectra for diamond (normalized) and also for sample (subtracted from diamond and
# area normalized)
diamondfinal[:,1] = fun(diamondfinal[:,1],corrY) # We correct the diamond spectra
sampleultimate = np.zeros(shape(samplefinal)) # we create the output array containing the good sample spectrum
sampleultimate[:,0] = samplefinal[:,0] 
sampleultimate[:,1] = samplefinal[:,1] - diamondfinal[:,1] #We subtract the second order diamond from sample

# Errors as sqrt(corrected y)
# would be relatively higher than initial errors and may take a part of the correction process into account
ese = np.zeros(np.shape(sampleultimate))
ese[:,0] = sampleultimate[:,0]
ese[:,1] = np.sqrt(abs(samplefinal[:,1]))/samplefinal[:,1]*sampleultimate[:,1]

# Now we can do whatever we want
# In the following we calculate the ratio between the D and H peaks
# Values have to be adjusted manally in the following, because the value that separate OH and OD signals
# varies! No automatic process also to find it (as for instance the minimum), because the D2 signal
# mess up the things here for such a process... 
# Here are the contrains for the calculation
lb = 2560 #no need to change
hb = 3800 #no need to change
mb1 = 2800 # good for fluid at 2860 (D2 doublet near 2900 to avoid), for melt at 2800, depends on temperature
mb2 = 2800 # set as mb1 for melt, but for fluid at 3000 it allows avoiding the D2 doublet near 2900 cm-1
cutfq = np.array([0.010,0.015])
sampleultimate = spectrafilter(sampleultimate,'bandstop',cutfq,1,np.array([1]))
peakOD = sampleultimate[np.where((ese[:,0]> lb) & (sampleultimate[:,0] < mb1))]
peakOH = sampleultimate[np.where((sampleultimate[:,0]> mb2) & (sampleultimate[:,0] < hb))]
esepeakOD = ese[np.where((ese[:,0]> lb) & (ese[:,0] < mb1))]
esepeakOH = ese[np.where((ese[:,0]> mb2) & (ese[:,0] < hb))]
AireOD = np.trapz(abs(peakOD[:,1]),peakOD[:,0])
AireOH = np.trapz(peakOH[:,1],peakOH[:,0])
eseAireOD = np.trapz(esepeakOD[:,1],esepeakOD[:,0])
eseAireOH = np.trapz(esepeakOH[:,1],esepeakOH[:,0])
ratioDH = AireOD/AireOH
eseratioDH = np.sqrt((1/AireOH)**2*eseAireOD**2+((AireOD-AireOH)/(AireOH**2))**2*eseAireOH**2)


figure(1,figsize=(10,6))
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

ax3.plot(sampleultimate[:,0],sampleultimate[:,1]/amax(peakOD[:,1])*100 ,'m-')

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

# data are out in this directory
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
savefilename = asksaveasfilename() # show an "Open" dialog box and return the path to the selected file
np.savetxt(savefilename,sampleultimate)

namesample = savefilename[0:savefilename.find('.')]
namefig = namesample+'.pdf'
savefig(namefig) # save the figure
 
