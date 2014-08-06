# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 17:52 2014

@author: charleslelosq
Carnegie Institution for Science

This script calculates pressure (MPa) in DAC experiments from Raman shift of 13C diamonds
Put it anywhere, you need however the lib-charles library as well as numpy, scipy, matplotlib, and Tkinter
that usually come with any python distribution
"""
import sys
sys.path.append("/Users/Celita/Desktop/RamPy-master/lib-charles/")
sys.path.append("/Users/Celita/Desktop/RamPy-master/lib-charles/gcvspl/")


import numpy as np
import scipy
import matplotlib
import matplotlib.gridspec as gridspec
from pylab import *

from spectratools import *
from IRTABS import *
import gcvspline

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfile
     
#### First thing: subtracting the second order of diamond from the experimental spectra

# Collecting the sample spectrum
tkMessageBox.showinfo(
            "Open file",
            "Please open the I spectrum")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
samplename = askopenfilename() # show an "Open" dialog box and return the path to the selected file


# Same for diamond
tkMessageBox.showinfo(
            "Open file",
            "Please open the Io spectrum")
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
backgroundname = askopenfilename() # show an "Open" dialog box and return the path to the selected file

# data are out in this directory
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
savefilename = asksaveasfile() # show an "Open" dialog box and return the path to the selected file

thick = 50*10**-6 # thickness in m

#### Now we use IRTABS
output = FTIRout(samplename,backgroundname,savefilename,thick,header_lines = 19,footer_lines = 37)

# Now we can subtract a 6-order polynomial baseline below the peaks in the region of interest
# Here the lower and upper bounds of this regions:
lb = 3900
hb = 6000

output[:,1] = output[:,1]/100 # because people give absorption in cm usually...
interestIR = output[np.where((output[:,0]>lb) & (output[:,0]<hb))] # we select data in the interest region

#### If needed their is here a filter
orisp = interestIR
cutfq = np.array([0.006])
interestIR = spectrafilter(interestIR,'low',cutfq,1,np.array([1]))
#cutfq = np.array([0.010,0.015])
#interestIR = spectrafilter(interestIR,'bandstop',cutfq,1,np.array([1]))

# we estimate errors as
errors= np.sqrt(interestIR[:,1])
dataset = np.zeros((len(interestIR),3))
dataset[:,0] = interestIR[:,0]
dataset[:,1] = interestIR[:,1]
dataset[:,2] = errors

#### BIR values:
b1 = lb
b2 = 4289
b3 = 4761
b4 = 4838
b5 = 5667
b6 = 6000

bir = np.array([(b1,b2),(b3,b4),(b5,b6)]) # BIR for constraining the baseline
corrIR, baseline, coeffs = linbaseline(dataset,bir,'gcvspline',0.004) # # Spline baseline with mode 3 of gcvspl.f
#corrIR, baseline, coeffs = linbaseline(interestIR,bir,'unispline',0.1) # Spline baseline with univariate spline of scipy

       
# Select interest areas for calculating the areas of OH and H2Omol peaks
intarea45 = corrIR[np.where((corrIR[:,0]> b2) & (corrIR[:,0]<b3))]
intarea52 = corrIR[np.where((corrIR[:,0]> b4) & (corrIR[:,0]<b5))]

area4500 = np.trapz(intarea45[:,1],intarea45[:,0])
area5200 = np.trapz(intarea52[:,1],intarea52[:,0])

# We assume that RELATIVE errors on areas are globally equal to 1/sqrt(Area)
esearea4500 = 1/sqrt(area4500)
esearea5200 = 1/sqrt(area5200)

ratioOH_H2O = area4500/area5200
eseratioOH_H2O = (esearea4500 + esearea5200)*ratioOH_H2O

figure(1,[10,10]) # And we do the figure
gs = matplotlib.gridspec.GridSpec(1, 3)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])

ax1.plot(output[:,0],output[:,1],'k-') # The global absorbance spectra
ax2.plot(interestIR[:,0],interestIR[:,1],'b-') # The spectra in the region of interest
ax2.plot(orisp[:,0],orisp[:,1],'k-')
ax2.plot(baseline[:,0],baseline[:,1],'g-') # The polynomial baseline

ax3.plot(corrIR[:,0],corrIR[:,1],'r-') # The spectra without the baseline

# Labels
ax1.set_ylabel("Absorbance, 1/cm", fontsize = 18, fontweight = "bold")
ax2.set_xlabel("Wavelength, cm$^{-1}$", fontsize = 18, fontweight = "bold")

#limits
ax2.set_xlim(lb,hb)
ax3.set_xlim(lb,hb)
ax2.set_ylim(np.min(interestIR[:,1])-0.1*np.min(interestIR[:,1]),np.max(interestIR[:,1])) # Automatic scaling
ax3.set_ylim(-2,np.max(corrIR[:,1])+0.2*np.max(corrIR[:,1]))

# We output the values directly on the plot
ax3.text(4000,np.max(intarea45[:,1])+0.03*np.max(intarea45[:,1]),('Area OH: \n'+'%.1f' % area4500),color='b',fontsize = 16)
ax3.text(4650,np.max(intarea52[:,1])+0.03*np.max(intarea52[:,1]),('Area H$_2$O$_{mol}$: \n'+ '%.1f' % area5200),color='b',fontsize = 16)
ax3.text(lb+100,np.mean((np.max(intarea45[:,1]),np.max(intarea52[:,1]))),('OH/H$_2$O$_{mol}$: \n'+'%.3f' % ratioOH_H2O+'\n+/-'+'%.3f' % eseratioOH_H2O),color='r',fontsize = 16)

nameout = savefilename.name
nameout3 = nameout[0:nameout.rfind('.')]
nameout3_data = (nameout3+'_bassub.txt')
nameout3_graph = (nameout3+'.pdf')
out2 = np.transpose(np.vstack((interestIR[:,0],interestIR[:,1],baseline[:,1],corrIR[:,1])))
np.savetxt(nameout3_data,out2) # Here a matrix containing the x, I/Io, baseline, corrected spectrum
savefig(nameout3_graph) # And we save the graph
