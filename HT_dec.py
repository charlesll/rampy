# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 16:26:34 2014

Modfied on Nov. 2014 for spectral differentiation and testing how
hydrogen bonding affects DH fractionation

Input => one line = paths of fluid and melts spectra with their 1000 ln(alpha) values and the initial DH values

@author: charleslelosq
"""
import numpy as np
import scipy
import matplotlib
import matplotlib.gridspec as gridspec
from pylab import *

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfile


# Collecting the list of spectra
tkMessageBox.showinfo(
            "Open ",
            "Please open the list of spectra")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
samplename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

# we import the information in an array, skipping the first line
dataliste = np.genfromtxt(samplename,dtype = 'string', skip_header=0,skip_footer=0)

pathfluid = (dataliste[:,0])
pathmelt = (dataliste[:,1])
alphas = np.genfromtxt(dataliste[:,2])
esealphas = np.genfromtxt(dataliste[:,3])
DHinit = dataliste[:,4]

xOH = np.arange(2200,3900,0.2)
ratio = np.zeros((len(alphas),9)) # Total aire diff sp, ese, positiv aire diff sp, ese

for i in range(len(alphas)): # We loop over in dataliste
    rawspfluid = np.genfromtxt(pathfluid[i]) # full path is expected
    rawspmelt = np.genfromtxt(pathmelt[i]) # full path is expected
    
    # resampling for having the same X axis
    rawfluid = np.zeros((len(xOH),3))
    rawfluid[:,1] = np.interp(xOH,rawspfluid[:,0],rawspfluid[:,1])
    rawfluid[:,0] = xOH  
    rawfluid[:,2] = sqrt(rawfluid[:,1])/rawfluid[:,1] #relative error  
    rawmelt = np.zeros((len(xOH),3))
    rawmelt[:,1] = np.interp(xOH,rawspmelt[:,0],rawspmelt[:,1])
    rawmelt[:,0] = xOH 
    rawmelt[:,2] = sqrt(rawmelt[:,1])/rawmelt[:,1] #relative error       
    
    # Boundaries for the OD and OH stretch peaks
    lbOH = 2900
    hbOH = 3800    
    lbOD = 2200
    hbOD = lbOH
    
    ODfluid = rawfluid[np.where((rawfluid[:,0]>lbOD) & (rawfluid[:,0] < hbOD))]
    ODmelt = rawmelt[np.where((rawmelt[:,0]>lbOD) & (rawmelt[:,0] < hbOD))]
    
    OHfluid = rawfluid[np.where((rawfluid[:,0]>lbOH) & (rawfluid[:,0] < hbOH))]
    OHmelt = rawmelt[np.where((rawmelt[:,0]>lbOH) & (rawmelt[:,0] < hbOH))]
      
    # Normalization to total area  
    aOHfluid = np.trapz(OHfluid[:,1],OHfluid[:,0])
    eseaOHfluid = sqrt(aOHfluid)
    aOHmelt = np.trapz(OHmelt[:,1],OHmelt[:,0])
    eseaOHmelt = sqrt(aOHmelt)
    OHfluid[:,1] = OHfluid[:,1]/aOHfluid
    OHmelt[:,1] = OHmelt[:,1]/aOHmelt
    
    aODfluid = np.trapz(ODfluid[:,1],ODfluid[:,0])
    eseaODfluid = sqrt(aODfluid)
    aODmelt = np.trapz(ODmelt[:,1],ODmelt[:,0])
    eseaODmelt = sqrt(aODmelt)
    ODfluid[:,1] = ODfluid[:,1]/aODfluid
    ODmelt[:,1] = ODmelt[:,1]/aODmelt
    
    diffOH = np.zeros(shape(OHfluid))
    diffOH[:,0] = OHmelt[:,0]
    diffOH[:,1] = OHmelt[:,1] - OHfluid[:,1]
    diffOH[:,2] = np.sqrt(OHmelt[:,2]**2 + OHfluid[:,2]**2)
    
    diffOD = np.zeros(shape(ODfluid))
    diffOD[:,0] = ODmelt[:,0]
    diffOD[:,1] = ODmelt[:,1] - ODfluid[:,1]
    diffOD[:,2] = np.sqrt(ODmelt[:,2]**2+ODfluid[:,2]**2)

    # Here we only quantify the total difference between the OH stretch of the melt and the fluid     
    ratio[i,0] = np.trapz(abs(diffOH[:,1]),diffOH[:,0])
    ratio[i,1] = np.trapz(abs(diffOH[:,2]),diffOH[:,0])
    
    # And here we quantify the difference between the mean O-D and O-H streching frequencies in the melt and the fluid
    # Mean frequencies of vibration
    OHFreqMelt = np.sum(OHmelt[:,0]*(OHmelt[:,1]/np.sum(OHmelt[:,1])))
    OHFreqFluid = np.sum(OHfluid[:,0]*(OHfluid[:,1]/np.sum(OHfluid[:,1])))
    ODFreqMelt = np.sum(ODmelt[:,0]*(ODmelt[:,1]/np.sum(ODmelt[:,1])))
    ODFreqFluid = np.sum(ODfluid[:,0]*(ODfluid[:,1]/np.sum(ODfluid[:,1])))

    pseudoDZPEmelt = OHFreqMelt - ODFreqMelt
    pseudoDZPEfluid = OHFreqFluid - ODFreqFluid
    diffPseudoDPE = pseudoDZPEmelt - pseudoDZPEfluid
    
    figure(i,figsize=(12,6))   
    gs = gridspec.GridSpec(1, 2)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
       
    ax1.plot(ODfluid[:,0],ODfluid[:,1],'b-')
    ax1.plot(ODmelt[:,0],ODmelt[:,1],'r-')
    ax1.plot(diffOD[:,0],diffOD[:,1],'g-')
    ax1.plot(np.array([ODFreqMelt,ODFreqMelt]),np.array([0,1]),'r-')    
    ax1.plot(np.array([ODFreqFluid,ODFreqFluid]),np.array([0,1]),'b-')    
    ax1.set_ylim(((1.10*np.min(diffOD[:,1])),(1.10*np.max(ODfluid[:,1]))))
    
    ax2.plot(OHfluid[:,0],OHfluid[:,1],'b-')
    ax2.plot(OHmelt[:,0],OHmelt[:,1],'r-')
    ax2.plot(diffOH[:,0],diffOH[:,1],'g-')
    ax2.plot(np.array([OHFreqMelt,OHFreqMelt]),np.array([0,1]),'r-')    
    ax2.plot(np.array([OHFreqFluid,OHFreqFluid]),np.array([0,1]),'b-')    
    ax2.set_ylim((1.10*min(diffOH[:,1]),1.10*max(OHfluid[:,1])))
    
    ratio[i,2] = OHFreqMelt    
    ratio[i,3] = OHFreqFluid
    ratio[i,4] = ODFreqMelt
    ratio[i,5] = ODFreqFluid
    ratio[i,6] = pseudoDZPEmelt
    ratio[i,7] = pseudoDZPEfluid
    ratio[i,8] = diffPseudoDPE    
     
figure()
plot(alphas,ratio[:,4],'ro')
 
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
savefilename = asksaveasfile() # show an "Open" dialog box and return the path to the selected file
out = vstack((alphas,esealphas,ratio[:,0],ratio[:,1],ratio[:,2],ratio[:,3],ratio[:,4],ratio[:,5],ratio[:,6],ratio[:,7],ratio[:,8])).T # Output matrix with 4 columns: Alphas, errors, Ratio Hydro bonds, errors
np.savetxt(savefilename,out)