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
sys.path.append("/Users/charleslelosq/anaconda/lib/python2.7/lib-charles")# my library


import csv
import numpy
import scipy
import matplotlib
import matplotlib.gridspec as gridspec
from pylab import *
from StringIO import StringIO
from scipy import interpolate

from spectratools import * 

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfile
import os

def removebas_melt(name):
    
    ##### COMMON X VALUE
    x = np.arange(280,1250,0.1)
    sample = np.zeros((len(x),3))
    rawsample = np.genfromtxt(name,skip_header=20, skip_footer=43) # WARNING adjust header and footer lines if needed
    tck = interpolate.splrep(rawsample[:,0],rawsample[:,1],s=0)
    y = interpolate.splev(x,tck,der=0)
    sample[:,0] = x
    sample[:,1] = y 
    sample[:,2] = np.sqrt(np.abs(y))
    
    #########################################################
    #### BASELINE INTERPOLATION REGIONS For melt
      
    # Automatic search  is kind of dangerous,
    # because if the background raises fastly you will just found the lower seach bond...
    # but if it works you can use this piece of code to do it:
    bebe = sample[np.where((sample[:,0] > 680)&(sample[:,0]<1000))] # if necessary adjust bonds
    minsp1 = bebe[np.where((bebe[:,1] == np.min(bebe[:,1])))]
        
    bebe4 = sample[np.where((sample[:,0] > 1100)&(sample[:,0]<1300))]
    minsp4 = bebe4[np.where((bebe4[:,1] == np.min(bebe4[:,1])))]    
    
    # For the melt we do not constrain to the beginning of spectra near 280 cm-1
    # only befire the Q1 and before the diamond peaks.
    bir = np.array([(680,730),(minsp4[:,0],1230)]) # if necessary adjust bonds manualy
    corrsample, baselineD, coeffsD = linbaseline(sample,bir,'spline',10000000) # if necessary adjust spline coeff
    
    # Plot the things to have a look at it
    figure()
    plot(x,sample[:,1],'k-')
    plot(x,baselineD[:,1],'r-')
    plot(x,corrsample[:,1]+(max(sample[:,1])/2),'g-') # corrected sample is shifted for better view
    xlim(280,1250)
    ylim(0,max(sample[:,1])) # automatic limit
    
    corrsample[:,2] = sample[:,2]/sample[:,1]*corrsample[:,1] # ese as dy/y * ycorr
    
    return corrsample

def removebas_fluid(name):

    # same comments as above    
    
    ##### COMMON X VALUE
    x = np.arange(280,1250,0.1)
    sample = np.zeros((len(x),3))
    rawsample = np.genfromtxt(name,skip_header=20, skip_footer=43)
    tck = interpolate.splrep(rawsample[:,0],rawsample[:,1],s=0)
    y = interpolate.splev(x,tck,der=0)
    sample[:,0] = x
    sample[:,1] = y 
    sample[:,2] = np.sqrt(np.abs(y))
     
    #### BASELINE INTERPOLATION REGIONS For fluid
    # same comment as above for the automatic seach
    bebe = sample[np.where((sample[:,0] > 450)&(sample[:,0]<470))] # if necessary adjust bonds
    minsp1 = bebe[np.where((bebe[:,1] == np.min(bebe[:,1])))]
    bebe2 = sample[np.where((sample[:,0] > 680)&(sample[:,0]<720))] # if necessary adjust bonds
    minsp2 = bebe2[np.where((bebe2[:,1] == np.min(bebe2[:,1])))]
    bebe3 = sample[np.where((sample[:,0] > 960)&(sample[:,0]<985))]
    minsp3 = bebe3[np.where((bebe3[:,1] == np.min(bebe3[:,1])))]
    bebe4 = sample[np.where((sample[:,0] > 1100)&(sample[:,0]<1300))]
    minsp4 = bebe4[np.where((bebe4[:,1] == np.min(bebe4[:,1])))]
    
    bir = np.array([(295,300),(450,460),(minsp3[:,0]-3,minsp3[:,0]+3),(1120,1250)]) # if necessary adjust bonds
    corrsample, baselineD, coeffsD = linbaseline(sample,bir,'spline',20000000) # if necessary adjust spline coeff
  
    figure()
    plot(x,sample[:,1],'k-')
    plot(x,baselineD[:,1],'r-')
    plot(x,corrsample[:,1]+(max(sample[:,1])/2),'g-')
    xlim(200,1250)
    ylim(0,max(sample[:,1]))
    
    corrsample[:,2] = sample[:,2]/sample[:,1]*corrsample[:,1] 
    
    return corrsample


##### NOW WE USE THE FUNCTIONS TO TREAT OUR DATA

tkMessageBox.showinfo(
            "Open file",
            "Please open the list of spectra")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
with open(filename) as inputfile:
    results = list(csv.reader(inputfile)) # we read the data list

#### WARNING WE ASSUME YOU HAVE A ./treated/ folder to record treated spectra !!!!!!
for lg in range(len(results)): # we do a loop of this data list
    name = str(results[lg]).strip('[]') # remove unwanted [] in the name of data file
    name = name[1:-1] # to remove unwanted ""
    out = removebas_fluid(name)
    name.rfind('/')
    nameout = name[name.rfind('/')+1::]   
    pathbeg = filename[0:filename.rfind('/')]
    pathint = str('/treated/') # WARNING OUTPUT FOLDER MUST EXIST
    pathout = pathbeg+pathint+nameout
    np.savetxt(pathout,out) # save the spectra
    # for figures
    nameoutfig = nameout[:nameout.rfind('.')]
    outfigpath = pathbeg+pathint+nameoutfig+'.pdf'
    savefig(outfigpath)
    
  
