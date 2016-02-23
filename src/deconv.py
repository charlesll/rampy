# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 11:28:55 2015

@author: closq

Library with deconvolution algorithms
Warning: deconvolution, not peak fitting!

Most have been found on the web or in several papers

The spectra have to be normalized to unit area prior to proceessing

"""
import numpy as np
import matplotlib
from pylab import *
import matplotlib.gridspec as gridspec
    
def figureresults(spectrum,deconvoluted, ok_save, reconstructed, mse):
    """
    This is an function for plotting the results of the deconvolution
    
    Spectrum: a 2 columns x y numpy array of length M
    ok1, ok_save, recons and mse are outputs from the deconvolution functions
        
    """    
    plt.figure() 
    gs = gridspec.GridSpec(2, 2)
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])    
    ax3 = plt.subplot(gs[1,0])
    ax4 = plt.subplot(gs[1,1])
    
    ax1.plot(spectrum[:,0],spectrum[:,1],'k-')
    ax1.plot(spectrum[:,0],reconstructed[:],'b-')
    
    ax2.plot(spectrum[:,0],deconvoluted[:,0],'r-')
    
    for i in range(5):
        ax3.plot(spectrum[:,0],spectrum[:,1],'k-',spectrum[:,0],ok_save[:,i],color = (0,i/5.,0,0))
    ax3.plot(spectrum[:,0],deconvoluted[:,0],'r-')
    
    ax4.plot(mse[:,0])
        
    ax1.set_title('Original and reconstruction')
    ax2.set_title('Deconvoluted')
    ax3.set_title('Intermediate spectra')
    ax4.set_title('MSE vs iterations')
    
    plt.tight_layout()
    
         
def Janson(spectrum,s,iterations,early,earlyit,ko,relaxB):
    """
    This is an Python implementation for 1D data of the Janson Van Cittert Algorithm
    as it is given in Janson, 1997, Deconvolution of Images and Spectra
    
    Spectrum: a 2 columns x y numpy array of length M
    s: M*M spread array
    Iterations: integer
    ko: relaxation coefficient, float
    relaxB: maximum bondary, float
    
    The relaxation parameter assumes that you put a restriction to 0. If you don't want, modify the relaxation equation.
    
    Return ok1, ok_save, recons, mse    
    """
    # Definition of the arrays  
    ok = np.zeros((len(spectrum),1))
    ok1 = np.zeros((len(spectrum),1))
    ok_best = np.zeros((len(spectrum),1))
    ok_save = np.zeros((len(spectrum),5)) # this matrix will record the spectra at 20, 40, 60 and 80% of the iteration process
    
    # Something to save the reconstruction error
    mse = np.zeros((iterations,1))
    
    #Something to count for early stopping
    count = 0
    i = 0
    
    # start iteration loop
    while i < iterations and count < earlyit:
        
        if i == 0: #initialisation
            ok[:,0] = spectrum[:,1]
            ok1[:,0] = spectrum[:,1]   
        else:
            ok[:,0] = ok1[:,0]
            
        for n in range(len(spectrum)):
        
            k_on = ko * (1.-2.*np.abs(ok[n,0]-relaxB/2.)/relaxB)  # relaxation parameter
        
            sum1 = np.sum(s[n,0:(n-1)]*ok1[0:(n-1),0]) # matrix indexing, much faster than creating a new loop!
            sum2 = np.sum(s[n,n::]*ok[n::,0])              
            ok1[n,0] = ok[n,0] + k_on/s[n,n] * (spectrum[n,1] - sum1 - sum2) #final equation
        
        mse[i,0] = np.sum((spectrum[:,1] - np.dot(s,ok1)[:])**2) 
        
        if i>0 and mse[i,0] < mse[i-1,0]:
            ok_best[:,0] = ok1[:,0]
            count = 0
            
        if i == int(iterations/5):
            ok_save[:,0] = ok1[:,0]
        elif i == int(iterations/5)*2-1:
            ok_save[:,1] = ok1[:,0]
        elif i == int(iterations/5)*3-1:
            ok_save[:,2] = ok1[:,0]
        elif i == int(iterations/5)*4-1:
            ok_save[:,3] = ok1[:,0]
        elif i == int(iterations/5)*5-1:
            ok_save[:,4] = ok1[:,0]
        
        i = i + 1
        if early == 1: #Early stopping if failing to improve mse for earlyit iterations
            count = count + 1
        
    recons = np.dot(s,ok_best) 
    
    return ok_best, ok_save, recons, mse
  
def richardsonlucy(spectrum,s,iterations):
    """
    This is an Python implementation for 1D data of the Richardson Lucy Algorithm
    as it is given in Janson, 1997, Deconvolution of Images and Spectra
    
    Spectrum: a 1 columns y numpy array of length M
    s: M*M spread array
    Iterations: integer
    
    Return ok1, ok_save, recons, mse    
    """
    # Definition of the arrays  
    ok = np.zeros((len(spectrum),1))
    ok1 = np.zeros((len(spectrum),1))
    c = np.zeros((len(spectrum),1))
    ok_best = np.zeros((len(spectrum),1))
    ok_save = np.zeros((len(spectrum),5)) # this matrix will record the spectra at 20, 40, 60 and 80% of the iteration process
    
    # Something to save the reconstruction error
    mse = np.zeros((iterations,1))

    # start iteration loop
    for i in range(iterations):
        
        if i == 0: #initialisation
            ok[:,0] = spectrum[:,0]
        else:
            ok[:,0] = ok1[:,0]
            
        for n in range(len(spectrum)):
            c[n,0] = np.sum(s[n,:]*ok[:,0])
            ok1[n,0] = ok[n,0]* np.sum(spectrum[n,0]/c[n,0]*s[n,:])         
        
        mse[i,0] = np.sum((spectrum[:,0] - np.dot(s,ok1)[:])**2) 
    
        if i == 0:
            msemin = mse[i,0]
        elif mse[i,0] < msemin:
            ok_best[:,0] = ok1[:,0]
            msemin = mse[i,0]
            
        if i == int(iterations/5):
            ok_save[:,0] = ok1[:,0]
        elif i == int(iterations/5)*2-1:
            ok_save[:,1] = ok1[:,0]
        elif i == int(iterations/5)*3-1:
            ok_save[:,2] = ok1[:,0]
        elif i == int(iterations/5)*4-1:
            ok_save[:,3] = ok1[:,0]
        elif i == int(iterations/5)*5-1:
            ok_save[:,4] = ok1[:,0]
            
    recons = np.dot(s,ok_best) 
    
    return ok_best, ok_save, recons, mse

