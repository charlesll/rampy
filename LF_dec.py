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
import numpy
import scipy
import matplotlib
import matplotlib.gridspec as gridspec
from pylab import *
from StringIO import StringIO
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import fmin as simplex
from scipy.optimize import fmin_powell as powell

from spectratools import * 
from GNinversion import PGNLQ

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename
import os

# On est obligé de faire une fonction locale adaptee pour fitter car curve_fit donne une liste de parametre et non une matrice...
def multigaussian_local(params,x): #attention il y a aussi un multigaussian dans spectra tools, attention aux noms ici
    nbpic = int(len(params)/3)
    a = np.zeros((1,nbpic))
    b = np.zeros((1,nbpic))
    c = np.zeros((1,nbpic))
    y = np.zeros((len(x),nbpic))
    for n in range(nbpic):
        m = 2*n # little trick for correct indexation
        a[0,n] = params[n+m]
        b[0,n] = params[n+m+1]
        c[0,n] = params[n+m+2]
        y[:,n] = a[0,n]*np.exp(-np.log(2)*((x[:]-b[0,n])/c[0,n])**2)
    ytot = sum(y,1)
    
    return ytot, y
     
def multigaussian_local_lfix(params,x): #attention il y a aussi un multigaussian dans spectra tools, attention aux noms ici
    nbpic = int((len(params)-2)/2)
    a = np.zeros((1,nbpic))
    b = np.zeros((1,nbpic))
    c = np.zeros((1,2)) # FWMH fixed and in first position of params
    c[0,0] = params[0]
    c[0,1] = params[1]
    b[0,:] = params[2:(nbpic+2)]
    a[0,:] = params[nbpic+2:(2*nbpic+2)]
    y = np.zeros((len(x),nbpic))
    for n in range(nbpic):
        if n == 0:
            y[:,n] = a[0,n]*np.exp(-np.log(2)*((x[:]-b[0,n])/c[0,0])**2)
        else:
            y[:,n] = a[0,n]*np.exp(-np.log(2)*((x[:]-b[0,n])/c[0,1])**2)
    ytot = sum(y,1)
    
    return ytot, y     
     
def func(params, X, Y, Err):
    nbpic = int(len(params)/3)
    a = np.zeros((1,nbpic))
    b = np.zeros((1,nbpic))
    c = np.zeros((1,nbpic))
    for n in range(nbpic):
        m = 2*n # little trick for correct indexation
        a[0,n] = params[n+m]
        b[0,n] = params[n+m+1]
        c[0,n] = params[n+m+2]    
    
    # compute chi-square
    chi2 = 0.0
    for n in range(len(X)):
        x = X[n]
        y = a *np.exp(-np.log(2)*((x-b)/c)**2)
        y = sum(y,1)  
        chi2 = chi2 + (Y[n] - y)*(Y[n] - y)/(Err[n]*Err[n])
    return chi2
        
#def derivGauss(params,x, y, e):
#    nbpic = int(len(params)/3)
#    fprime = np.zeros(len(params)) #trois parametres * n gaussians    
#    
#    for n in range(nbpic):
#        m = 2*n # little trick for correct indexation
#        a1 = params[n+m]
#        b1 = params[n+m+1]
#        c1 = params[n+m+2]
#        
#        fprime[n+m]= sum(2*(a1*exp(-(-b1 + x)**2*log(2)/(c1**2)) - y)*exp(-(-b1 + x)**2*log(2)/(c1**2))/(e**2))
#        fprime[n+1+m]= sum(-2*a1*(2*b1 - 2*x)*(a1*exp(-(-b1 + x)**2*log(2)/(c1**2)) - y)*exp(-(-b1 + x)**2*log(2)/(c1**2))*log(2)/(c1**2*e**2))
#        fprime[n+2+m]= sum(4*a1*(-b1 + x)**2*(a1*exp(-(-b1 + x)**2*log(2)/(c1**2)) - y)*exp(-(-b1 + x)**2*log(2)/(c1**2))*log(2)/(c1**3*e**2))
    

    #fprime = np.concatenate((grad1,grad2,grad3),1)    
#    
#    return fprime
    
def gauss_lsq(params,x): 
    nbpic = int(len(params)/3)
    a = np.zeros((1,nbpic))
    b = np.zeros((1,nbpic))
    c = np.zeros((1,nbpic))
    y = np.zeros((len(x),nbpic))
    for n in range(nbpic):
        m = 2*n # little trick for correct indexation
        a[0,n] = params[n+m]
        b[0,n] = params[n+m+1]
        c[0,n] = params[n+m+2]
        y[:,n] = a[0,n]*np.exp(-np.log(2)*((x[:]-b[0,n])/c[0,n])**2)
    ytot = sum(y,1)
    
    return ytot
 
def gauss_lsq_lfix(params,x):
    nbpic = int((len(params)-2)/2)
    a = np.zeros((1,nbpic))
    b = np.zeros((1,nbpic))
    c = np.zeros((1,2)) # FWMH fixed and in first position of params
    c[0,0] = params[0]
    c[0,1] = params[1]
    b[0,:] = params[2:(nbpic+2)]
    a[0,:] = params[nbpic+2:(2*nbpic+2)]
    y = np.zeros((len(x),nbpic))
    for n in range(nbpic):
        if n == 0:
            y[:,n] = a[0,n]*np.exp(-np.log(2)*((x[:]-b[0,n])/c[0,0])**2)
        else:
            y[:,n] = a[0,n]*np.exp(-np.log(2)*((x[:]-b[0,n])/c[0,1])**2)
    ytot = sum(y,1)
    
    return ytot   
 
###### ERROR FUNCTION FOR THE FIT      
errfunc = lambda p, x, y: gauss_lsq(p, x) - y # Distance to the target function that is chi2

def deconvolution(name):
    # get the sample to deconvolute
    sample = np.genfromtxt(name)

    ## Fit of spectra
    # we search for peak minimum near 700 cm-1 for melt (work in our case)
    bebe = sample[np.where((sample[:,0] > 650)&(sample[:,0]<750))]
    minsp1 = bebe[np.where((bebe[:,1] == np.min(bebe[:,1])))]
    minsp1[:,0] = 720 # for fluid

    interestspectra = sample[np.where((sample[:,0] > minsp1[:,0])&(sample[:,0] < 1250))]
    ese0 = interestspectra[:,2]/abs(interestspectra[:,1]) #take ese  as a percentage
    interestspectra[:,1] = interestspectra[:,1]/np.amax(interestspectra[:,1])*100 # normalise spectra to maximum, easier to handle after 
    sigma = np.zeros((len(interestspectra[:,0])))+1#np.abs(ese0*interestspectra[:,1]) #calculate good ese
    #sigma = None

    xfit = interestspectra[:,0]
    d = interestspectra[:,1]

    #Initial guess.
    x0 = np.array([(21,763,19,20,830,35,27,890,45,30,1000,41,83,1065,32)]) # for melt lfree
    cons = [(0,50),(740,790),(0,60),(0,50),(800,850),(20,50),(0,50),(850,930),(20,50),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None)]   
    #x0 = np.array([(21,558,19,20,641,35,100,760,27,30,820,41,83,890,32,10,1053,30)]) # for fluid lfree
    #x0 = np.array([(19,40,763,832,891,1000,1065,30,40,50,60,90)]) # for melt lfix
    #x0 = np.array([(20,40,760,558,641,820,890,1053,100,20,20,20,20,20)]) # for fluid lfix
    # Apply algorythm.
    xopt, f, dic = fmin_l_bfgs_b(func, x0, args=(xfit, d, sigma),approx_grad=1,bounds = cons, epsilon=1e-05) #approx_grad=1
    #xopt, covar = scipy.optimize.leastsq(errfunc,x0,args = (xfit,d))
    


    # Gauss Newton Algorythm
    # Not good yet.... Problems ith matrix and arrays....
#    params = np.array((33, 761,17,17,815,31,43,871,45,41,1015,47,77,1065,37)) #amplitude, frequency, FWMH
#    sigparams = np.array((0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
#    xopt = PGNLQ(params,sigparams,xfit,d,sigma,100,2)  
 

    ycalc, pics = multigaussian_local(xopt,xfit)
    out = np.column_stack((xfit,d,ycalc,pics))     
    
    
#    ####################################################
#    ## 3. COMPUTE THE FIT AND FIT ERRORS USING bootstrap
#    ####################################################        
#
#    # Estimate the confidence interval of the fitted parameter using
#    # the bootstrap Monte-Carlo method
#    # http://phe.rockefeller.edu/LogletLab/whitepaper/node17.html
#    residuals = errfunc( xopt, xfit, d)
#    s_res = numpy.std(residuals)
#    ps = []
#    # 100 random data sets are generated and fitted
#    for i in range(300):
#      randomDelta = numpy.random.normal(0., s_res, len(d))
#      randomdataY = d + randomDelta   
#      randomfit, randomcov = \
#      scipy.optimize.leastsq( errfunc, x0, args=(xfit, randomdataY),\
#                        full_output=0)
#      ps.append( randomfit ) 
#
#    ps = numpy.array(ps)
#    mean_pfit = numpy.mean(ps,0)
#    Nsigma = 1. # 1sigma gets approximately the same as methods above
#                # 1sigma corresponds to 68.3% confidence interval
#                # 2sigma corresponds to 95.44% confidence interval
#    err_pfit = Nsigma * numpy.std(ps,0) 

    fig = figure()
    plot(interestspectra[:,0],interestspectra[:,1],'k-')
    plot(xfit,ycalc,'r-')
    for i in range(int(len(xopt))/3):
         plot(xfit,pics[:,i],'b-')
    xlim(700,1230)
    ylim(-5,105)
    xlabel("Raman shift, cm$^{-1}$", fontsize = 18, fontweight = "bold")
    ylabel("Normalized intensity, a. u.", fontsize = 18, fontweight = "bold")
    text(650,110,name,fontsize = 18,color='r')     
    return out, xopt, fig




tkMessageBox.showinfo(
            "Open file",
            "Please open the list of spectra")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
with open(filename) as inputfile:
    results = list(csv.reader(inputfile)) # we read the data list

for lg in range(4):
    name = str(results[lg]).strip('[]')
    name = name[1:-1] # to remove unwanted ""
    # get the sample to deconvolute
    sample = np.genfromtxt(name)

    ## Fit of spectra
    # we search for peak minimum near 700 cm-1 for melt (work in our case)
    bebe = sample[np.where((sample[:,0] > 650)&(sample[:,0]<750))]
    minsp1 = bebe[np.where((bebe[:,1] == np.min(bebe[:,1])))]
    minsp1[:,0] = 720 # for fluid

    interestspectra = sample[np.where((sample[:,0] > minsp1[:,0])&(sample[:,0] < 1250))]
    ese0 = interestspectra[:,2]/abs(interestspectra[:,1]) #take ese  as a percentage
    interestspectra[:,1] = interestspectra[:,1]/np.amax(interestspectra[:,1])*100 # normalise spectra to maximum, easier to handle after 
    sigma = np.zeros((len(interestspectra[:,0])))+1#np.abs(ese0*interestspectra[:,1]) #calculate good ese
    #sigma = None

    xfit = interestspectra[:,0]
    d = interestspectra[:,1]

    #Initial guess.
    x0 = np.array([(21,763,19,20,830,35,27,890,45,30,1000,41,83,1065,32)]) # for melt lfree
    cons = [(0,50),(740,790),(0,60),(0,50),(800,850),(20,50),(0,50),(850,930),(20,50),(None,None),(None,None),(None,None),(None,None),(None,None),(None,None)]   
    #x0 = np.array([(21,558,19,20,641,35,100,760,27,30,820,41,83,890,32,10,1053,30)]) # for fluid lfree
    #x0 = np.array([(19,40,763,832,891,1000,1065,30,40,50,60,90)]) # for melt lfix
    #x0 = np.array([(20,40,760,558,641,820,890,1053,100,20,20,20,20,20)]) # for fluid lfix
    # Apply algorythm.
    xopt, f, dic = fmin_l_bfgs_b(func, x0, args=(xfit, d, sigma),approx_grad=1,bounds = cons, epsilon=1e-05) #approx_grad=1
    ycalc, pics = multigaussian_local(xopt,xfit)
    fig = figure()
    plot(interestspectra[:,0],interestspectra[:,1],'k-')
    plot(xfit,ycalc,'r-')
    for i in range(int(len(xopt))/3):
         plot(xfit,pics[:,i],'b-')
    xlim(700,1230)
    ylim(-5,105)
    xlabel("Raman shift, cm$^{-1}$", fontsize = 18, fontweight = "bold")
    ylabel("Normalized intensity, a. u.", fontsize = 18, fontweight = "bold")
    text(650,110,name,fontsize = 18,color='r')     

    #yout, model, fig = deconvolution(name)
#    name.rfind('/')
#    nameout = name[name.rfind('/')+1::]
#    namesample = nameout[0:nameout.find('.')]
#    pathbeg = filename[0:filename.rfind('/')]
#    pathint = str('/deconv/')
#    ext1 = '_ydec.txt'
#    ext2 = '_params.txt'
#    #ext3 = '_paramserrors.txt'
#    ext4 = '.pdf'
#    pathout1 = pathbeg+pathint+namesample+ext1
#    pathout2 = pathbeg+pathint+namesample+ext2
#    #pathout3 = pathbeg+pathint+namesample+ext3
#    pathout4 = pathbeg+pathint+namesample+ext4
#    np.savetxt(pathout1,yout)
#    np.savetxt(pathout2,model)
#    #np.savetxt(pathout3,errors)
#    savefig(pathout4)
    
#### IF YOU WANT TO ONLY TREAT ON FILE? PLEASE UNCOMMENT THE FOLLOWING AND COMMENT BEFORE
## data are out in this directory
#Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
#filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
#out = removebas(filename)
#Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
#savefilename = asksaveasfile() # show an "Open" dialog box and return the path to the selected file
#np.savetxt(savefilename,out)