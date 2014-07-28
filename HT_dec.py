# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 16:26:34 2014

@author: charleslelosq
"""
import numpy as np
import scipy

from scipy.optimize import curve_fit
from scipy.optimize import fmin_bfgs
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import fmin as simplex
from scipy.optimize import fmin_powell as powell

from Tkinter import *
import tkMessageBox
from tkFileDialog import askopenfilename
from tkFileDialog import asksaveasfile


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
    
def derivGauss(params,x, y, e):
    nbpic = int(len(params)/3)
    fprime = np.zeros(len(params)) #trois parametres * n gaussians    
    
    for n in range(nbpic):
        m = 2*n # little trick for correct indexation
        a1 = params[n+m]
        b1 = params[n+m+1]
        c1 = params[n+m+2]
        
        fprime[n+m]= sum(2*(a1*exp(-(-b1 + x)**2*log(2)/(c1**2)) - y)*exp(-(-b1 + x)**2*log(2)/(c1**2))/(e**2))
        fprime[n+1+m]= sum(-2*a1*(2*b1 - 2*x)*(a1*exp(-(-b1 + x)**2*log(2)/(c1**2)) - y)*exp(-(-b1 + x)**2*log(2)/(c1**2))*log(2)/(c1**2*e**2))
        fprime[n+2+m]= sum(4*a1*(-b1 + x)**2*(a1*exp(-(-b1 + x)**2*log(2)/(c1**2)) - y)*exp(-(-b1 + x)**2*log(2)/(c1**2))*log(2)/(c1**3*e**2))
    
    #fprime = np.concatenate((grad1,grad2,grad3),1)    
    
    return fprime
    
def gauss_lsq(params,x): #attention il y a aussi un multigaussian dans spectra tools, attention aux noms ici
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
    
errfunc = lambda p, x, y: gauss_lsq(p, x) - y # Distance to the target function

# This software is intended to do online fitting
# So here is for try on one spectra
# But will implement a way to treat a list of spectra soon
# Collecting the sample spectrum
tkMessageBox.showinfo(
            "Open file",
            "Please open the sample spectrum")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
sample = np.genfromtxt(filename,skip_header=20, skip_footer=43) # Skipping lines from JASCO files


# If you want the following fit the OD peak
#Fitting the OD peak with the multigaussian_local, DO NOT FORGET TO CHANGE THE FUNCTION IF YOU CHANGE NUMBER OF PEAKS
#params_OD_out, cov_OD_out = curve_fit(multigaussian_local, peakOD[:,0], peakOD[:,1],params_in)
peakOD = sample[np.where((sample[:,0]> 2200) & (sample[:,0] < 2900))]
ese = np.sqrt(np.abs(peakOD[:,1]))#/peakOD[:,1] #taken as a positive percentage
peakOD[:,1] = peakOD[:,1]/amax(peakOD[:,1])*100 
#ese = ese#*peakOD[:,1]

#For inversion
x = peakOD[:,0]
d = peakOD[:,1]

x0 = np.array((14.10698991e+00,   2509,   3.07119981e+01,
         60,   2.59747641e+03,   35,
         35,   2.64934712e+03,   3.07536403e+01))
xopt = fmin_bfgs(func, x0, args=(x, d, ese),fprime=derivGauss,gtol = 1e-15) #approx_grad=1
#xopt, cov_xopt = scipy.optimize.leastsq(errfunc,x0,args = (x,y))
#xopt = powell(func, x0, args=(peakOD[:,0], peakOD[:,1], ese))
ytot, y= multigaussian_local(xopt,x)


###### If not anted comment the following section
##### Now we calculate the relationship between vOH and vOD and O-O distances
#tkMessageBox.showinfo(
#            "Open file",
#            "Please open the OO vs OD/OH distance file")
#Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
#filename3 = askopenfilename() # show an "Open" dialog box and return the path to the selected file
#distances = np.genfromtxt(filename3,skip_header=1)
#vOD = distances[:,0]/distances[:,1]
#oo = distances[:,2]
#spline_vDH = UnivariateSpline(vOD,oo, s=1) 

figure()

plot(peakOD[:,0],peakOD[:,1],'k-')
plot(peakOD[:,0],ytot,'r-')
plot(peakOD[:,0],y[:,0],'m-')
plot(peakOD[:,0],y[:,1],'m-')
plot(peakOD[:,0],y[:,2],'m-')

#plot(peakOD[:,0],y[:,3],'m-')
