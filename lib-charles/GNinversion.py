# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 11:18:06 2014
NOT GOOD FOR NOW.... PROBLEM WIth ARRAY AND MATRIX tHAT ARE NOT GOOD....
@author: charleslelosq
Functions for making GaussNewton inversion of spectra, with fitting Gaussian bands
"""
import numpy as np
import scipy # faster than numpy for linalg
from numpy import linalg
from numpy import matrix

def derivGauss(params,x,variance):
    nbpic = int(len(params)/3)
    A = np.zeros((nbpic))
    F = np.zeros((nbpic))
    L = np.zeros((nbpic))
    for n in range(nbpic):
        m = 2*n # little trick for correct indexation
        A[n] = params[n+m]
        F[n] = params[n+m+1]
        L[n] = params[n+m+2]       
    
    pic = np.zeros((len(x),nbpic))
    for i in range(nbpic):
        pic[:,i] = A[i] * np.exp(-np.log(2)*((x[:]-F[i])/L[i])**2)
    
    ycalc = np.sum(pic,1)

    dA = np.zeros((len(x),nbpic))
    dAsig = np.zeros((len(x),nbpic))
    dF = np.zeros((len(x),nbpic))
    dFsig = np.zeros((len(x),nbpic))
    dL = np.zeros((len(x),nbpic))
    dLsig = np.zeros((len(x),nbpic))
     
    for m in range(nbpic):
            dA[:,m] = np.exp(-(-F[m] + x)**2*np.log(2)/L[m]**2)
            dAsig[:,m] = dA[:,m]/variance
            dF[:,m] = -A[m]*(2*F[m] - 2*x)*np.exp(-(-F[m] + x)**2*np.log(2)/L[m]**2)*np.log(2)/L[m]**2
            dFsig[:,m] = dF[:,m]/variance
            dL[:,m] = 2*A[m]*(-F[m] + x)**2*np.exp(-(-F[m] + x)**2*np.log(2)/L[m]**2)*np.log(2)/L[m]**3
            dLsig[:,m] = dL[:,m]/variance
            
#    G = np.zeros((len(x),3*nbpic))     
#    GT = np.zeros((len(x),3*nbpic)) 
    G = np.concatenate((dA,dF,dL),1)
    GT = np.concatenate((dAsig,dFsig,dLsig),1)
    GT = GT.T
     
    return ycalc, G, GT


def PGNLQ(params,sigparams,x,y,sig,numberiter,relax):
# ######################## Newton INVERSION ############################

# Pseudo Newton algorhythm, maybe
# mathematically the more rigourous according to Tarantolla.
# See the Tarantola book of 2005, Inverse Problem Theory and Methods for
# Model Parameter Estimation, page69. You can download it for free if you tip that
# into Google Scholar. Don't worry, Albert Tarantola gave it without asking
# for money and the book is free. And the following algorithm is the one
# that I used during the Albert Tarantola's lessons. I just modify it because of the high number of data in Raman spectra sometimes which could overflow the computer buffers...
# So we do not have the intervention of the CD and ICD matrix into the
# algorithm, because of the division of GT matrix by the data errors. 

    # we use matrix class of numpy
    params = matrix(params)
    x = matrix(x)
    y = matrix(y)
    sig = matrix(sig)

    CM = matrix.diagonal(sigparams**2) # Here its the matrix of the prior model uncertainty...
    ICM = CM.I # and its inverse
                
    var = sig**2 # calculating the variance as the square of 1sig errors  
    # we do not need CD and ICD if we devideGT by the variance in the called function, we do that because too much numbers if we do the direct calculation with matrix...
    #CD = np.diag(var)    
    #CD = matrix(CD)
    #ICD = CD.I

    mprior = params # single vector of parameters
    mcurrent = mprior # the algorithm is initialized at the a priori model values
    
    #boucle
    iteration = 1 # Just one variable to see the number of iteration 
    #dchi2 = 1; # One parameter for the chi2 calculation
    
    while iteration < numberiter: # So iterations will continue untill the difference of the chi2 is more than 1e-10. You can adjust that to your problem too
                   
        ycalc,G,GT = derivGauss(mcurrent,x,var) # Here it is the more important part. I call the derivGauss function which calculate the derivatives        
        dataresiduals = ycalc-y # as the name variable says...
        modelresiduals = mcurrent-mprior # same, note that at the first iteration mcurrent = mprior
        #chi2 = sum((y-ycalc)**2/var) # The "old" chi2
                      
        gradient = GT*dataresiduals + ICM*modelresiduals
        d1 = GT*G+ICM
        direction = d1.I*gradient# The direction 
        
        # Here we "relax" just the direction, in order to avoid the model to
        # converge too fast and sometime explode. Note that you can change the
        # 2 by another values such as 1 if you have a well constrained
        # problem with only one solution (one local minimal), but maybe 4 or 6 if you have
        # many local minima and so "intermediate" solutions:
        mnew = mcurrent - direction/relax
        
        #sfunction = sum((mnew-mcurrent)**2./mcurrent) # The misfit function as defined into Tarantola, 2005. just in case...
                     
        mcurrent = mnew # Here I set the current model to be equal to the new model
    
        #ycalcnew, G, GT = derivGauss(params,x,sig) # Pour calcul du chi2
        #chi2new = np.sum((y-ycalc)**2/var) # The "new" chi2
        #dchi2 = chi2-chi2new;
        #dchi2=dchi2+1;
        iteration = iteration + 1;
        
    return mcurrent  