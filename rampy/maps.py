#!/usr/bin/env python
#-*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

from rampy import peak_shapes

def read_renishaw(file):
    #Renishaw file reading
    df=pd.read_csv(file,names=['x','y','lambda','int'],sep='\t')
    table=df.loc[:,'int'].values
    lambdas=df.loc[:,'lambda'].values
    lambda_0=lambdas[0]
    lambdas_one=lambdas[: (np.argwhere(lambdas==lambda_0)[1])[0]]
    
    X=df.iloc[::lambdas_one.shape[0],0].values
    Y=df.iloc[::lambdas_one.shape[0],1].values
      
    intensities=np.transpose(np.reshape(table,(X.shape[0],lambdas_one.shape[0])))
    lambdas=np.transpose(np.reshape(lambdas,(X.shape[0],lambdas_one.shape[0])))
    
    return X, Y, lambdas_one,intensities

def read_horiba(file):
    #Horiba file reading
    df=pd.read_csv(file,sep='\t')
    intensities=df.iloc[:,2:].values
    lambdas=df.columns.values[2:]
    X=df.iloc[:,0].values
    Y=df.iloc[:,1].values
    return X, Y, lambdas,intensities.T

def peak(X, Y, lambdas,intensities,function,Xrange,amp,Xmean,sigma,y0,A):
    #fitting
    if function=='gauss':
        fun=peak_shapes.create_gauss()
    elif function=='lorenz':
        fun=peak_shapes.create_lorenz()    
    results=np.empty(5)
    for d in np.ndindex(intensities.shape[1]):
        try:
            popt, pcov = curve_fit(fun, lambdas[Xrange[0]:Xrange[1]], 
                                   np.squeeze(intensities[Xrange[0]:Xrange[1],d]),
                                   p0=(amp,Xmean,sigma,y0,A))
        except RuntimeError:
            print("Error - curve_fit failed")
        results=np.vstack((results,popt))
    #maps
    if X[0]!=X[1]:
        n_X0=np.argwhere(X!=X[0])[0,0] # while main axis in x
        n_X1=int(X.shape[0]/n_X0)
    else:
        n_X1=np.argwhere(Y!=Y[0])[0,0] # while main axis in y
        n_X0=int(Y.shape[0]/n_X1)
    
    
    rmap=np.empty([n_X0,n_X1])

    for d in np.ndindex(results.shape[1]):
        rmap=np.dstack((rmap,results[1:,d].reshape(n_X0,n_X1)))
    
    return results,rmap
