import pandas as pd
import numpy as np
import functions
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def read_renishaw(file):
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

def peak(data,function,Xrange,Xmean,sigma,amp,y0,A):
    (X, Y, lambdas,intensities)=data
    if function=='gauss':
        fun=functions.create_gauss()
    elif function=='lorenz':
        fun=functions.create_lorenz()    
    results=np.empty(5)
    for d in np.ndindex(intensities.shape[1]):
        try:
            popt, pcov = curve_fit(fun, lambdas[Xrange[0]:Xrange[1]], 
                                   np.squeeze(intensities[Xrange[0]:Xrange[1],d]),
                                   p0=(amp,Xmean,sigma,y0,A))
        except RuntimeError:
            print("Error - curve_fit failed")
        results=np.vstack((results,popt))
    return results

#data=read_renishaw(file)

