import pandas as pd
import numpy as np
import functions
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

file='4338-NDC1-(31)_4X4.txt'

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
    
    plt.plot(lambdas[:,0],intensities[:,0])
    
    return X, Y, lambdas_one,intensities

def peak(data,function,energies):
    (X, Y, lambdas,intensities)=data
    if function=='gauss':
        fun=functions.create_gauss()
        fit_fun=lambda x: curve_fit(fun, lambdas, x)
    
    data=pd.DataFrame(intensities)
    #data.apply(lambda x: curve_fit(fun, lambdas, x))
    #plt.plot(lambdas[:,0],intensities[:,0])
    np.apply_along_axis(fit_fun, 0, intensities)
    
    #return curve_fit(fun, lambdas, intensities)
    #return np.array([curve_fit(fun, lambdas_one, xi) for xi in intensities],axis=1)
    


data=read_renishaw(file)
result=peak(data,'gauss',(100,200))
