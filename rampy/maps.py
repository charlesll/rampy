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
    
    plt.plot(lambdas_one,intensities[:,0])
    
    return X, Y, lambdas_one,intensities

def peak(data,function,Xrange,Xmean,sigma,amp,y0,A):
    (X, Y, lambdas,intensities)=data
    if function=='gauss':
        fun=functions.create_gauss()
    elif function=='lorenz':
        fun=functions.create_lorenz()    
        print('correct')
        
    fit_fun=lambda x: curve_fit(fun, lambdas[Xrange[0]:Xrange[1]], x[Xrange[0]:Xrange[1]],
                                p0=(amp,Xmean,sigma,y0,A))
    draw_fun=lambda x: plt.plot(lambdas[Xrange[0]:Xrange[1]], x[Xrange[0]:Xrange[1]])
    
    #df=pd.DataFrame(intensities)
    #data=df.apply(lambda x: curve_fit(fun, lambdas[Xrange[0]:Xrange[1]], x[Xrange[0]:Xrange[1]],
                       #         p0=(amp,Xmean,sigma,y0,A)))
    #return data

    for d in np.ndindex(intensities.shape[0]):
        plt.plot(lambdas[Xrange[0]:Xrange[1]],intensities[Xrange[0]:Xrange[1],d])
        print(lambdas[Xrange[0]:Xrange[1]].shape)
        print(intensities[Xrange[0]:Xrange[1],d].shape)
        try:
            popt, pcov = curve_fit(fun, lambdas[Xrange[0]:Xrange[1]], 
                                   np.squeeze(intensities[Xrange[0]:Xrange[1],d]),
                                   p0=(amp,Xmean,sigma,y0,A))
        except RuntimeError:
            print("Error - curve_fit failed")
    #data=pd.DataFrame(intensities)
    #data.apply(lambda x: curve_fit(fun, lambdas, x))
    #plt.plot(lambdas[Xrange[0]:Xrange[1]],intensities[Xrange[0]:Xrange[1],0])
    #plt.plot(lambdas[Xrange[0]:Xrange[1]],fun(lambdas[Xrange[0]:Xrange[1]],
    #                                          amp,Xmean,sigma,y0,A))
    #data1,data2=np.apply_along_axis(fit_fun, 0, intensities)
    #np.apply_along_axis(draw_fun, 0, intensities)
    
    #return curve_fit(fun, lambdas, intensities)
   # return np.array([curve_fit(fun, lambdas_one, xi) for xi in intensities],axis=1)
    


data=read_renishaw(file)
result=peak(data,'lorenz',[0,200],3700,40,2000,6000,0)
