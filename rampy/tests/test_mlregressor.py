import unittest

import numpy as np
np.random.seed(42)
import scipy

from scipy.stats import norm

import rampy as rp

class TestML(unittest.TestCase):

    def test_mlregressor(self):

        x = np.arange(0,600,1.0)
        nb_samples = 100 # number of samples in our dataset

        # true partial spectra
        S_1 = norm.pdf(x,loc=200.,scale=130.)
        S_2 = norm.pdf(x,loc=400,scale=70)
        S_true = np.vstack((S_1,S_2))

        #60 samples with random concentrations between 0 and 1
        C_ = np.random.rand(nb_samples) 
        C_true = np.vstack((C_,(1-C_))).T

        # We make some observations with random noise
        Obs = np.dot(C_true,S_true) + np.random.randn(nb_samples,len(x))*1e-4

        # new observations
        C_new_ = np.random.rand(10) #10 samples with random concentrations between 0 and 1
        C_new_true = np.vstack((C_new_,(1-C_new_))).T
        
        noise_new = np.random.randn(len(x))*1e-4
        Obs_new = np.dot(C_new_true,S_true) + noise_new

        model = rp.mlregressor(Obs,C_true[:,0].reshape(-1,1))

        for i in ["KernelRidge", "SVM", "LinearRegression", "Lasso", "ElasticNet", "NeuralNet", "BaggingNeuralNet"]:
            model.algorithm = i
            model.user_kernel = 'poly'
            model.fit()

            C_new_predicted = model.predict(Obs_new)

if __name__ == '__main__':
    unittest.main()
