import unittest

import numpy as np
np.random.seed(42)
import scipy

from scipy.stats import norm

import rampy as rp

class TestMLC(unittest.TestCase):

    def test_mlclassificator(self):

        # The X axis
        x = np.arange(0, 1000, 1.0)

        # The perfect 5 signals
        spectra_1 = rp.gaussian(x, 10.0, 300., 25.) + rp.lorentzian(x, 15., 650., 50.)
        spectra_2 = rp.gaussian(x, 20.0, 350., 25.) + rp.gaussian(x, 25.0, 380., 20.) + rp.lorentzian(x, 15., 630., 50.)
        spectra_3 = rp.gaussian(x, 10.0, 500., 50.) + rp.lorentzian(x, 15.0, 520., 10.) + rp.gaussian(x, 25., 530., 3.)
        spectra_4 = rp.gaussian(x, 10.0, 100., 5.) + rp.lorentzian(x, 30.0, 110., 3.) + rp.gaussian(x, 5., 900., 10.)
        spectra_5 = rp.gaussian(x, 10.0, 600., 200.)

        # the number of observations of each signal
        number_of_spectra = 20

        # generating a dataset (will be shuffled later during the train-test split)
        dataset = np.hstack((np.ones((len(x),number_of_spectra))*spectra_1.reshape(-1,1),
                             np.ones((len(x),number_of_spectra))*spectra_2.reshape(-1,1),
                             np.ones((len(x),number_of_spectra))*spectra_3.reshape(-1,1),
                             np.ones((len(x),number_of_spectra))*spectra_4.reshape(-1,1),
                             np.ones((len(x),number_of_spectra))*spectra_5.reshape(-1,1)
                            )).T

        # add noise
        noise_scale = 2.0
        dataset = dataset + np.random.normal(scale=noise_scale,size=(len(dataset),len(x)))

        # create numeric labels
        labels =  np.vstack((np.tile(np.array([1]).reshape(-1,1),number_of_spectra),
                             np.tile(np.array([2]).reshape(-1,1),number_of_spectra),
                             np.tile(np.array([3]).reshape(-1,1),number_of_spectra),
                             np.tile(np.array([4]).reshape(-1,1),number_of_spectra),
                             np.tile(np.array([5]).reshape(-1,1),number_of_spectra),
                            )).reshape(-1,1)

        # new observations
        C_new_ = np.random.rand(10) #10 samples with random concentrations between 0 and 1
        C_new_true = np.vstack((C_new_,(1-C_new_))).T

        noise_new = np.random.randn(len(x))*1e-4
        Obs_new = np.dot(C_new_true,S_true) + noise_new

        explo = rp.mlexplorer(Obs)

        # we just test that it runs for now
        explo.algorithm = 'NMF'
        explo.nb_compo = 2
        explo.test_size = 0.3
        explo.scaler = "MinMax"
        explo.fit()
        explo.refit()

if __name__ == '__main__':
    unittest.main()
