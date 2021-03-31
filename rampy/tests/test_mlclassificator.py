import unittest

import numpy as np
np.random.seed(42)

import scipy
from scipy.stats import norm

from sklearn.utils import shuffle

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

        X = dataset
        y = labels
        
        names = ["Nearest Neighbors", "Linear SVM", "RBF SVM", "Gaussian Process",
                 "Decision Tree", "Random Forest", "Neural Net", "AdaBoost",
                 "Naive Bayes", "QDA"]

        # we just test that it runs for now
        # initiate model
        MLC = rp.mlclassificator(X,y,scaling=False,test_size=0.33)

        # iterate over classifiers
        for name in names:
            MLC.algorithm = name
            MLC.fit()
            score = MLC.model.score(MLC.X_test, MLC.y_test)

if __name__ == '__main__':
    unittest.main()
