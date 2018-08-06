import unittest

import numpy as np
import scipy

import rampy as rp

class TestSmooth(unittest.TestCase):

    def test_mixing(self):

        # dummy gaussians
        x = np.arange(0,100,1.0) # a dummy x axis
        ref1 = 50.0*np.exp(-1/2*((x-40)/20)**2) + np.random.randn(len(x)) # a gaussian with added noise
        ref2 = 70.0*np.exp(-1/2*((x-60)/15)**2) + np.random.randn(len(x)) # a gaussian with added noise

        # mixed signals
        F1_true = np.array([0.80,0.60,0.40,0.20])
        obs = np.dot(ref1.reshape(-1,1),F1_true.reshape(1,-1)) + np.dot(ref2.reshape(-1,1),(1-F1_true.reshape(1,-1)))

        # calculation
        F1_meas = rp.mixing_sp(obs,ref1,ref2)

        # assertion
        np.testing.assert_almost_equal(F1_true,F1_meas,decimal=3)

if __name__ == '__main__':
    unittest.main()
