import unittest

import numpy as np
import scipy

import rampy as rp

#from matplotlib import pyplot as plt

class TestSmooth(unittest.TestCase):

    def test_smooth(self):

        nb_points  = 200
        x = np.linspace(50, 600, nb_points)

        # gaussian peaks
        p1 = 20.0 * np.exp(-np.log(2) * ((x-150.0)/15.0)**2)
        p2 = 100.0 * np.exp(-np.log(2) * ((x-250.0)/20.0)**2)
        p3 = 50.0 * np.exp(-np.log(2) * ((x-450.0)/50.0)**2)
        p4 = 20.0 * np.exp(-np.log(2) * ((x-350.0)/300)**2)
        p5 = 30.0 * np.exp(-np.log(2) * ((x-460.0)/5.0)**2)

        # background: a large gaussian + linear
        bkg = 60.0 * np.exp(-np.log(2) * ((x-250.0)/200.0)**2) + 0.1*x

        #noise
        noise = 2.0 * np.random.normal(size=nb_points)

        #observation
        y = p1 + p2 + p3 + p4 + p5 + noise +bkg

        # calculating the baselines
        #y_smo_1 = rp.smooth(x,y,method="GCVSmoothedNSpline")
        #y_smo_2 = rp.smooth(x,y,method="DOFSmoothedNSpline")
        #y_smo_3 = rp.smooth(x,y,method="MSESmoothedNSpline")
        y_smo_4 = rp.smooth(x,y,method="savgol",window_length=5,polyorder=2)
        y_smo_5 = rp.smooth(x,y,method="whittaker",Lambda=10**0.5)
        y_smo_6 = rp.smooth(x,y,method="flat",window_length=5)
        y_smo_7 = rp.smooth(x,y,method="hanning",window_length=5)
        y_smo_8 = rp.smooth(x,y,method="hamming",window_length=5)
        y_smo_9 = rp.smooth(x,y,method="bartlett",window_length=5)
        y_smo_10 = rp.smooth(x,y,method="blackman",window_length=5)

        # Testing the shapes
        #np.testing.assert_equal(y_smo_1.shape,y.shape)
        #np.testing.assert_equal(y_smo_2.shape,y.shape)
        #np.testing.assert_equal(y_smo_3.shape,y.shape)
        np.testing.assert_equal(y_smo_4.shape,y.shape)
        np.testing.assert_equal(y_smo_5.shape,y.shape)
        np.testing.assert_equal(y_smo_6.shape,y.shape)
        np.testing.assert_equal(y_smo_7.shape,y.shape)
        np.testing.assert_equal(y_smo_8.shape,y.shape)
        np.testing.assert_equal(y_smo_9.shape,y.shape)
        np.testing.assert_equal(y_smo_10.shape,y.shape)

        #testing the y values, difference should be less than a percent
        #self.assertTrue(np.sum(np.abs(y-y_smo_1))/np.sum(y)<0.1)
        #self.assertTrue(np.sum(np.abs(y-y_smo_2))/np.sum(y)<0.1)
        #self.assertTrue(np.sum(np.abs(y-y_smo_3))/np.sum(y)<0.1)
        self.assertTrue(np.sum(np.abs(y-y_smo_4))/np.sum(y)<0.1)
        self.assertTrue(np.sum(np.abs(y-y_smo_5))/np.sum(y)<0.1)
        self.assertTrue(np.sum(np.abs(y-y_smo_6))/np.sum(y)<0.1)
        self.assertTrue(np.sum(np.abs(y-y_smo_7))/np.sum(y)<0.1)
        self.assertTrue(np.sum(np.abs(y-y_smo_8))/np.sum(y)<0.1)
        self.assertTrue(np.sum(np.abs(y-y_smo_9))/np.sum(y)<0.1)
        self.assertTrue(np.sum(np.abs(y-y_smo_10))/np.sum(y)<0.1)

if __name__ == '__main__':
    unittest.main()
