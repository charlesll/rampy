import unittest

import numpy as np
import scipy
from scipy.constants import c, h, k

import rampy as rp

class Testtlcorrection(unittest.TestCase):

    def test_tlcorrection(self):
        # testing long correction function
        x_for_long = np.array([20.,21.,22.,23.,24.,25.])
        y_for_long = np.array([1.0,1.0,1.0,1.0,1.0,1.0])

        nu0 = 1.0/(514.532)*1e9 #laser wavenumber at 514.532
        nu = 100.0*x_for_long # cm-1 to m-1
        T = 23.0+273.15 # the temperature in K

        x_long,long_res,eselong = rp.tlcorrection(x_for_long, y_for_long,23.0,514.532,correction = 'long',normalisation='area') # using the function
        t0 = nu0**3.0*nu/((nu0-nu)**4)
        t1= 1.0 - np.exp(-h*c*nu/(k*T)) # c in m/s  : t1 dimensionless
        long_calc= y_for_long*t0*t1 # pour les y
        long_calc = long_calc/np.trapz(long_calc,x_for_long) # area normalisation

        np.testing.assert_equal(long_res,long_calc)
        np.testing.assert_equal(x_for_long,x_long)

        x_long,long_res,eselong = rp.tlcorrection(x_for_long, y_for_long,23.0,514.532,correction = 'long',normalisation='no') # using the function
        t0 = nu0**3.0*nu/((nu0-nu)**4)
        t1= 1.0 - np.exp(-h*c*nu/(k*T)) # c in m/s  : t1 dimensionless
        long_calc= y_for_long*t0*t1 # pour les y

        np.testing.assert_equal(long_res,long_calc)
        np.testing.assert_equal(x_for_long,x_long)

if __name__ == '__main__':
    unittest.main()
