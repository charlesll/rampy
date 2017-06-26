import unittest

import sys
sys.path.append("../../")

import numpy as np
import scipy

import gcvspline, rampy

class TestBaseline(unittest.TestCase):
    
    def baseline_fit(self):

        nb_points  =500
        x = np.linspace(100, 200, nb_points)
        y_true = 10.0 * np.exp(-np.log(2) * ((x-150.0)/15.0)**2) 

        noise = 0.1 * np.random.normal(size=nb_points)

        y_obs = y_true + noise + 0.1*x
    
        spectrum = np.transpose(np.vstack((x,y_obs)))
        roi = np.array([[100,110],[190,200]])
        ycalc, base, coefs = rampy.baseline(spectrum,roi,'gcvspline',2.0 )
        
        self.assertAlmostEqual((y_true+noise)-ycalc)
        
if __name__ == '__main__':
    unittest.main()




