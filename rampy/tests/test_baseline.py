import unittest

import sys
sys.path.append("../../")

import numpy as np
import scipy

import gcvspline, rampy

class TestBaseline(unittest.TestCase):
    
    def test_baseline(self):

        x2 = np.arange(1,100,0.5)
        y2 = 2*x2 + 20.0 * np.exp(-np.log(2) * ((x2-50.0)/10.0)**2)

        # need to define some fitting regions for the spline
        roi2 = np.array([[1,30],[80,100]])

        # calculating the baselines
        ycalc1, base1 = rampy.baseline(x2,y2,roi2,'poly',polynomial_order=1)
        ycalc2, base2 = rampy.baseline(x2,y2,roi2,'gcvspline',s=0.1 )
        ycalc3, base3 = rampy.baseline(x2,y2,roi2,'unispline',s=1e0)
        ycalc4, base4 = rampy.baseline(x2,y2,roi2,'als',lam=10**7,p=0.05)
        ycalc5, base5 = rampy.baseline(x2,y2,roi2,'arPLS',lam=10**7s,ratio=0.1)
        
        
        np.testing.assert_almost_equal(2*x2,base1[:,0],0)
        np.testing.assert_almost_equal(2*x2,base2[:,0],0)
        np.testing.assert_almost_equal(2*x2,base3[:,0],0)
        np.testing.assert_almost_equal(2*x2,base4[:,0],0)
        np.testing.assert_almost_equal(2*x2,base5[:,0],0)
        
if __name__ == '__main__':
    unittest.main()




