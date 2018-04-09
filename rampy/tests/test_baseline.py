import unittest

import numpy as np
np.random.seed(42) # fixing the random seed
import scipy

import rampy

#from matplotlib import pyplot as plt

class TestBaseline(unittest.TestCase):

    def test_baseline(self):

        x2 = np.arange(1,100,0.5)
        base_ori = 0.001*x2
        base_exp = rampy.funexp(x2,0.1,0.05,50.)
        base_log = rampy.funlog(x2,1.,1.,1.,1.)


        y_ori = 1.0 * np.exp(-np.log(2) * ((x2-50.0)/10.0)**2) + 0.05*np.random.randn(len(x2))

        y2 = base_ori + y_ori
        y_exp = base_exp+y_ori
        y_log = base_log+y_ori


        # need to define some fitting regions for the spline
        roi2 = np.array([[1,20],[80,100]])

        # calculating the baselines
        ycalc1, base1 = rampy.baseline(x2,y2,roi2,'poly',polynomial_order=1)
        #ycalc2, base2 = rampy.baseline(x2,y2,roi2,'gcvspline',s=0.1 )
        ycalc3, base3 = rampy.baseline(x2,y2,roi2,'unispline',s=1e0)
        ycalc4, base4 = rampy.baseline(x2,y2,roi2,'als',lam=10**7,p=0.05)
        ycalc5, base5 = rampy.baseline(x2,y2,roi2,'arPLS',lam=10**7,ratio=0.1)
        ycalc6, base6 = rampy.baseline(x2,y2,roi2,'exp',p0_exp=[0.1,0.1,45])

        # Testing the shapes
        np.testing.assert_equal(ycalc1.shape,base1.shape)
        #np.testing.assert_equal(ycalc2.shape,base2.shape)
        np.testing.assert_equal(ycalc3.shape,base3.shape)
        np.testing.assert_equal(ycalc4.shape,base4.shape)
        np.testing.assert_equal(ycalc5.shape,base5.shape)
        np.testing.assert_equal(ycalc6.shape,base6.shape)
        #np.testing.assert_equal(ycalc7.shape,base7.shape)

        # testing the baselines
        np.testing.assert_almost_equal(base_ori,base1[:,0],0)
        #np.testing.assert_almost_equal(base_ori,base2[:,0],0)
        np.testing.assert_almost_equal(base_ori,base3[:,0],0)
        np.testing.assert_almost_equal(base_ori,base4[:,0],0)
        np.testing.assert_almost_equal(base_ori,base5[:,0],0)
        #exp-log cases
        np.testing.assert_almost_equal(base_exp,base6[:,0],0)
        #np.testing.assert_almost_equal(base_log,base7[:,0],0)

        #testing the corrected data
        np.testing.assert_almost_equal(y_ori,ycalc1[:,0],1)
        #np.testing.assert_almost_equal(y_ori,ycalc2[:,0],0)
        np.testing.assert_almost_equal(y_ori,ycalc3[:,0],0)
        np.testing.assert_almost_equal(y_ori,ycalc4[:,0],0)
        np.testing.assert_almost_equal(y_ori,ycalc5[:,0],0)
        np.testing.assert_almost_equal(y_ori,ycalc6[:,0],0)

if __name__ == '__main__':
    unittest.main()




