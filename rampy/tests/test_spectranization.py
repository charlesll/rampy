# -*- coding: utf-8 -*-
import unittest

import numpy as np
import scipy

import gcvspline, rampy

#from matplotlib import pyplot as plt

class TestSpectranization(unittest.TestCase):
    
    def test_flipsp(self):
        
        x = np.arange(1,100,1)
        y = 2*x
        
        sp = np.transpose(np.vstack([x,y]))
        sp_flipped = np.flip(sp,0)
        sp_unflipped = rampy.flipsp(sp_flipped)
        # Testing 
        np.testing.assert_equal(sp_unflipped,sp)
        
    def test_resample(self):
        
        x = np.arange(1,100,1)
        x2 = np.arange(4,95,0.5)
        y = 2*x
        
        f = scipy.interpolate.interp1d(x,y)
        y_new_1 = f(x2)
        
        y_new_2 = rampy.resample(x,y,x2)
        # Testing 
        np.testing.assert_equal(y_new_1,y_new_2)
        
                
if __name__ == '__main__':
    unittest.main()




