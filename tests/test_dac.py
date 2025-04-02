import unittest
import numpy as np
import rampy

class TestPressureFunctions(unittest.TestCase):
    
    def test_ruby_T_corr(self):
        """Test temperature correction for Ruby pressure scale."""
        T_K = np.array([296., 300.0])
        expected_corrections = np.array([0.0, 0.0297924])  # Calculated manually
        corrections = rampy.ruby_T_corr(T_K)
        
        # Assert that corrections match expected values
        np.testing.assert_almost_equal(corrections, expected_corrections, decimal=5)

    def test_borate_T_corr(self):
        """Test temperature correction for Borate pressure scale."""
        T_K = np.array([296.0, 300.])
        expected_corrections = np.array([0.0, -0.00027423])  # Calculated manually
        corrections = rampy.borate_T_corr(T_K)
        
        # Assert that corrections match expected values
        np.testing.assert_almost_equal(corrections, expected_corrections, decimal=6)

    def test_pressure_sensor_ruby(self):
        """Test pressure calculation for Ruby scale."""
        lbd_nm = 695.22
        T_K = 300.
        expected_pressure = 1.87e3 * ((lbd_nm - rampy.ruby_T_corr(T_K) - 694.22) / 694.22) * (
            1 + 5.63 * ((lbd_nm - rampy.ruby_T_corr(T_K) - 694.22) / 694.22)
        )
        
        # Calculate pressure using function
        pressure = rampy.pressure_sensor(lbd_nm=lbd_nm, T_K=T_K)
        
        # Assert that calculated pressure matches expected value
        self.assertAlmostEqual(pressure, expected_pressure, places=5)

    def test_pressure_sensor_borate_datchi(self):
        """Test pressure calculation for Borate scale using Datchi."""
        lbd_nm = 686.44
        T_K = 300.0
        A, B, C = 3.989, 0.006915, 0.0166
        dlbd = lbd_nm - rampy.borate_T_corr(T_K) - 685.44
        
        # Expected pressure calculation
        expected_pressure = A * dlbd * ((1.0 + B * dlbd) / (1.0 + C * dlbd))
        
        # Calculate pressure using function
        pressure = rampy.pressure_sensor(lbd_nm=lbd_nm, T_K=T_K, scale="borate", borate_scale="datchi2007")
        
        # Assert that calculated pressure matches expected value
        self.assertAlmostEqual(pressure, expected_pressure, places=5)

    def test_pressure_sensor_invalid_scale(self):
        """Test invalid scale raises ValueError."""
        
        with self.assertRaises(ValueError):
            rampy.pressure_sensor(lbd_nm=695.22, T_K=300, scale="invalid")