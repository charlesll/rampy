import unittest
import numpy as np
import rampy
from scipy.integrate import simpson

class TestPeakArea(unittest.TestCase):
    def test_gaussian_peak(self):
        """Test Gaussian peak area calculation."""
        amp = 10.0
        HWHM = 2.0
        # Analytical formula for Gaussian peak area
        expected_area = np.sqrt(np.pi / np.log(2)) * amp * HWHM
        area, ese_area = rampy.peakarea("gaussian", amp=amp, HWHM=HWHM)
        self.assertAlmostEqual(area, expected_area, places=6)
        self.assertIsNone(ese_area)

    def test_gaussian_with_uncertainties(self):
        """Test Gaussian peak area calculation with uncertainties."""
        amp = 10.0
        HWHM = 2.0
        ese_amp = 0.5
        ese_HWHM = 0.1
        expected_area = np.sqrt(np.pi / np.log(2)) * amp * HWHM
        expected_ese_area = np.sqrt(
            ((np.sqrt(np.pi / np.log(2)) * HWHM) ** 2) * ese_amp**2 +
            ((np.sqrt(np.pi / np.log(2)) * amp) ** 2) * ese_HWHM**2
        )
        area, ese_area = rampy.peakarea("gaussian", amp=amp, HWHM=HWHM, ese_amp=ese_amp, ese_HWHM=ese_HWHM)
        self.assertAlmostEqual(area, expected_area, places=6)
        self.assertAlmostEqual(ese_area, expected_ese_area, places=6)

    def test_lorentzian_peak(self):
        """Test Lorentzian peak area calculation using numerical integration."""
        amp = 10.0
        HWHM = 2.0
        pos = 5.0
        x_int = np.linspace(pos - 10 * HWHM, pos + 10 * HWHM, 10000)
        y_int = rampy.lorentzian(x_int, amp=amp, freq=pos, HWHM=HWHM)
        expected_area = simpson(y_int, x=x_int)
        
        area, _ = rampy.peakarea("lorentzian", amp=amp, HWHM=HWHM, pos=pos)
        self.assertAlmostEqual(area, expected_area, places=6)

    def test_pseudovoigt_peak(self):
        """Test pseudo-Voigt peak area calculation using numerical integration."""
        amp = 10.0
        HWHM = 2.0
        pos = 5.0
        L_ratio = 0.5
        x_int = np.linspace(pos - 10 * HWHM, pos + 10 * HWHM, 10000)
        y_int = rampy.pseudovoigt(x_int, amp=amp, freq=pos, HWHM=HWHM, L_ratio=L_ratio)
        expected_area = simpson(y_int, x=x_int)

        area, _ = rampy.peakarea("pseudovoigt", amp=amp, HWHM=HWHM, pos=pos, L_ratio=L_ratio)
        self.assertAlmostEqual(area, expected_area, places=6)

    def test_pearson7_peak(self):
        """Test Pearson VII peak area calculation using numerical integration."""
        amp = 10.0
        HWHM = 2.0
        pos = 5.0
        a3 = 1.5
        x_int = np.linspace(pos - 10 * HWHM, pos + 10 * HWHM, 10000)
        y_int = rampy.pearson7(x_int, a0=amp, a1=pos, a2=HWHM, a3=a3)
        expected_area = simpson(y_int, x=x_int)

        area, _ = rampy.peakarea("pearson7", amp=amp, HWHM=HWHM, pos=pos, a3=a3)
        self.assertAlmostEqual(area, expected_area, places=6)

    def test_invalid_shape(self):
        """Test invalid shape raises NotImplementedError."""
        with self.assertRaises(NotImplementedError):
            rampy.peakarea("invalid_shape", amp=10.0, HWHM=2.0)

    def test_missing_parameters(self):
        """Test missing parameters raise ValueError."""
        
        # Missing 'pos' for non-Gaussian shapes
        with self.assertRaises(ValueError):
            rampy.peakarea("lorentzian", amp=10.0, HWHM=2.0)

    def test_invalid_L_ratio(self):
        """Test invalid L_ratio raises ValueError."""
        
        # L_ratio out of bounds for pseudo-Voigt peaks
        with self.assertRaises(ValueError):
            rampy.peakarea("pseudovoigt", amp=10.0, HWHM=2.0, pos=5.0, L_ratio=-1)

    def test_array_inputs(self):
        """Test array inputs for multiple peaks."""
        
        # Array inputs for Gaussian peaks
        amps = np.array([10.0, 15.0])
        HWHMs = np.array([2.0, 3.0])
        
        areas_expected = np.sqrt(np.pi / np.log(2)) * amps * HWHMs
        
        areas_calculated, _ = rampy.peakarea("gaussian", amp=amps.tolist(), HWHM=HWHMs.tolist())
        
        np.testing.assert_allclose(areas_calculated.flatten(), areas_expected.flatten(), rtol=1e-6)