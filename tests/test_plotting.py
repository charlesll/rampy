import unittest
import numpy as np
from plotly.graph_objects import Figure
from rampy import plot_spectrum

class TestPlotSpectra(unittest.TestCase):
    def setUp(self):
        """Set up test data."""
        # Spectrum
        self.x1 = np.linspace(0, 10, 500)
        self.y1 = np.sin(self.x1) + np.random.normal(0, 0.1, len(self.x1))
        self.baseline1_1 = 0.5 * np.ones_like(self.x1)
        self.baseline1_2 = 0.3 * np.sin(0.5 * self.x1) + 0.5
        self.smoothed1_1 = np.sin(self.x1)
        self.smoothed1_2 = np.sin(self.x1)


    def test_plot_spectra_no_baselines_or_smoothed_signals(self):
        """Test plotting spectra without baselines or smoothed signals."""
        fig = plot_spectrum(
            x=self.x1,
            y=self.y1,
            label="Spectrum"
        )
        
        # Check that the returned object is a Plotly figure
        self.assertIsInstance(fig, Figure)
        
        # Check that the correct number of traces is added (original spectra only)
        self.assertEqual(len(fig.data), 1)

    def test_plot_spectra_with_baselines_and_smoothed_signals(self):
        """Test plotting spectra with baselines and smoothed signals."""
        fig = plot_spectrum(
            x=self.x1,
            y=self.y1,
            baselines=[self.baseline1_1, self.baseline1_2],
            baseline_labels=["Baseline A", "Baseline B"],
            smoothed_signals=[self.smoothed1_1, self.smoothed1_2],
            smoothed_labels=["Smoothed A", "Smoothed B"],
            label=["Spectrum 1"]
        )
        
        # Check that the returned object is a Plotly figure
        self.assertIsInstance(fig, Figure)
        
        # Check that the correct number of traces is added (original + baselines + smoothers)
        total_traces = (
            len([self.y1]) + len([self.baseline1_1, self.baseline1_2]) + len([self.smoothed1_1, self.smoothed1_2])
        )
        self.assertEqual(len(fig.data), total_traces)

    def test_plot_spectra_with_missing_labels(self):
        """Test plotting spectra with missing labels for baselines and smoothers."""
        fig = plot_spectrum(
            x=self.x1,
            y=self.y1,
            baselines=[self.baseline1_1,],
            smoothed_signals=[self.smoothed1_1,],
            label="Spectrum"
        )
        
        # Check that the returned object is a Plotly figure
        self.assertIsInstance(fig, Figure)
        
        # Check that default labels are used for baselines and smoothers
        trace_names = [trace.name for trace in fig.data]
        expected_names = ["Spectrum - Original", "Spectrum - Baseline 1", "Spectrum - Smoothed 1"]
        
        self.assertEqual(trace_names[:len(expected_names)], expected_names)

if __name__ == "__main__":
    unittest.main()
