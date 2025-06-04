Maps
====

Rampy offers the possiblity to analyse maps from Raman spectrometers using the ``rampy.maps`` module.

Line and 2D Maps saved in CSV or TXT format from Horiba and Renishaw spectrometers can be loaded, and various data treatments can be performed afterwards.

Below you will find the documentation of the ``rampy.maps`` module, and have a look at the example notebooks too for further examples of use: :ref:`RST Notebooks`

.. module:: rampy.maps
   :synopsis: Module for treating maps of Raman spectra

.. class:: maps(file, spectrometer_type="horiba", map_type="2D", despiking=False, neigh=4, threshold=3)

   A class for treating maps of Raman spectra.
   
   :param file: Filename including path
   :type file: str
   :param spectrometer_type: Type of spectrometer, choose between "horiba" or "renishaw"
   :type spectrometer_type: str, optional, default="horiba"
   :param map_type: Type of map, choose between "2D" or "1D"
   :type map_type: str, optional, default="2D"
   :param despiking: Whether to apply despiking to remove spikes from the spectra
   :type despiking: bool, optional, default=False
   :param neigh: Number of points around spikes to select for calculating average value for despiking
   :type neigh: int, optional, default=4
   :param threshold: Multiplier of sigma (RMSE) for spike detection
   :type threshold: int, optional, default=3
   
   .. note::
      The class automatically flips the frequency axis if it's decreasing to ensure it's always increasing.

   .. method:: background(bir, method="poly", **kwargs)
   
      Correct a background from the initial signal I on a map using rampy.baseline.
      
      :param bir: Arrays of the background interpolation regions
      :type bir: ndarray
      :param method: Method for baseline correction, see rampy.baseline documentation for available methods
      :type method: str, optional, default="poly"
      :param **kwargs: Additional parameters for rampy.baseline()
      
      :return: None, but sets self.I_background and self.I_corrected attributes
      
   .. method:: normalise(y, method="intensity")
   
      Normalise the spectra to their max intensity, their area, or min-max normalisation.
      
      :param y: The intensities to normalise (e.g., self.I_corrected)
      :type y: ndarray
      :param method: Method used, choose between "area", "intensity", "minmax"
      :type method: str, optional, default="intensity"
      
      :return: None, but sets self.I_normalised attribute
      
   .. method:: smooth(y, method="whittaker", **kwargs)
   
      Uses the smooth function of rampy to smooth the signals.
      
      :param y: The intensities to smooth (e.g., self.I_corrected)
      :type y: ndarray
      :param method: Method for smoothing the signal
      :type method: str, optional, default="whittaker"
      :param **kwargs: Additional parameters for the smoothing method
      
      :return: None, but sets self.y_smoothed attribute
      
   .. method:: centroid(y, region_to_investigate)
   
      Calculate the centroid in a given region of interest.
      
      :param y: The intensities to analyze (e.g., self.I_normalised)
      :type y: ndarray
      :param region_to_investigate: The x values of regions where the centroid will be measured
      :type region_to_investigate: 1x2 array
      
      :return: None, but sets self.centroid_position attribute
      
   .. method:: intensity(y, region_to_investigate)
   
      Get the maximum intensity in the region to investigate.
      
      :param y: The intensities to consider (e.g., self.I_normalised)
      :type y: ndarray
      :param region_to_investigate: The x values of regions where the intensity will be measured
      :type region_to_investigate: 1x2 array
      
      :return: None, but sets self.I_max attribute
      
   .. method:: area(y, region_to_investigate)
   
      Get the area under the curve in the region to investigate.
      
      :param y: The intensities to consider (e.g., self.I_normalised)
      :type y: ndarray
      :param region_to_investigate: The x values of regions where the area will be measured
      :type region_to_investigate: 1x2 array
      
      :return: None, but sets self.A attribute
      
   .. method:: intensity_ratio(y, region_to_investigate)
   
      Get the intensity ratio between two regions of interest.
      
      :param y: The intensities to consider (e.g., self.I_normalised)
      :type y: ndarray
      :param region_to_investigate: The x values of regions where the intensity ratios will be measured
      :type region_to_investigate: 2x2 array
      
      :return: None, but sets self.I_ratio attribute
      
   .. method:: area_ratio(y, region_to_investigate)
   
      Get the area ratio between two regions of interest.
      
      :param y: The intensities to consider (e.g., self.I_normalised)
      :type y: ndarray
      :param region_to_investigate: The x values of regions where the areas will be measured
      :type region_to_investigate: 2x2 array
      
      :return: None, but sets self.A_ratio attribute

.. function:: read_renishaw(file)

   Read Renishaw csv maps.
   
   :param file: Filename including path
   :type file: str
   
   :return: Tuple containing (X, Y, lambdas, intensities)
   :rtype: tuple
   
   - X: m by n array of X positions
   - Y: m by n array of Y positions
   - lambdas: m by n array of Raman shifts
   - intensities: m by n array of intensities

.. function:: read_horiba(file, map_type="2D")

   Read Horiba csv maps (1D, 2D).
   
   :param file: Filename including path
   :type file: str
   :param map_type: 1D map (line) or 2D map
   :type map_type: str, optional, default="2D"
   
   :return: For 2D maps: Tuple containing (X, Y, lambdas, intensities)
   :return: For 1D maps: Tuple containing (X, lambdas, intensities)
   :rtype: tuple
   
   :raises ValueError: If map_type is not '1D' or '2D'

.. function:: peak(X, Y, lambdas, intensities, function, Xrange, amp, Xmean, sigma, y0, A)

   Fit peaks in a map. Work in progress.
   
   :param X: X positions
   :type X: ndarray
   :param Y: Y positions
   :type Y: ndarray
   :param lambdas: Raman shifts
   :type lambdas: ndarray
   :param intensities: Intensities
   :type intensities: ndarray
   :param function: Type of function to fit ('gauss' or 'lorenz')
   :type function: str
   :param Xrange: Range of X values to fit
   :type Xrange: list or ndarray
   :param amp: Initial amplitude
   :type amp: float
   :param Xmean: Initial mean position
   :type Xmean: float
   :param sigma: Initial width
   :type sigma: float
   :param y0: Initial y-offset
   :type y0: float
   :param A: Initial area
   :type A: float
   
   :return: Tuple containing (results, rmap)
   :rtype: tuple
