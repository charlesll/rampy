Preprocessing
=============

Rampy offers handful functions to preprocess your spectra. Do not hesitate to propose/ask for new functionalities!

In the following, we assume you already performed library importation as

.. code-block:: python

  # importing rampy
  import rampy as rp
  # and for numpy we will respect the usual name:
  import numpy as np
  # for matplotlib
  import matplotlib.pyplot as plt

Flip the X axis
---------------

Some spectra come with decreasing X values. Rampy offers a simple function to flip them. This can be necessary to resample them (as interpolation algorithms usually require the X values to increase). Here for our spectrum we can do:

.. code-block:: python

  spectrum_increasing = rp.flipsp(spectrum)

Remove spikes
-------------

Spikes can be removed via the rp.despiking() function. It takes as input X and Y values of a spectrum and a threshold. 
The threshold is the number of standard deviation above the mean noise value that a point must be to be considered as a spike. 
For instance, if the threshold is 3, then a point will be considered as a spike if it is 3 standard deviation above the mean of the noise. 
The function will then replace the spike by the mean of *k* points before and after the spike. 

The function returns the new Y values.

Resampling a spectrum
---------------------

We need sometime to resample a spectrum with a new X axis. ``rampy.resample()`` offers such ability. For instance,
we have a spectrum that has a X axis from 400 to 1300 cm-1, with points each 0.9 cm-1. We want the same but with an X axis with a value each cm-1. We can do for our spectrum:

.. code-block:: python

  x_new = np.arange(400., 1300., 1.0) # we generate the new X values with numpy.arange()
  y_new = rp.resample(spectrum[:,0], spectrum[:,1], x_new)

Sometime, this can fail because we have a new X point at higher or lower values than the original X axis.
This is simply solved by asking ``rampy.resample()`` to extrapolate:

.. code-block:: python

  x_new = np.arange(400., 1300., 1.0) # we generate the new X values with numpy.arange()
  y_new = rp.resample(spectrum[:,0], spectrum[:,1], x_new, fill_value="extrapolate")

``rampy.resample`` uses ``scipy.interpolate.interp1d``. It can takes all the arguments that the latter function takes, so see its `doc <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html>`_ for details!

We can eventually create a new array with the x_new and y_new like:

.. code-block:: python

  spectrum_resample = np.vstack((x_new,y_new)).T

Normalisation
-------------

Rampy provides the ``rampy.normalisation()`` function to normalise the Y values of a spectrum to

- the maximum intensity
- the trapezoidal area under the curve
- to min-max values of intensities

See the `notebook on Github for examples <https://github.com/charlesll/rampy/blob/master/examples/Normalisation.ipynb>`_

Smoothing
---------

Rampy provides the ``rampy.smooth()`` function to smooth a spectrum. 10 different algorithms are available, see the `notebook on Github <https://github.com/charlesll/rampy/blob/master/examples/Smoothing.ipynb>`_ for examples of use.

.. image:: ./images/smooth.png
  :width: 400

Baseline
--------

Rampy allows you to fit polynomial, spline, generalized cross-validated spline, logarithms, exponential, ALS, arPLS, drPLS and rubberband baselines to your spectra, in order to remove the background.

See `this notebook <https://github.com/charlesll/rampy/blob/master/examples/Baseline.ipynb>`_ for an example of use, and the help of ``rampy.baseline()`` for details on each algorithm.

.. image:: ./images/baseline.png
  :width: 500

Temperature and excitation line effects
---------------------------------------

Raman spectra may need correction from temperature and excitation line effects. See the review of `Brooker et al. 1988 <http://onlinelibrary.wiley.com/doi/10.1002/jrs.1250190202/abstract>`_ for details. rampy offers several way to do so with the ``rampy.tlexcitation()`` function. See its documentation for details.
