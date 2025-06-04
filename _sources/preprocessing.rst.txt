.. _RST Preprocessing:

Preprocessing
=============

Rampy offers handful functions to preprocess your spectra. Do not hesitate to propose/ask for new functionalities!

Below you will find the documentation of the relevant functions. They are often used in the different example notebooks: :ref:`RST Notebooks`

Flip the X axis
---------------

Some spectra come with decreasing X values. Rampy offers a simple function to flip them. This can be necessary to e.g. resample them (as interpolation algorithms usually require the X values to increase). The following function does this. It even allows even arbitrary X value positions. It returns a sorted array with an  increasing x axis.

.. autoclass:: rampy.spectranization.flipsp
   :members:
   :undoc-members:
   :show-inheritance:

Shift the X axis 
----------------

You can shift the X axis from a given value using the following function: 

.. autoclass:: rampy.spectranization.shiftsp
   :members:
   :undoc-members:
   :show-inheritance:

Extract a portion / portions of a signal
----------------------------------------

You can use the function ``rampy.extract_signal()`` to do that (old version: ``rampy.get_portion_interest``)

.. autoclass:: rampy.baseline.extract_signal
   :members:
   :undoc-members:
   :show-inheritance:

Remove spikes
-------------

Spikes can be removed via the ``rampy.despiking()`` function. It takes as input X and Y values of a spectrum and a threshold. The threshold is the number of standard deviation above the mean noise value that a point must be to be considered as a spike. For instance, if the threshold is 3, then a point will be considered as a spike if it is 3 standard deviation above the mean of the noise. The function will then replace the spike by the mean of *k* points before and after the spike. 

.. autoclass:: rampy.spectranization.despiking
   :members:
   :undoc-members:
   :show-inheritance:

Resampling a spectrum
---------------------

We need sometime to resample a spectrum with a new X axis. ``rampy.resample()`` offers such ability. For instance,
we have a spectrum that has a X axis from 400 to 1300 cm-1, with points each 0.9 cm-1. We want the same but with an X axis with a value each cm-1. We can do for our spectrum:

.. autoclass:: rampy.spectranization.resample
   :members:
   :undoc-members:
   :show-inheritance:

Normalisation
-------------

Rampy provides the ``rampy.normalisation()`` function to normalise the Y values of a spectrum to

- the maximum intensity
- the trapezoidal area under the curve
- to min-max values of intensities

.. autoclass:: rampy.spectranization.normalise
   :members:
   :undoc-members:
   :show-inheritance:

Temperature and excitation line effects
---------------------------------------

Raman spectra may need correction from temperature and excitation line effects. See the review of `Brooker et al. 1988 <http://onlinelibrary.wiley.com/doi/10.1002/jrs.1250190202/abstract>`_ for details. rampy offers several way to do so with the ``rampy.tlexcitation()`` function.

.. autoclass:: rampy.tlcorrection.tlcorrection
   :members:
   :undoc-members:
   :show-inheritance:

Wavelength-wavenumber convertion
--------------------------------

The ``convert_x_units()`` function allows to convert your X values in nm in inverse cm, or the opposite! Do not hesitate to propose new ways to enrich it!

.. autoclass:: rampy.spectranization.convert_x_units
   :members:
   :undoc-members:
   :show-inheritance:

