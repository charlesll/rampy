Peak fitting
=============

Rampy does not offer a dedicated function for peak fitting. Instead, we invite users to use lmfit or scipy.optimize to perform peak fitting, which is basically the action to fit a model (sum of peaks) to your data.

Rampy offers functions for various peak shapes, including:
- gaussian peaks > ``rampy.gaussian``
- lorentzian peaks > ``rampy.lorentzian``
- pseudo-voigt peaks > ``rampy.pseudovoigt``
- pearson7 peaks > ``rampy.pearson7``

See the help for each function on this website. Those can be used to easily create a model that will be fitted to your spectrum.

We will upload soon an example of Bayesian peak fitting with a function integrated to rampy.

lmfit
-----

We provide an example of peak-fitting with the lmfit for instance. See `this notebook <https://github.com/charlesll/rampy/blob/master/examples/Raman_spectrum_fitting.ipynb>`_ for an example of use.

Calculate peak areas
--------------------

Peak area can be calculated using the ``rampy.peakarea`` function. For instance, for a Gaussian peak at 1100 cm-1 with an amplitude of 1.0 and a half-width at half maximum of 25.0, we can do:

.. code-block:: python

  position = 1100.
  amplitude = 1.0
  hwhm = 25.0

  area = rp.peakarea("gaussian", amp=amplitude, HWHM=hwhm)

Propagate uncertainties
-----------------------

The best way to propagate the uncertainties of your model is to directly use the ``uncertainties`` package, see the docs `here <https://pythonhosted.org/uncertainties/>`_.
