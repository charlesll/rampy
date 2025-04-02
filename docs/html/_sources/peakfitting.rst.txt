.. _RST Fitting:

Peak fitting
=============

Rampy does not (yet) offer a dedicated function for peak fitting. Instead, we invite users to use lmfit or scipy.optimize to perform peak fitting, which is basically the action to fit a model (sum of peaks) to your data.

We provide an example of peak-fitting with the lmfit for instance here: :ref:`RST Notebooks` .

Peak shapes
--------------------

Rampy offers functions for various peak shapes, including:

* gaussian peaks > ``rampy.gaussian``
* lorentzian peaks > ``rampy.lorentzian``
* pseudo-voigt peaks > ``rampy.pseudovoigt``
* pearson7 peaks > ``rampy.pearson7``

.. automodule:: rampy.peak_shapes
   :members:
   :undoc-members:
   :show-inheritance:

Peak areas
----------

Peak area can be calculated using the ``rampy.peakarea`` function. 

.. automodule:: rampy.peak_area
   :members:
   :undoc-members:
   :show-inheritance:

Propagate uncertainties
-----------------------

The best way to propagate the uncertainties of your model is to directly use the ``uncertainties`` package, see the docs `here <https://pythonhosted.org/uncertainties/>`_.
