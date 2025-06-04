.. _RST Fitting:

Peak fitting
=============

Rampy will soon offer a dedicated function for peak fitting: ``rampy.fit_peaks`` will use `curve_fit` from `scipy.optimize` to fit your spectrum with a model that is a sum of peaks.

Users can also use dedicated interfaces such as lmfit. We provide an example of peak-fitting with the lmfit for instance here: :ref:`RST Notebooks` .

Peak shapes
-------------

Rampy offers functions for various peak shapes, including:

* gaussian peaks > ``rampy.gaussian``
* lorentzian peaks > ``rampy.lorentzian``
* pseudo-voigt peaks > ``rampy.pseudovoigt``
* pearson7 peaks > ``rampy.pearson7``

.. autoclass:: rampy.peak_shapes.gaussian
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: rampy.peak_shapes.lorentzian
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: rampy.peak_shapes.pseudovoigt
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: rampy.peak_shapes.pearson7
   :members:
   :undoc-members:
   :show-inheritance:

Peak areas
----------

Peak area can be calculated using the ``rampy.area_peak`` function, using analytical functions.

Note that the old function ``rampy.peakarea`` is still available but it will be removed in a near future. It uses numerical integration and is thus less precise than the ``rampy.area_peak`` function.

.. automodule:: rampy.peak_area.area_peak
   :members:
   :undoc-members:
   :show-inheritance:

Propagate uncertainties
-----------------------

A good way to propagate the uncertainties of your model is to directly use the ``uncertainties`` package, see the docs `here <https://pythonhosted.org/uncertainties/>`_.
