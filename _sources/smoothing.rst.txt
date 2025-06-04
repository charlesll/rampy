Smoothing and filtering
========================

Functions are provided to smooth and filter spectra. Smoothing is useful to remove noise from spectra, while filtering can be used to remove specific frequency components. 

Below you will find the documentation of the relevant functions, and have a look at the example notebooks too: :ref:`RST Notebooks`


Smoothing
---------

.. image:: ./images/smooth.png
  :width: 400

Rampy provides the ``rampy.smooth()`` function to smooth a spectrum. 10 different algorithms are available, see the `notebook on Github <https://github.com/charlesll/rampy/blob/master/examples/Smoothing.ipynb>`_ for examples of use.

You can also directly use the Whittaker smoother. This can be useful if you want to use weights, which can allow you to avoid smoothing some regions for instance, or it can offer you more control.

You will find example notebooks below as well as the docstrings of the functions!

.. autoclass:: rampy.filters.smooth
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: rampy.filters.whittaker
   :members:
   :undoc-members:
   :show-inheritance:

Filtering
---------

The ``rampy.spectrafilter()`` function allows you to filter your spectra using a Butterworth filter. Low, high, bandstop and bandpass filters are possible. Can be useful to remove wavelets from FTIR spectra for instance!

.. autoclass:: rampy.filters.spectrafilter
   :members:
   :undoc-members:
   :show-inheritance:

