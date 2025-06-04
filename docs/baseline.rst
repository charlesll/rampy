Baseline
========

Rampy allows you to fit polynomial, spline, generalized cross-validated spline, logarithms, exponential, ALS, arPLS, drPLS, rubberband and whittaker baselines to your spectra, in order to remove the background.

.. image:: ./images/baseline.png
  :width: 500

Below you will find the documentation of the relevant functions, and have a look at the example notebooks too: :ref:`RST Notebooks`

.. autofunction:: rampy.baseline

For some baseline types, individual baseline functions also can be used, if you are interested:

.. autofunction:: rampy.rubberband_baseline

.. autofunction:: rampy.als_baseline

.. autofunction:: rampy.arPLS_baseline

.. autofunction:: rampy.drPLS_baseline