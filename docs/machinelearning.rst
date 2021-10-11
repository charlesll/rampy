Machine learning
================

Rampy offers three classes for performing classification, regression or unsupervised ML exploration of a set of spectra. Those are helpful ways to automatically perform usual scalling, test-train splitting and ML training of popular ML algorithms, and use `scikit-learn <https://scikit-learn.org/stable/>`_ in the background.

Those functions work (``rp.regressor`` was used in `this publication <https://doi.org/10.2138/am-2019-6887>`_) but still may evolve in the futur. For advanced, ML, I suggest using directly scikit-learn or other ML libraries.

Do not hesitate to ask for new features depending on your needs !

In the following, we assume you already performed library importation as

.. code-block:: python

  # importing rampy
  import rampy as rp
  # and for numpy we will respect the usual name:
  import numpy as np
  # for matplotlib
  import matplotlib.pyplot as plt

ML Exploration
--------------

The ``rampy.mlexplorer`` class allows performing PCA and NMF at the moment (work in progress). See the `example notebook <https://github.com/charlesll/rampy/blob/master/examples/Machine%20Learning%20Exploration.ipynb>`_.

ML regression
-------------

The ``rampy.regressor`` class allows linking variations in spectra to some known variations like the composition of a sample. See the `example notebook <https://github.com/charlesll/rampy/blob/master/examples/Machine%20Learning%20Regression.ipynb>`_.

ML classification
-----------------

The ``rampy.classificator`` class allows automatic recognition of substances/stuffs from spectra. See the `example notebook here <https://github.com/charlesll/rampy/blob/master/examples/MachineLearning_Classification.ipynb>`_.

Linear mixture
--------------

If you have two endmember spectra, you can use the ``rampy.mixing()`` function. See the `mixing example notebook here <https://github.com/charlesll/rampy/blob/master/examples/Mixing_spectra.ipynb>`_.

If you do not know the endmember spectra, then you may be interested in using directly the PyMCR library, see the documentation `here <https://pages.nist.gov/pyMCR/>`_ and an example notebook `here <https://github.com/usnistgov/pyMCR/blob/master/Examples/Demo.ipynb>`_. We used it in `this publication <https://doi.org/10.2138/am-2019-6887>`_, see the code `here <https://github.com/charlesll/rampy/blob/master/examples/Iron_AmMin_paper/Iron_MORB_code.ipynb>`_.
