.. RamPy documentation master file, created by
   sphinx-quickstart on Fri Oct  8 15:22:27 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RamPy's documentation!
=================================

Copyright (2015-2024) C. Le Losq and co.

Charles Le Losq, Institut de physique du globe de Paris, University of Paris
lelosq@ipgp.fr

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1168729.svg
   :target: https://doi.org/10.5281/zenodo.1168729

The Rampy project
-----------------

Rampy is a Python library that aims at helping processing spectroscopic data, such as Raman, Infrared or XAS spectra.
It offers, for instance, functions to subtract baselines as well as to stack, resample or smooth spectra.
It aims at facilitating the use of Python in processing spectroscopic data.

The project is hosted on Github (`here! <https://github.com/charlesll/rampy>`_), please do not hesitate to propose new features or simply ideas!

Integration with other packages
-------------------------------

Rampy integrates within a workflow that uses Numpy/Scipy/Matplotlib as well as optimisation libraries such as lmfit, emcee, PyMC3 for instance.

Rampy can be used also to analyse the output of the `RADIS <https://radis.readthedocs.io/en/latest/>`_ package. See the discussion about this `here <https://github.com/charlesll/rampy/issues/13>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   philosophy
   installation
   firststeps
   preprocessing
   peakfitting
   machinelearning
   rampy
