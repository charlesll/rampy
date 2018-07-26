# RamPy
=====

Copyright (2015-2018) C. Le Losq.

charles.lelosq@anu.edu.au

[![Build Status](https://travis-ci.org/charlesll/rampy.svg?branch=master)](https://travis-ci.org/charlesll/rampy)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1168730.svg)](https://doi.org/10.5281/zenodo.1168730)

Rampy is a Python library that aims at helping processing spectroscopic data, such as Raman, Infrared or XAS spectra. It offers, for instance, functions to subtract baselines as well as to stack, resample or smooth spectra. It aims at facilitating the use of Python in processing spectroscopic data. It integrates within a workflow that uses Numpy/Scipy as well as optimisation libraries such as lmfit or emcee, for instance.

The /examples/ folder contain various examples.

# REQUIREMENTS

Rampy is tested on Python 2.7 and 3.6 (see Travis badge; no garantee that it works on other Python versions)

The following libraries are required and indicated in setup.py:

- Scipy
- Numpy >= 1.12
- sklearn
- pandas
- cvxpy

Optional dependencies:

- gcvspline (you need a working FORTRAN compiler for its installation. To avoid this problem under Windows, wheels for Python 2.7, 3.4 and 3.6 are provided for 64 bit Windows, and a wheel for Python 3.6 is provided for Windows 32 bits. If installation fails, please check if is due to a fortran compiler issue.)

*Installation of gcvspline is necessary for use of the rampy.rameau() class.*

Additional libraries for model fitting may be wanted:

- lmfit & aeval (http://cars9.uchicago.edu/software/python/lmfit/)
- emcee

# INSTALLATION

Install with pip:

  `pip install rampy` 

If you want to use gcvspline, also install it:

  `pip install gcvspline`
 
# EXAMPLES

See the /example folder.

Updated July 2018





