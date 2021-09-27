# RamPy
=======

Copyright (2015-2021) C. Le Losq and co.

lelosq@ipgp.fr

[![Build Status](https://travis-ci.org/charlesll/rampy.svg?branch=master)](https://travis-ci.org/charlesll/rampy) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1168729.svg)](https://doi.org/10.5281/zenodo.1168729) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/charlesll/rampy.git/master)

Rampy is a Python library that aims at helping processing spectroscopic data, such as Raman, Infrared or XAS spectra. It offers, for instance, functions to subtract baselines as well as to stack, resample or smooth spectra. It aims at facilitating the use of Python in processing spectroscopic data. It integrates within a workflow that uses Numpy/Scipy as well as optimisation libraries such as lmfit or emcee, for instance.

The /examples/ folder contain various examples.

# REQUIREMENTS

Rampy is tested on Python 3.8 (see Travis badge; no garantee that it works on other Python versions)

The following libraries are required and indicated in setup.cfg:

- Scipy
- Numpy >= 1.12
- sklearn
- pandas & xlrd

Optional dependencies:

- gcvspline (you need a working FORTRAN compiler for its installation. To avoid this problem under Windows, wheels for Python 2.7, 3.4 and 3.6 are provided for 64 bit Windows, and a wheel for Python 3.6 is provided for Windows 32 bits. If installation fails, please check if is due to a fortran compiler issue.)

- xlrd and matplotlib

*Installation of gcvspline as well as matplotlib and xlrd are necessary for use of the `rampy.rameau()` class.*

- cvxpy v 1.1 or higher. As for gcvspline, the installation of cvxpy can cause problems for Windows users due to missing compiler. See instructions from cvxpy in this case.

*Installation of cvxpy is necessary for use of the `rampy.mixing()` class.*

Additional libraries for model fitting may be wanted:

- lmfit & aeval (http://cars9.uchicago.edu/software/python/lmfit/)
- emcee

# INSTALLATION

Install with pip:

  `pip install rampy`

If you want to use gcvspline and cvxpy, also install it:

  `pip install gcvspline`

  `pip install cvxpy`

# EXAMPLES

Given a signal [x y] containing a peak, and recorded in a text file myspectrum.txt.

You can import it, remove a automatic background, plot the result, and print the centroid of the peak as:

```
import matplotlib.pyplot as plt
import numpy as np
import rampy as rp

spectrum = np.genfromtxt("myspectrum.txt")

bir = np.array([[0,100., 200., 1000]]) # the frequency regions devoid of signal, used by rp.baseline()
y_corrected, background = rp.baseline(spectrum[:,0],spectrum[:,1],bir,"arPLS",lam=10**10)

plt.figure()
plt.plot(spectrum[:,0],spectrum[:,1],"k",label="raw data")
plt.plot(spectrum[:,0],background,"k",label="background")
plt.plot(spectrum[:,0],y_corrected,"k",label="corrected signal")
plt.show()

print("Signal centroid is %.2f" % rp.centroid(spectrum[:,0],y_corrected))
```

See the /example folder for further examples.

# Other packages

rampy can be used also to analyse the output of the [RADIS](https://radis.readthedocs.io/en/latest/) package.

See for instance https://github.com/charlesll/rampy/issues/13

Updated September 2021
