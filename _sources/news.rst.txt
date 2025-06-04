Rampy News
==========

Copyright (c) 2014-2025 Dr. Charles Le Losq et al.

email: lelosq@ipgp.fr

Licence: see LICENCE.md

Wanted features (planned)
-------------------------

- Addition of more unsupervised machine learning techniques

- Classification via 1D CNN and other algorithms using Keras or Pytorch

- peak fitting (in particular for maps): working on a Scikit-Learn API

0.6.2 (stable)
----------------

  - Fix the documentation.

0.6.1 (stable)
----------------

Additions:

  - area_peak(): function to calculate analytically the area of peaks in a spectrum.

Modifications:

  - Minimal Python version is now 3.10
  - normalise(): "area" normalisation now uses the simpson() function from scipy.integrate, instead of the trapezoidal rule. This should be more precise.

Improvements:

  - maps(): documentation of the class has been improved.

Fix:

  - baseline(): a problem with the roi has been raised (issue #34). This is now fixed, roi argument is completely optional.
  - resample(): fill_value was not working properly (not the good default option). It is now fixed. resample() will extrapolate automatically.

0.6.0
----------------

Improvements:

  - baseline(): we use the make_smoothing_spline function from Scipy for the GCVSmoothedNSpline smoothing mode (smooth function) and the gcvspline baseline, in replaced of the call to the gcvspline library.  make_smoothing_spline is a reimplenetation of the Woltring Fortran code in gcvspline library, so modulo numerical errors, it returns the same results. This offers now the gcvspline baseline to all, without problem of Fortran compilation.
  - baseline(): bir are now an optional argument in rampy.baseline (because they do not appear in als, arPLS and drPLS algorithms)
  - peakarea(): function is improved and a test is added.
  - mlclassificator(): this class was badly handling model parameters. This is now fixed. The code has been improved.
  - various improvements of the code formatting
  - tests & examples updated
  - docstrings re-written in Google Style
  - docs updated

Additions:

  - plot_spectrum(): makes an interactive plot of a signal, and also possibly of added baselines and smoothed signals.
  - whittaker(): it is now possible to pass weights
  - baseline(): Whittaker and Gaussian process baselines going through regions of interest is also available in the baseline() function 
  - ruby_T_corr(): Computes the temperature correction for the Ruby pressure scale.
  - borate_T_corr(): Computes the temperature correction for the SrB4O7:Sm2+ pressure scale.
  - pressure_sensor(): Converts the wavelength of fluorescence lines into pressure in GPa.

Deprecation:

  - get_portion_interest() becomes extract_signal(). get_portion_interest() still works but it will be removed in a future release

0.5.3
-----

  - Correction in despiking function (thanks Kevin Yuan)

0.5.2
-----

  - fix installation : something went bad as we switched to pyproject.toml in 0.5.1

0.5.1
-----

  - fix installation dependencies

0.5.0
-----

  - map() function is now maps() to avoid any conflict with built-in Python map() function
  - add the possiblity of importing lines of spectra for Horiba spectrometers
  - add a despiking() function
  - add external calibration in rameau
  - add the posibility to calibrate K coefficient(s) on another dataset of spectra
  - correction of a typo in mlclassificator function

0.4.9
-----

- quick fix of a bug in read_horiba function (from rampy.maps) 

0.4.8
-----

  - area and area ratio calculations available for maps
  - read_horiba function (from rampy.maps) corrected (Github issue #25)
  - rubberband baseline was fixed, test added, fit_baseline example update (Github issue #1, thanks @sjfraser05)

0.4.7
-----

  - urgent correction of a bug in rp.pseudovoigt() > float entries were resulting in an error message...
  - map treatment available for HORIBA and RENISHAW spectrometer (see example folder). Please report any bug!

0.4.6
-----

  - doc improvements
  - cvxpy 1.1 or higher is an optional dependency.
  - switch to setup.cfg
  - add shiftsp() function

0.4.5
-----

  - update docs
  - enhancing the pseudovoigt section for its use with arrays
  - adding an experimental "peakarea" function
  - drPLS algorithm added in baselines (thanks @Snijderfrey)

0.4.4
-----

  - Cleaning the code
  - Correction of a bug in tlcorrection
  - add lorentzian() and pearson7() peak shape functions
  - add mlexplorator function for easy PCA/NMF on spectroscopic data
  - update docs

0.4.3
-----

  - Improvements in documentation of mlregressor.
  - Improvements of rampy.normalise and rampy.centroid. Those functions can treat arrays of spectra now.
  - Correction of a bug in rampy.normalize that caused the "area" method to not work when entering x.
  - Better tests

0.4.2
-----

  - Removing dependency to cvxpy that does not build well in Windows... It affects the use of the rampy.mixing function

0.4.1
-----

  - Correction of an error in `mlregressor` which made impossible to import X_test datasets.

0.4.0
-----

  - BREAKING CHANGE: `mlregressor` is now a class and not a function anymore. You can provide directly sklearn arguments to the algorithms through dictionaries.
    The use of the class simplifies the use of `mlregressor`, as the created objects saves everything!
    It also makes it very easy to change the algorithm and try something else. See the example in the example folder!
  - addition of the `centroid()` function, that calculates the centroid of a signal.
  - addition of tests and examples for the `mlregressor()` class, the `resample()` and `flipsp()` functions.
  - `chemical_splitting()` allows one to select the random seed.
  - Correction of the `mixing_sp()` function, rampy is now compatible with cvxpy v1.0.
  - arguments can be provided to `resample()` to use different techniques of interpolation in `scipy.interpolate.interp1d`.
  - Various documentation improvements

0.3.6
-----

  - Correction of the `tlcorrection()` function: the 'hehlen' correction was missing a frequency term to be complete (eq. 2 and 3 in Hehlen 2010 J. Phys. Condes. Matter 22: 025401).

0.3.5
-----

  - Addition of the `rampy.mixing_sp()` function. See `help(rampy.mixing_sp())`, as well as the example folder.

0.3.4
-----

  - gcvspline is not a requirement anymore. Error messages will outputs when trying to use it, inviting to install it manually. This is implemented to avoid problems with FORTRAN compilation for people not interested in using gcvspline.
  - Add early stopping in` mlregressor` neural networks.

0.3.3
-----

  - Minor dependency correction

0.3.2
-----

  - Adding the names in the `rameau` object
  - Improvements of the `mlregressor` function, with addition of neural nets and bagging neural nets algorithms

0.3.1
-----

  - Rameau is now an object-oriented interface
  - `smooth()` function updated; 10 algorithms are available.
  - updated example of peak fitting

0.3.0
-----

  - Documentation improvements
  - Python 3 compatible
  - Addition of the Rameau function

0.2.9
-----

  - addition of gaussian baseline
  - addition of `flipsp()` to flip spectra along the axis 0
  - addition of `resample()` to resample spectra with scipy.interpolate
  - improvement of `tlcorrection`

0.2.8
-----

  - addition of the `tlcorrection()` function to replace the Long function

0.2.7
-----

  - Minor correction of the `baseline()` documentation and removing a `print()` command.

0.2.6
-----

  - standardizatin in the `baseline()` function is now included (improve polynomial fits);
  - addition of the arPLS algorithm from Baek et al. (2015) to automatic fit the baseline;
  - addition of the whittaker smoother to fit the baseline (Eiler 2003);
  - addition of the ALS algorithm (Eilers and Boelens 2005).
