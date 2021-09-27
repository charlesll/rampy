# Rampy News

Copyright (c) 2014-2021 Dr. Charles Le Losq

email: lelosq@ipgp.fr

Licence: see LICENCE.md

# Wanted features (planned)

- Addition of more unsupervised machine learning techniques

- Classification via 1D CNN and other algorithms using Keras or Pytorch

- Addition of external calibration method in rameau

- peak fitting for maps

# 0.4.9 (stable)

- quick fix of a bug in read_horiba function (from rampy.maps) 

# 0.4.8

- area and area ratio calculations available for maps

- read_horiba function (from rampy.maps) corrected (Github issue #25)

- rubberband baseline was fixed, test added, fit_baseline example update (Github issue #1, thanks @sjfraser05)

# 0.4.7

- urgent correction of a bug in rp.pseudovoigt() > float entries were resulting in an error message...

- map treatment available for HORIBA and RENISHAW spectrometer (see example folder). Please report any bug!

# 0.4.6

- doc improvements

- cvxpy 1.1 or higher is an optional dependency.

- switch to setup.cfg

- add shiftsp() function

# 0.4.5

- update docs

- enhancing the pseudovoigt section for its use with arrays

- adding an experimental "peakarea" function

- drPLS algorithm added in baselines (thanks @Snijderfrey)

# 0.4.4

- Cleaning the code

- Correction of a bug in tlcorrection

- add lorentzian() and pearson7() peak shape functions

- add mlexplorator function for easy PCA/NMF on spectroscopic data

- update docs

# 0.4.3

- Improvements in documentation of mlregressor.

- Improvements of rampy.normalise and rampy.centroid. Those functions can treat arrays of spectra now.

- Correction of a bug in rampy.normalize that caused the "area" method to not work when entering x.

- Better tests

# 0.4.2

- Removing dependency to cvxpy that does not build well in Windows... It affects the use of the rampy.mixing function

# 0.4.1

- Correction of an error in `mlregressor` which made impossible to import X_test datasets.

# 0.4.0

- BREAKING CHANGE: `mlregressor` is now a class and not a function anymore. You can provide directly sklearn arguments to the algorithms through dictionaries.
  The use of the class simplifies the use of `mlregressor`, as the created objects saves everything!
  It also makes it very easy to change the algorithm and try something else. See the example in the example folder!

- addition of the `centroid()` function, that calculates the centroid of a signal.

- addition of tests and examples for the `mlregressor()` class, the `resample()` and `flipsp()` functions.

- `chemical_splitting()` allows one to select the random seed.

- Correction of the `mixing_sp()` function, rampy is now compatible with cvxpy v1.0.

- arguments can be provided to `resample()` to use different techniques of interpolation in `scipy.interpolate.interp1d`.

- Various documentation improvements

# 0.3.6

- Correction of the `tlcorrection()` function: the 'hehlen' correction was missing a frequency term to be complete (eq. 2 and 3 in Hehlen 2010 J. Phys. Condes. Matter 22: 025401).

# 0.3.5

- Addition of the `rampy.mixing_sp()` function. See `help(rampy.mixing_sp())`, as well as the example folder.

# 0.3.4

- gcvspline is not a requirement anymore. Error messages will outputs when trying to use it, inviting to install it manually. This is implemented to avoid problems with FORTRAN compilation for people not interested in using gcvspline.

- Add early stopping in` mlregressor` neural networks.

# 0.3.3

- Minor dependency correction

# 0.3.2

- Adding the names in the `rameau` object
- Improvements of the `mlregressor` function, with addition of neural nets and bagging neural nets algorithms

# 0.3.1

- Rameau is now an object-oriented interface

- `smooth()` function updated; 10 algorithms are available.

- updated example of peak fitting

# 0.3.0

- Documentation improvements

- Python 3 compatible

- Addition of the Rameau function

# 0.2.9

- addition of gaussian baseline

- addition of `flipsp()` to flip spectra along the axis 0

- addition of `resample()` to resample spectra with scipy.interpolate

- improvement of `tlcorrection`

# 0.2.8

- addition of the `tlcorrection()` function to replace the Long function

# 0.2.7

- Minor correction of the `baseline()` documentation and removing a `print()` command.

# 0.2.6

- standardizatin in the `baseline()` function is now included (improve polynomial fits);

- addition of the arPLS algorithm from Baek et al. (2015) to automatic fit the baseline;

- addition of the whittaker smoother to fit the baseline (Eiler 2003);

- addition of the ALS algorithm (Eilers and Boelens 2005).
