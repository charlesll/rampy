# Rampy News

Copyright (c) 2014-2018 Dr. Charles Le Losq

email: charles.lelosq@anu.edu.au

Licence: see LICENCE.md

As Rampy starts to grow, I will summarise changes in this file starting at version 0.2.6

# 0.4

- BREAKING CHANGE: mlregressor is now a class and not a function anymore. You can provide directly sklearn arguments to the algorithms through dictionnaries.
The use of the class simplifies the use of mlregressor, as the created objects saves everything!
It also makes it very easy to change the algorithm and try something else.

- rampy.chemical_splitting() allows one to select the random seed.

- addition of tests and examples of the mlregressor() class and of resample() and flipsp() functions.

- Correction of the rp.mixing_sp() function, rampy is now compatible with cvxpy v1.0.

- arguments can be provided to resample() to use different techniques of interpolation in scipy.interpolate.interp1d.

- addition of the centroid() function, that calculates the centroid of a signal.

# 0.3.6

- Correction of the tlcorrection() function: the 'hehlen' correction was missing a frequency term to be complete (eq. 2 and 3 in Hehlen 2010 J. Phys. Condes. Matter 22: 025401).

# 0.3.5

- Addition of the rampy.mixing_sp() function. See help(rampy.mixing_sp()), as well as the example folder.

# 0.3.4

- gcvspline is not a requirement anymore. Error messages will outputs when trying to use it, inviting to install it manually. This is implemented to avoid problems with FORTRAN compilation for people not interested in using gcvspline.

- Add early stopping in mlregressor neural networks.

# 0.3.3

- Minor dependency correction

# 0.3.2

- Adding the names in the rameau object
- Improvements of the mlregressor function, with addition of neural nets and bagging neural nets algorithms

# 0.3.1

- Rameau is now an object-oriented interface

- smooth() function updated; 10 algorithms are available.

- updated example of peak fitting

# 0.3.0

- Documentation improvements

- Python 3 compatible

- Addition of the Rameau function

# 0.2.9

- addition of gaussian baseline

- addition of flipsp() to flip spectra along the axis 0

- addition of resample() to resample spectra with scipy.interpolate

- improvement of tlcorrection

# 0.2.8

- addition of the tlcorrection() function to replace the Long function

# 0.2.7

- Minor correction of the baseline() documentation and removing a print() command.

# 0.2.6

- standardizatin in the baseline() function is now included (improve polynomial fits);

- addition of the arPLS algorithm from Baek et al. (2015) to automatic fit the baseline;

- addition of the whittaker smoother to fit the baseline (Eiler 2003);

- addition of the ALS algorithm (Eilers and Boelens 2005).
