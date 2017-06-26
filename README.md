#RamPy
=====

This is the skeleton of a module to treat spectroscopic data in Python.

The core is in the spectratool.py module, with various functions to help organising the data, subtracting baseline (polynomial, linear, cubic splines), with access to linear, polynomial or gaussian/lorentzian function that can be used with peak fitting routine for instance. 

the /examples/ folder contain various examples.

Copyright (2015-2017) C. Le Losq.

charles.lelosq@anu.edu.au

# DISCLAIMER

Rampy contains a set of functions that I developped for my personal research use. 
I currently am focusing on the development of Spectra.jl in Julia for various spectroscopic treatment, 
but I still wanted to provide various functions I developped in Python in a clean module. 

Therefore, rampy is far less complete than Spectra.jl such that I will advocate to try the later package in Julia. 

However, for small needs as baseline fitting, Rampy should be good enougth and easy to integrate in a Python workflow.

If you are willing to help developping it, please feel free to fork it and create pull requests!

# REQUIREMENTS

For Python:
Scipy
Numpy
Matplotlib
lmfit & aeval (http://cars9.uchicago.edu/software/python/lmfit/)
gcvvspline

# INSTALLATION

For now, just clone this repo and run 

	python setup.py

will install rampy. In a short future, you will need to run only:

	pip install rampy 

# EXAMPLES

See the example folder

Updated June 2017





