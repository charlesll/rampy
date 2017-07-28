RamPy
=====

This is the skeleton of a module to treat spectroscopic data in Python.

The baseline function implements various interesting algorithms for baseline subtraction.

The spectratool.py part of the library offers various functions to help organising the data, subtracting baseline (polynomial, linear, cubic splines), with access to linear, polynomial or gaussian/lorentzian function that can be used with peak fitting routine for instance. 

the /examples/ folder contain various examples.

Copyright (2015-2017) C. Le Losq.

charles.lelosq@anu.edu.au

# DISCLAIMER

Rampy contains a set of functions that I developped for my personal research use. 
I currently am focusing on the development of Spectra.jl in Julia for various spectroscopic treatment, 
but I still wanted to provide various functions I developped in Python in a clean module. 

Rampy is less compelte than Spectra.jl but I aim at maintaining the baseline function efficiently, which may become the backend of the baseline function in Spectra.jl through a PyCall instance.

Therefore, part of rampy are not complete the the baseline function and others definitely work. For more functionalities, I invite to check Spectra.jl in Julia, or to push commits to add functionalities to RamPy!


# REQUIREMENTS

For Python:
Scipy
Numpy
Matplotlib
lmfit & aeval (http://cars9.uchicago.edu/software/python/lmfit/)
gcvvspline (you need a working FORTRAN compiler for its installation... Warning Windows users! Check you FORTRAN compiler!)

# INSTALLATION

Install with pip:

	pip install rampy 

# EXAMPLES

See the example folder

Updated July 2017





