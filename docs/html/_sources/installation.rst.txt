Installation
============

General preparation
-------------------

Rampy runs with a traditional Python stack.

If you are not familiar with Python, you can first have a look at the `scipy lecture notes <https://scipy-lectures.org/>`_,
a set of tutorials for the beginner.

You can install `Anaconda Python <https://www.anaconda.com/products/individual>`_ to get a running Python distribution. See the documentation of Anaconda for those steps.

Rampy installation
------------------

Install with pip in the command line:

 ``pip install rampy``

Optional dependencies
---------------------

If you want to use the "gcvspline" option from the ``baseline()`` function, the ``mixing()`` function, or the ``rampy.rameau()`` class. please
install the following libraries:

 ``pip install gcvspline``

 ``pip install cvxpy``

Please note that gcvspline installation requires gfortran, see `gcvspline documentation <https://charlesll.github.io/gcvspline/>`_.
For Windows, wheels for Python 2.7, 3.4 and 3.6 are provided for 64 bit Windows, and a wheel for Python 3.6 is provided for Windows 32 bits. If installation fails, please check if is due to a fortran compiler issue.

As for gcvspline, the installation of cvxpy can cause problems for Windows users due to missing compiler. See instructions from cvxpy in this case.

Additional libraries for model fitting may be wanted:

- lmfit & aeval (http://cars9.uchicago.edu/software/python/lmfit/)
- emcee
