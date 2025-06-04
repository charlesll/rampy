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

If you want to use the "MSESmoothedNSpline" or "DOFSmoothedNSpline" options from the ``smooth()`` function, please install the following library:

 ``pip install gcvspline``

Please note that gcvspline installation requires gfortran, see `gcvspline documentation <https://charlesll.github.io/gcvspline/>`_.

Installation is done? Have a look here if you are new to Python: :ref:`RST First`