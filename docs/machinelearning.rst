Machine learning
================

Rampy offers three classes for performing classification, regression or unsupervised ML exploration of a set of spectra. Those are helpful ways to automatically perform usual scalling, test-train splitting and ML training of popular ML algorithms, and use `scikit-learn <https://scikit-learn.org/stable/>`_ in the background.

Those functions work (``rp.regressor`` was used in `this publication <https://doi.org/10.2138/am-2019-6887>`_) but still may evolve in the futur. For advanced, ML, I suggest using directly scikit-learn or other ML libraries.

Below you will find the documentation of the relevant functions, and have a look at the example notebooks too: :ref:`RST Notebooks`

Do not hesitate to ask for new features depending on your needs !

Machine learning classification
-------------------------------

Based on a set of spectra and their labels, the ``rampy.ml_classification`` module allows you to perform a classification of the spectra using a supervised ML algorithm. The class will take care of splitting the data into training and test sets, scaling the data, and training the model. You can then use the trained model to predict the labels of new spectra.

.. automodule:: rampy.ml_classification
   :members:
   :undoc-members:
   :show-inheritance:

Machine learning exploration
----------------------------

The ``rampy.ml_exploration`` module allows you to perform unsupervised ML exploration of a set of spectra. The class will take care of scaling the data and training the model. You can then use the trained model to explore the data and find patterns.

.. automodule:: rampy.ml_exploration
   :members:
   :undoc-members:
   :show-inheritance:

Machine learning regression
---------------------------


Based on a set of spectra and their labels, the ``rampy.ml_regressor`` module allows you to perform a regression using the spectra and a supervised ML algorithm. The class will take care of splitting the data into training and test sets, scaling the data, and training the model. You can then use the trained model to predict the new values of your target from new spectra.

.. automodule:: rampy.ml_regressor
   :members:
   :undoc-members:
   :show-inheritance:

Linear mixture
--------------

This function helps you solve a simple problem: you have spectra that are obtained by a linear combination of two endmember spectra.

If you have the two endmember spectra, you can use the ``rampy.mixing()`` function to know the fraction of each endmember in the mixture.

If you do not know the endmember spectra, then you may be interested in using directly the PyMCR library, see the documentation `here <https://pages.nist.gov/pyMCR/>`_ and an example notebook `here <https://github.com/usnistgov/pyMCR/blob/master/Examples/Demo.ipynb>`_. We used it in `this publication <https://doi.org/10.2138/am-2019-6887>`_, see the code `here <https://github.com/charlesll/rampy/blob/master/examples/Iron_AmMin_paper/Iron_MORB_code.ipynb>`_.

.. autoclass:: rampy.mixing.mixing_sp
   :members:
   :undoc-members:
   :show-inheritance: