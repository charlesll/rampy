Signal parameters
========================

I will enrich rampy with new functions to calculate the characteristics/parameters of signals. This is not urgent as scipy.signal allows us to do many things already (e.g. with find_peaks...), and there is no need to reinvent the wheel. One function is often useful and that is not in scipy.signal:

.. autofunction:: rampy.centroid

