# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import numpy as np
import cvxpy

def mixing_sp(y_fit: np.ndarray, ref1: np.ndarray, ref2: np.ndarray) -> np.ndarray:
    """
    Mixes two reference spectra to match given experimental signals.

    This function calculates the fractions of the first reference spectrum (`ref1`) 
    in a linear combination of `ref1` and `ref2` that best matches the provided signals (`y_fit`). 
    The calculation minimizes the sum of the least absolute values of the objective function:
    \( \text{obj} = \sum \left| y_{\text{fit}} - (F_1 \cdot \text{ref1} + (1 - F_1) \cdot \text{ref2}) \right| \).

    Args:
        y_fit (np.ndarray): Array containing the experimental signals with shape `(m, n)`, 
            where `m` is the number of data points and `n` is the number of experiments.
        ref1 (np.ndarray): Array containing the first reference signal with shape `(m,)`.
        ref2 (np.ndarray): Array containing the second reference signal with shape `(m,)`.

    Returns:
        np.ndarray: Array of shape `(n,)` containing the fractions of `ref1` in the mix. 
        Values range between 0 and 1.

    Notes:
        - The calculation is performed using `cvxpy` for optimization.
        - Ensure that `y_fit`, `ref1`, and `ref2` have compatible dimensions.

    Example:

        >>> import numpy as np
        >>> y_fit = np.array([[0.5, 0.6], [0.4, 0.5], [0.3, 0.4]])
        >>> ref1 = np.array([0.5, 0.4, 0.3])
        >>> ref2 = np.array([0.2, 0.3, 0.4])
        >>> fractions = mixing_sp(y_fit, ref1, ref2)
        >>> print(fractions)

    """

    ref1 = ref1.reshape(1,-1)
    ref2 = ref2.reshape(1,-1)

    F1 = cvxpy.Variable(shape=(y_fit.shape[1],1))

    objective = cvxpy.Minimize(cvxpy.sum(cvxpy.abs(cvxpy.multiply(F1,ref1) + cvxpy.multiply((1-F1),ref2) - y_fit.T))) 

    constraints = [0 <= F1, F1 <= 1]

    prob = cvxpy.Problem(objective, constraints)
    prob.solve()
    return np.asarray(F1.value).reshape(-1)
