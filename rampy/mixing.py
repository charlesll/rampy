import cvxpy
import numpy as np

def mixing_sp(y_fit,ref1,ref2):
    """mix two reference spectra to match the given ones

    Parameters
    ----------
    y_fit : ndarray, shape m * n
        an array containing the signals with m datapoints and n experiments
    ref1 : ndarray, shape m
        an array containing the first reference signal
    ref2 : ndarray, shape m
        an array containing the second reference signal

    Returns
    -------
    out : ndarray, shape n
        the fractions of ref1 in the mix

    Notes
    -----
    Performs the calculation by minimizing the sum of the least absolute value of the objective function:
        obj = sum(abs(y_fit-(ref1*F1 + ref2*(1-F1))))

    Uses cvxpy to perform this calculation
    """

    ref1 = ref1.reshape(1,-1)
    ref2 = ref2.reshape(1,-1)
    
    F1 = cvxpy.Variable(shape=(y_fit.shape[1],1))

    objective = cvxpy.Minimize(cvxpy.sum(cvxpy.abs(F1*ref1 + (1-F1)*ref2 - y_fit.T))) 

    constraints = [0 <= F1, F1 <= 1]

    prob = cvxpy.Problem(objective, constraints)
    prob.solve()
    return np.asarray(F1.value).reshape(-1)
