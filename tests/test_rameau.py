# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.spatial import ConvexHull
import scipy.sparse as sparse
from numpy.linalg import norm
from sklearn import preprocessing

import rampy as rp

