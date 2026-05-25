# -*- coding: utf-8 -*-
"""

"""
import numpy as np
from scipy.interpolate import interp1d

lmbd = [8.25, 11.6, 28.9, 52.5, 61.5]
exponent = [-1.0, 0.0, 1.0, 2.0, 3.0]


lmbd = np.array(lmbd)
exponent = np.array(exponent)

f = interp1d(lmbd, exponent)

def tlak_perlit(lmbd):

    exp = f(lmbd*1000)
    p = 10**(exp)
    return p #Pa


tlak_perlit(0.04)
