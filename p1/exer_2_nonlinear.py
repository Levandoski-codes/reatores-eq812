#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 18:46:06 2020

@author: levandoski
"""

from scipy.optimize import fsolve
import numpy as np

def func(x):
    rho = 0.9
    x3 = 0.97
    k = 0.8
    V, x1, x2 = x
    ma = 133700
    v0 = ma/rho
    f = [v0*(x1)/(k*(1-x1)) - V,
         v0*(x2 - x1)/(k*(1-x2)) - V,
         v0*(x3 - x2)/(k*(1-x3)) - V]
    
    
    return f

root = fsolve(func, [1, 0, 0]) # chute inicial
print(root)

print(np.isclose(func(root), [0.0, 0.0, 0.0], atol=1e-5)) # usa o resultado encontrado pra verificar consegue zerar
