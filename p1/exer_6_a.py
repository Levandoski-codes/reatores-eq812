#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 20:56:38 2020

@author: levandoski
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def bat(t, x):
    Ca0 = 2
    T = 300
    k = 353105*np.exp(-5050/T)
    print(k)
    dtdx = (k*Ca0*(1-x)**2)**(-1)
    
    return dtdx


x = np.linspace(0, 0.95, 100)
t0 = 0

sol = odeint(bat, t0, x)

print(max(sol))

plt.plot(sol, x)