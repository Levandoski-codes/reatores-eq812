#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 00:32:39 2020

@author: levandoski
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def exponential_decay(t, y): return -0.5 * y

t_e = np.linspace(0, 10, 10)
t_e = None

sol = solve_ivp(exponential_decay, [0, 10], (4, ), t_eval=t_e)

t = sol.t
print(sol.t)
y = sol.y[0, :]
#print(sol.y)

plt.plot(t, y)
plt.show()