#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 15:03:08 2020

@author: levandoski
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def BM(w, X):
    
    k = 0.01
    fa0 = 10
    kc = 4242
    ca0 = 0.001
    ca = ca0*(1 - X)
    ra = -(k*kc*ca)/(k + kc)
    dydw = -ra/fa0 
    return 1/dydw


#y = X, T
w0 = 0 
x = np.linspace(0, 0.6, 100)

sol = odeint(BM, w0, x)
w = sol[:, 0]



plt.figure(dpi=100, figsize=(9, 4))
plt.plot(w, x, 'b')
plt.xlabel('Massa de catalisador')
plt.ylabel('Conversao (X)')
plt.title('Reator X vs n')
plt.grid()
# plt.savefig('convert_n.png')
plt.show()