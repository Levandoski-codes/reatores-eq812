#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 13:40:20 2020

@author: levandoski
"""
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def PFR(y, t):
    V = y
    X = t
    
    k = 0.8
    m = 133700
    rho = 0.9
    vo = m/(rho*1000)
    dVdX = vo/(k*(1-X))
    
    dydt = dVdX
    return dydt

y0 = 0
x = np.linspace(0, 0.97, 100) # h

sol = odeint(PFR, y0, x)
V = sol

plt.figure(dpi=100, figsize=(9, 4))
plt.plot(V, x, 'b', label='Sem Q, X(t)')
plt.legend(loc='best')
plt.xlabel('V (L)')
plt.ylabel('Conversão (X)')
plt.title('Reator batelada X vs V')
plt.grid()
plt.show()