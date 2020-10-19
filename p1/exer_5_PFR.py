#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 13:40:20 2020

@author: levandoski
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def PFR(y, v):
    T, X, Ta = y
    
    v0 = 0.1 #m³/s
    T0 = 361 # K
    Pt = 5 * 101325 # atm ==> Pa
    
    Pa0 = Pt * 1 # Pa
    R = 8.314 # J /mol K
    
    Fa0 = (Pa0 * v0)/(R*T0) # mol de a / s
    
    Cpa = 85 # J/(kg K)
    deltaHpad = -18590 #J/mol
    
    T = ((-X*deltaHpad) + Cpa*T0)/Cpa
    
    k = 2485*np.exp(-1665/T) / 60 #1/h ==> 1/s
    
    Pa = Pa0*(1 - X)
    ra = k*Pa
    dVdX = Fa0/ra

    return dTdV, dXdV, dTadV

V0 = 0
X = np.linspace(0, 0.5, 1000) # h

sol = odeint(PFR, V0, X)
V = sol

plt.figure(dpi=100, figsize=(9, 4))
plt.plot(V, X, 'b')
plt.xlabel('V (m³)')
plt.ylabel('Conversão (X)')
plt.title('Reator batelada X vs V')
plt.grid()
plt.show()


