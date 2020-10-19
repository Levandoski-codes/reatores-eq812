#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 10:52:30 2020

@author: levandoski
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def PFR(y, V):
    T, X = y
    
    v0 = 0.1 #m³/s
    
    Cpa = 85 # J/(kg K)
    Cpa = Cpa * 250 # J/(mol K)
    deltaHpad = -18590 #J/mol
    
    k = 2485*np.exp(-1665/T) / 60 #1/h ==> 1/s
    
    dTdv = (-k*(1-X)*deltaHpad)/(v0*Cpa)
    dXdv = (k*(1-X))/(v0)
    
    dydV = [dTdv, dXdv]
    return dydV

y0 = [361, 0]
V = np.linspace(0, 1, 1000) # h

sol = odeint(PFR, y0, V)
T, X = sol[:, 0], sol[:, 1]

plt.figure(dpi=100, figsize=(9, 4))
plt.plot(V, X, 'b')
plt.xlabel('V (m³)')
plt.ylabel('Conversão (X)')
plt.title('Reator batelada X vs V')
plt.grid()
plt.show()