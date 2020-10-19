#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 21:44:06 2020

@author: levandoski
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def bat(y, t):
    X, T = y
           
    Ca0 = 2
    V = 1200
    Cpa = Cpb = 20
    Na0 = Ca0*V
    deltaHpad = -10000
    k = 353105*np.exp(-5050/T)
    ra = k*(Ca0**2)*((1-X)**2)
    dXdt = k*Ca0*(1-X)**2
    dTdt = (-deltaHpad * ra * V)/(Na0*(Cpa + Cpb))
        
    
    dydt = [dXdt, dTdt]
    return dydt


#y = X, T
y0 = [0, 300] 
t = np.linspace(0, 4, 1000) # h

sol = odeint(bat, y0, t)
X, T = sol[:, 0], sol[:, 1]

plt.figure(dpi=100, figsize=(9, 4))
plt.plot(t, X, 'b', label='Sem Q, X(t)')
plt.legend(loc='best')
plt.xlabel('t (min)')
plt.ylabel('Conversão (X)')
plt.title('Reator batelada X vs t')
plt.grid()
plt.show()


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(t, T - 273.15, 'b', label='Sem Q, T(t)')
plt.legend(loc='best')
plt.xlabel('t (h)')
plt.ylabel('Temperatura (ºC)')
plt.title('Reator batelada T vs t')
plt.grid()
plt.show()


for i, _ in enumerate(X):
#    print(f"{t[i]} : {X[i]}")
    if (T[i] > 190 + 273.15 and T[i] < 201 + 273.15):
        print(f"{t[i]} : {T[i] - 273.15}")