#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 13:40:20 2020

@author: levandoski
"""
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def bat(y, t):
    X, T = y
    
    T1 = 163 + 273.15# C
    Ea = 28960 # cal/g mol
    deltaHpad = -83 #cal/g
    Cpa = 0.5 #cal/g C
    k1 = 0.8 # 1/h
    R = 1.987 # cal /mol C
    k = k1*np.exp((Ea/R)*(1/T1 - 1/T))
    
    dTdt = -deltaHpad*k*(1-X)/Cpa
    dXdt = k*(1-X)
    
    dydt = [dXdt, dTdt]
    return dydt


#y = X, T
y0 = [0, 163 + 273.15] 
t = np.linspace(0, 0.1172, 1000) # h

sol = odeint(bat, y0, t)
X, T = sol[:, 0], sol[:, 1]

plt.figure(dpi=100, figsize=(9, 4))
plt.plot(t, X, 'b', label='Sem Q, X(t)')
plt.legend(loc='best')
plt.xlabel('t (h)')
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
    if (X[i] > 0.96 and X[i] < 1):
        print(f"{t[i]} : {X[i]}")