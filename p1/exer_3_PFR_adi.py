#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 13:40:20 2020

@author: levandoski
"""
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def function(X, T, T0):
    
    
    m = 133700# g/h
    Fa0 = m/250 # g/h / g/mol  = mol/h
    
    T1 = 163 + 273.15 # K    
    Ea = 28960 # cal/g mol
    deltaHpad = -83 #cal/g
    Cpa = 0.5 #cal/g C
    k1 = 0.8 # 1/h
    R = 1.987 # cal /mol C
    
    T = ((-X * deltaHpad + Cpa*T0)/Cpa) + 273.15
    
    k = k1*np.exp((Ea/R)*(1/T1 - 1/T))
    
    rho = 0.9 * 1000 # g/cm³ kg/m³
    vo = m/rho
    vo = vo
    
    ra = k*((Fa0/vo)*(1-X))
    Ca = (Fa0/vo)*(1-X)
    
    
    return Fa0/ra, T, ra, k, Ca

xf = 0.97

x = np.linspace(0, xf, 1000)
T = T0 = 20 # C
T0 = 163 #<======================================== DEVIA SER 20

T_axis = []
re_axis = []
reac_axis = []
k_axis = []
Ca_axis = []

for i in x:
    reac, T, re, k, Ca = function(i, T, T0)
    reac_axis.append(reac)
    T_axis.append(T)
    re_axis.append(re)
    k_axis.append(k)
    Ca_axis.append(Ca)
    
    
Vol = simps(reac_axis, x)


fig = plt.figure(dpi=100)
plt.title("Gráfico de Levenspiel")
plt.xlabel("X (Conversão)")
plt.ylabel(r"$\frac{F_{A0}}{-r_a}$")
plt.plot(x, reac_axis)
plt.show()

fig = plt.figure(dpi=100)
plt.title("Temperatura")
plt.xlabel("X (Conversão)")
plt.ylabel("Temperatura (ºC)")
plt.plot(x, np.asanyarray(T_axis) - 273.15)
plt.show()

fig = plt.figure(dpi=100)
plt.title("k")
plt.xlabel("X (Conversão)")
plt.ylabel("k (1/h)")
plt.plot(x, k_axis)
plt.show()


fig = plt.figure(dpi=100)
plt.title("Taxa de reação")
plt.xlabel("X (Conversão)")
plt.ylabel("$-r_a$")
plt.plot(x, re_axis)
plt.show()


fig = plt.figure(dpi=100)
plt.title("Taxa de reação")
plt.xlabel("T (K)")
plt.ylabel("$-r_a$")
plt.plot(T_axis, re_axis)
plt.show()


fig = plt.figure(dpi=100)
plt.title("Ca")
plt.xlabel("X (Conversão)")
plt.ylabel("$C_a$")
plt.plot(x, Ca_axis)
plt.show()



def PFR(y, t):
    V = y
    X = t
    
    k = 0.8
    m = 133700
    rho = 0.9
    vo = m/(rho*1000)
    T0 = 20
    T0 = 163 #<======================================== DEVERIA SER 20
    
    T1 = 163 + 273.15 # K
    Ea = 28960 # cal/g mol
    deltaHpad = -83 #cal/g
    Cpa = 0.5 #cal/g C
    k1 = 0.8 # 1/h
    R = 1.987 # cal /mol C
    
    T = (-X * deltaHpad + Cpa*T0)/Cpa + 273.15
    
    k = k1*np.exp((Ea/R)*(1/T1 - 1/T))
    
    
    dVdX = vo/(k*(1-X))
    
    dydt = dVdX
    return dydt

y0 = 0
x = np.linspace(0, xf, 1000) # h

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




