#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 10:52:30 2020

@author: levandoski
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def PFR(V, X):
    
    v0 = 0.1 #m³/s
    T0 = 361 # K
    Pt = 5 * 101325 # atm ==> Pa
    
    Pa0 = Pt * 1 # Pa
    R = 8.314 # J /mol K
    
    Fa0 = (Pa0 * v0)/(R*T0) # mol de a / s
    
    Cpa = 85 # J/(kg K)
    Cpa = Cpa * 250 # J/(mol K)
    deltaHpad = -18590 #J/mol
    
    T = ((-X*deltaHpad) + Cpa*T0)/Cpa
    
    k = 2485*np.exp(-1665/T) / 60 #1/h ==> 1/s
    
    Pa = Pa0*(1 - X)
    ra = k*Pa
    dVdX = Fa0/ra

    return dVdX

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



#from scipy.integrate import simps
#
#def function(X, T, T0):
#    
#    
#    v0 = 0.1 #m³/s
#    
#    Pt = 5 * 101325 # atm ==> Pa
#    
#    Pa0 = Pt * 1 # Pa
#    R = 8.314 # J /mol K
#    
#    Fa0 = (Pa0 * v0)/(R*T0) # mol de a / s
#    
#    Cpa = 85 # J/(kg K)
#    deltaHpad = -18590 #J/mol
#    
#    T = ((-X*deltaHpad) + Cpa*T0)/Cpa
#    
#    k = 2485*np.exp(-1665/T) / 60 #1/h ==> 1/s
#    
#    Pa = Pa0*(1 - X)
#    ra = k*Pa
#    Pb = Pa0*X
#    
#    
#    return Fa0/ra, T, ra, k, Pa, Pb
#
#x = np.linspace(0, 0.5, 1000)
#T = T0 = 361 
#
#
#T_axis = []
#re_axis = []
#reac_axis = []
#k_axis = []
#Pa_axis = []
#Pb_axis = []
#
#for i in x:
#    reac, T, re, k, Pa, Pb = function(i, T, T0)
#    reac_axis.append(reac)
#    T_axis.append(T)
#    re_axis.append(re)
#    k_axis.append(k)
#    Pa_axis.append(Pa)
#    Pb_axis.append(Pb)
#    
#    
#Vol = simps(reac_axis, x)
#
#
#fig = plt.figure(dpi=100)
#plt.title("Gráfico de Levenspiel")
#plt.xlabel("X (Conversão)")
#plt.ylabel(r"$\frac{F_{A0}}{-r_a}$")
#plt.plot(x, reac_axis)
#plt.show()
#
#fig = plt.figure(dpi=100)
#plt.title("Temperatura")
#plt.xlabel("X (Conversão)")
#plt.ylabel("Temperatura (ºC)")
#plt.plot(x, np.asanyarray(T_axis))
#plt.show()
#
#fig = plt.figure(dpi=100)
#plt.title("k")
#plt.xlabel("X (Conversão)")
#plt.ylabel("k (1/h)")
#plt.plot(x, k_axis)
#plt.show()
#
#fig = plt.figure(dpi=100)
#plt.title("Taxa de reação")
#plt.xlabel("X (Conversão)")
#plt.ylabel("$-r_a$")
#plt.plot(x, re_axis)
#plt.show()
#
#fig = plt.figure(dpi=100)
#plt.title("Taxa de reação")
#plt.xlabel("T (K)")
#plt.ylabel("$-r_a$")
#plt.plot(T_axis, re_axis)
#plt.show()
#
#fig = plt.figure(dpi=100)
#plt.title("Pressão")
#plt.xlabel("V (m³)")
#plt.ylabel("$P_a \ \ (Pa)$")
#plt.plot(V, Pa_axis, 'r')
#plt.plot(V, Pb_axis, 'b')
#plt.show()