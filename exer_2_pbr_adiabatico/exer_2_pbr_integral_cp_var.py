

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 00:32:39 2020

@author: levandoski
"""

import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt


def function(X, T, T0, P):
    TR = 298 # K    
    
    FE0 = 945894.82 # mol/h
    FI = 2417286.76 # mol/h
    
#    yE0 = FE0 / (FE0 + FI)
    
    FE = FE0*(1 - X)
    FH = FE0*X
    FA = FE0*X
    FT = (FE + FI + FH + FA)
    
    yE = FE/FT
    yH = FH/FT
    yA = FA/FT
    
    PE = yE*P # atm
    PH = yH*P # atm
    PA = yA*P # atm
    
    R = 1.987 # cal/mol*K
    R_u = 8.314 #J/mol*K
    
    dGA = -133000 # J/mol
    dGE = -167900 # J/mol
    dGH = 0 # J/mol
    
    dG_padrao = dGA + dGH - dGE # J/mol
    
    dHrx = 82500 # J/mol*K    
    
    CpA = 8.314*(1.693+0.017978*T - 0.000006158*T**2)
    CpE = 8.314*(3.518+0.020001*T - 0.000006002*T**2)
    CpH = 8.314*(3.249+0.000422*T + (0.083*10**5 / (T**2)))
    CpI = 8.314*(3.47+0.00145*T + (0.121*10**5 / (T**2)))
    
    dCp = CpA + CpH - CpE
    
    T_num = X * (-dHrx) + (CpE + (FI/FE0)*CpI)*T0 + X*TR*dCp
    T_den =  (CpE + (FI/FE0)*CpI) + X*dCp
    T = T_num/T_den    
    
    K = np.exp(-dG_padrao/(R_u*T))
    
    KE = np.exp(5560/(R*T) - 5.97) # atm^-1
    KA = np.exp(11070/(R*T) - 9.40) # atm^-1
    KH = np.exp(6850/(R*T) - 7.18) # atm^-1
    k = np.exp(-15300/(R*T) + 15.32) #mol/ g*h
    
    rE = k*KE*(PE - (PA*PH)/K) / (1 + KE*PE + KA*PA + KH*PH)**2
    
    
#    print("T {}\n X {}\n FE {}\n FI {}\n FA {}\n rE {}".format(T, X, FE, FI, FA, rE))
#    print("CpA {}\n CpH {}\n CpE {}\n CpI {}\n dCp {}\n".format(CpA, CpH, CpE, CpI, dCp))
    
    return FE0/-rE, T, -rE


x = np.linspace(0, 0.75, 25)
T = T0 = 750 # K
P = 50 # atm

sol = []
T_axis = []
re_axis = []
reac_axis = []

for i in x:
    reac, T, re = function(i, T, T0, P)
    reac_axis.append(reac)
    T_axis.append(T)
    re_axis.append(re)

print(f"Massa de Catalisador sem P_var {simps(reac_axis, x)/1000:5.7} kg")

#
#fig = plt.figure(dpi=100)
#plt.title("Gráfico de Levenspiel")
#plt.xlabel("X (Conversão)")
#plt.ylabel(r"$\frac{F_{E0}}{-r_E}$")
#plt.plot(x, reac_axis, marker=".")
#plt.show()
#
#fig = plt.figure(dpi=100)
#plt.title("Temperatura")
#plt.xlabel("X (Conversão)")
#plt.ylabel("Temperatura (ºC)")
#plt.plot(x, T_axis)
#plt.show()
#
#
#fig = plt.figure(dpi=100)
#plt.title("Taxa de reação")
#plt.xlabel("X (Conversão)")
#plt.ylabel("$r_e$")
#plt.plot(x, re_axis)
#plt.show()
#
#
#fig = plt.figure(dpi=100)
#plt.title("Taxa de reação")
#plt.xlabel("T (Conversão)")
#plt.ylabel("$r_e$")
#plt.plot(T_axis, re_axis)
#plt.show()
