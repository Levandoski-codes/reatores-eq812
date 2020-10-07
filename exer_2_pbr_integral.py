#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 00:32:39 2020

@author: levandoski
"""

import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt


def function(X): 
    P = 20 # atm
    T0 = 750 #500 # K
    TR = 298 # K
    
    FE0 = 945894.82 # mol/h
    FI = 2417286.76 # mol/h
    
#    FE0 = 904493.1626 # mol/h
#    FI = 1156647.238 * 2# mol/h
    
#    yE0 = FE0 / (FE0 + FI)
    
    FE = FE0*(1 - X)
    FH = FE0*X
    FA = FE0*X
    FT = (FE + FI + FH + FA)
    
    yE = FE/FT
    yH = FH/FT
    yA = FA/FT
    
    PE = yE*P
    PH = yH*P
    PA = yA*P
    
    R = 1.987 # cal/mol*K
    R_u = 8.314 #J/mol*K
#    R = R_u
    
    dGA = -133000
    dGE = -167900
    dGH = 0
    
    dG_padrao = dGA + dGH - dGE
#    dG_padrao = 34900
    
    dHrx = 82500 # J/mol*K
    #Ref: https://webbook.nist.gov/chemistry/name-ser/
#    CpE = 95.78 # J/mol*K a 500K
#    CpH = 29.26 # J/mol*K
#    CpA = 76.68 # J/mol*K 
#    CpI = 35.22 # J/mol*K aproximação
    
    CpE = 118.83 # J/mol*K a 750
    CpH = 29.50 # J/mol*K
    CpA = 96.04 # J/mol*K 
    CpI = 38.00 # J/mol*K aproximação

    dCp = CpA + CpH - CpE
    
    T_num = X * (-dHrx) + (CpE + (FI/FE0)*CpI)*T0 + X*TR*dCp
    T_den =  (CpE + (FI/FE0)*CpI) + X*dCp
    T = T_num/T_den    
    
    K = np.exp(-dG_padrao/(R_u*T))
#    K = 2.26e-4
    KE = np.exp(5560/(R*T) - 5.97) # atm^-1
    KA = np.exp(11070/(R*T) - 9.40) # atm^-1
    KH = np.exp(6850/(R*T) - 7.18) # atm^-1
    k = np.exp(-15300/(R*T) + 15.32) #mol/ g*h
    
    rE = k*KE*(PE - (PA*PH)/K) / (1 + KE*PE + KA*PA + KH*PH)**2
    
#    dXdW = rE/FE0
    
    print("T {}\n X {}\n FE {}\n FI {}\n FA {}\n rE {}".format(T, X, FE, FI, FA, rE))
    
    return FE0/-rE


x = np.linspace(0, 0.75, 20)

#w_eval = None
#w = [0, 1e4] # g

sol = function(x)

plt.xlabel("x")
plt.ylabel("FE0/-rE")
plt.plot(x, sol)

print(trapz(sol, x))

#CpA = 8.314(1.693+0.017978*T - 0.000006158*T**2)
#CpE = 8.314(3.518+0.020001*T - 0.000006002*T**2)
#CpH = 8.314(3.249+0.000422*T + (0.083*10**5 / (T**2)))
#CpI = 8.314(3.47+0.00145*T + (0.121*10**5 / (T**2)))