#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 00:32:39 2020

@author: levandoski
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def function(W, X): 
    P = 20
    T0 = 500
    TR = 298
    
    FE0 = 945894.82
    FI = 2417286.76
    
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
    
    R = 1.987
    R_u = 8.314
    R = R_u
    
    dGA = -133000
    dGE = -167900
    dGH = 0
    
    dG_padrao = dGA + dGH - dGE
    dHrx = 82500
    
    CpE = 95.78
    CpH = 29.26
    CpA = 76.68
    CpI = 36.016

    dCp = CpA + CpH - CpE
    
    T_num = X * (-dHrx) + (CpE + (FI/FE0)*CpI)*T0 + X*TR*dCp
    T_den =  (CpE + (FI/FE0)*CpI) + X*dCp
    T = T_num/T_den
    
    K = np.exp(-dG_padrao/(R_u*T))
    KE = np.exp(5560/(R*T) - 5.97)
    KA = np.exp(11070/(R*T) - 9.40)
    KH = np.exp(6850/(R*T) - 7.18)
    k = np.exp(-15300/(R*T) - 15.32)
    
    rE = k*KE*(PE - (PA*PH)/K) / (1 + KE*PE + KH*PH + KA*PA)**2
    
    dXdW = -rE/FE0
    
    
    return dXdW

w = [0, 10000]
sol = solve_ivp(function, w, [0, ])

x = sol.y[0, :]
w = sol.t
plt.xlabel("x")
plt.ylabel("w")
plt.plot(x, w)
