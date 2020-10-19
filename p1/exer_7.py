#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 23:06:28 2020

@author: levandoski
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint

global a_CO
global a_CM
global a_CP

def pbr(y, w):
    global a_CO
    global a_CM
    global a_CP
    T, FM, FO, FP = y
    
    CT0 = 2 #mol/dm3
    deltaH1 = -1800 #J/mol o-xileno
    deltaH3 = -1100 #J/mol o-xileno
    Kc=10*np.exp(4.8*(430/T-1.5))
    k1 = 0.5*np.exp(2*(1-320/T)) #dm3/kgcat.min
    k2 = k1/Kc
    k3=0.005*np.exp((4.6*(1-(460/T)))) #dm3/kgcat.min
    T0 = 330 #K
    Ta=500 #K
    Uarho = 16 #J/kg cal.min.°C
    
    FT = FO + FM + FP
    
    CO = CT0*(FO/FT)*(T0/T)
    CM = CT0*(FM/FT)*(T0/T)
    CP = CT0*(FP/FT)*(T0/T)
    a_CO.append(CO)
    a_CM.append(CM)
    a_CP.append(CP)
    
    rM = k1*CO - k2*CM
    rP = k3*CO
    rO = -(rM + rP)
    
    CpO = 100 #J/mol.K
    CpM = 100 #J/mol.K
    CpP = 100 #J/mol.K
    
    soma1 = -deltaH1*rM - deltaH3*rP
    soma2 = FO*CpO + FM*CpM + FP*CpP
    
    dTdW = (Uarho*(Ta - T) + soma1)/soma2
    dFMdW = rM    
    dFOdW = rO    
    dFPdW = rP
    
    dydw = [dTdW, dFMdW, dFOdW, dFPdW]
    return dydw

a_CO = []
a_CM = []
a_CP = []

#y0 = T, FM, FO, FP
y0 = [330, 2, 2, 0]
w = np.linspace(0, 100, 100)

sol = odeint(pbr, y0, w)

a_CO = np.asanyarray(a_CO)
a_CM = np.asanyarray(a_CM)
a_CP = np.asanyarray(a_CP)

T = sol[:, 0]
FM = sol[:, 1]
FO = sol[:, 2]
FP = sol[:, 3]


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(w, FM, label='FM')
plt.plot(w, FO, label='FO')
plt.plot(w, FP, label='FP')
plt.legend()
plt.xlabel('W')
plt.ylabel('Fluxo')
plt.grid()
plt.show()

plt.figure(dpi=100, figsize=(9, 4))
plt.plot(w, T, label='T')
plt.legend()
plt.xlabel('W')
plt.ylabel('Temperatura (K)')
plt.grid()
plt.show()