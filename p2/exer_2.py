#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 15:03:08 2020

@author: levandoski
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def BM(y, w):
    X, T = y
    
    Tr = T0 = 300
    Ea = 4000
    R = 1.987
    kr = 0.01
    fa0 = fi0 = 10
    fb0 = 0
    kc = 4242
    ca0 = 0.001
    cpa = 25
    cpb = 25
    cpi = 75
    tetaa = 1
    tetab = fb0/fa0
    tetai = fi0/fa0
    somacps = tetaa*cpa +tetab*cpb + tetai*cpi
    deltacp = cpb - cpa
    
    deltaH = -10000
    # T = 10*X*10000/1000 + T0
    
    ca = ca0*(1 - X)*(T0/T)
    k = kr*np.exp(Ea/R*(1/Tr - 1/T))
    # k = kr
    ra = -(k*kc*ca)/(k + kc)
    dxdw = -ra/fa0 
        
    dtdw = ra*deltaH/(fa0*(somacps + deltacp*X))
    dydw = [dxdw, dtdw]
    
    return dydw


#y = X, T
y0 = [0, 300]
m = 545000
# m = 1020e3
w = np.linspace(0, m, 100)

sol = odeint(BM, y0, w)

x, T = sol[:, 0], sol[:, 1]


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(w, x, 'b')
plt.xlabel('Massa de catalisador')
plt.ylabel('Conversao (X)')
plt.title('Reator X vs n')
plt.grid()
# plt.savefig('convert_n.png')
plt.show()


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(w, T, 'b')
plt.xlabel('Numero de Peneiras')
plt.ylabel('Temperatura (K)')
plt.title('Reator T vs n')
plt.grid()
# plt.savefig('temp_n.png')
plt.show()