#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:53:07 2020

@author: levandoski
"""
import numpy as np 

def bat(y, t):
    X, T = y
    T = 570
    TR = 298 # K
    
    A = 1.94e15 # min^-1
    Ea = 186000 #J mol^-1
    R = 8.314 # J mol^-1 K^-1
    
    k = A*np.exp(-Ea/(R*T)) # min^-1
    
    deltaHpad = 60000 #J mol^-1
    
    Cpo = 966.15 # J mol^−1 K^−1
    Cpa = 123 # J mol^−1 K^−1
    Cps = 843.15 # J mol^−1 K^−1
    
    
    Mo = 384.6 # g mol^-1
    
#    # dados trocador
#    U = 230 #W m^-2 K^−1    
#    mc = 10 #kg min^-1
#    Cpc = 20 # J g^-1 K^-1    
#    Ta1 = 590 # K  
#    d = 1.2 # m    
#    h = d*((5**(1/2)) - 2)/2 # m    
#    Acin = 2*np.pi*d**2    
#    Acal = (np.pi*h**2)/2    
#    Area = Acin + Acal   
    
    Q = 0 # U*Area*(T - Tf)
    
    rho_o = 900 # kg m^-3
    
    V = 2.5 # m^3
    
    deltaCps = Cps + Cpa - Cpo
    
    
    deltaH = deltaHpad + deltaCps*(T - TR)
    
    mo = V * rho_o # kg
    
    Mo = 0.3846 # kg mol^-1
    Na0 = mo/Mo
    
    
    ra = -k*Na0*(1 - X)/V
    
    dTdt = (Q + (-deltaH)*(-ra*V))/(Na0 * (Cpo + deltaCps*X))
    
    dXdt = (-ra*V)/Na0    
    
    dydt = [dXdt, dTdt]
    return dydt

#y = X, T
y0 = [0, 615] 
t = np.linspace(0, 70, 200) # Será gerada uma solução com 101 amostras uniformemente espaçadas no intervalo 0 <= t <= 10

from scipy.integrate import odeint
sol = odeint(bat, y0, t) # Chamada de odeint para gerar a solução. Para fornecer b e c para pend é usado a argumento args
print(sol)


import matplotlib.pyplot as plt
plt.plot(t, sol[:, 0], 'b', label='X(t)')
#plt.plot(t, sol[:, 1], 'g', label='T(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()


plt.plot(t, sol[:, 1], 'g', label='T(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()



