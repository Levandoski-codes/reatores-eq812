#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:53:07 2020

@author: levandoski
"""
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def bat(y, t, troca):
    X, T = y
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
    
#    # dados trocador         <----------- ########################
    if troca:
        U = 230 #W m^-2 K^−1    
        mc = 10 #kg min^-1
        Cpc = 20 # J g^-1 K^-1    
        Ta1 = 5900 # K  
        d = 1.2 # m    
        h = d*((5**(1/2)) - 2)/2 # m    
        Acin = 2*np.pi*d**2    
        Acal = np.pi*(d/2)*h    
        Area = Acin + Acal   
        
        Q = mc*Cpc*(T - Ta1)*(1 - np.exp(-U*Area/(mc*Cpc)))
    else:
        Q = 0
    
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
t = np.linspace(0, 31, 200)

troca = True

sol = odeint(bat, y0, t, args=(troca, )) # Chamada de odeint para gerar a solução.
#print(sol[-20:])
X_q, T_q = sol[:, 0], sol[:, 1]

troca = False
sol = odeint(bat, y0, t, args=(troca, )) # Chamada de odeint para gerar a solução.
#print(sol[-20:])
X, T = sol[:, 0], sol[:, 1]

plt.figure(dpi=100)
plt.plot(t, X, 'b', label='X(t)')
plt.plot(t, X_q, 'r', label='Com Q, X(t)')
plt.legend(loc='best')
plt.xlabel('t (min)')
plt.grid()
plt.show()


plt.figure(dpi=100)
plt.plot(t, T, 'g', label='T(t)')
plt.plot(t, T_q, 'r', label='Com Q, T(t)')
plt.legend(loc='best')
plt.xlabel('t (min)')
plt.grid()
plt.show()


#X_q, T_q = sol[:, 0], sol[:, 1]


