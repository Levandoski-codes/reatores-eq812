# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def const(T, Ea, A):
    R = 8.314
    k = A*np.exp(-Ea/(R*T))
#    k = A* T**2 + Ea *T**3
    return k


ks = np.array([0.0567, 0.0806, 0.1580, 0.380])
Ts = np.array([10, 15, 25, 40]) + 273.15

pars, _ = curve_fit(const, Ts, ks)

Tas = np.linspace(10, 50, 100) + 273.15

Ea, A = pars[0], pars[1]

plt.scatter(Ts - 273.15, ks)
plt.plot(Tas - 273.15, const(Tas, Ea, A))


def bat(y, t, Ea, A):
    X, T = y
    
    Ca0 = 0.000216 #mol/cm3
    Cps = 3.762 #j/g K
    V = 200e3 # cm3
    
    MM = 102.0886 #g/mol
    Cps = Cps * MM # J/mol K
#    Cps = 4184 # J/kg K
    
    deltaH = -209000 #j/mol
    
    Na0 = Ca0*V
    
    k = const(T, Ea, A)
#    print(k)
#    k = np.mean([0.0567, 0.0806, 0.1580, 0.380])
#    k = np.mean([0.380])
    ra = k*Ca0*(1-X)
    
    dTdt = -deltaH*k*(1-X)/Cps
    dXdt = k*(1-X)
    
    dydt = [dXdt, dTdt]
    return dydt


#y = X, T
y0 = [0, 10 + 273.15] 
t = np.linspace(0, 10, 100) # h

sol = odeint(bat, y0, t, args=(Ea, A))
X, T = sol[:, 0], sol[:, 1]


#k = const(Tas, Ea, A)


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(t, X, 'b')
plt.xlabel('t (min)')
plt.ylabel('Conversão (X)')
plt.title('Reator batelada X vs t')
plt.grid()
plt.show()


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(t, T - 273.15, 'b')
plt.xlabel('t (min)')
plt.ylabel('Temperatura (ºC)')
plt.title('Reator batelada T vs t')
plt.grid()
plt.show()