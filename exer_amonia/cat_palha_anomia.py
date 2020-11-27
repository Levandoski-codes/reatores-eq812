#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 13:56:33 2020

@author: levandoski
"""
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def palha(y, v):
    X, T = y
    pi = np.pi
    inch = 39.3701
    
    D = 1 #m peneira
    d = 0.044e-3 #m diâmetro do fio
    N = 250 # mesh
    R = 8.2057e-5 # constante gases atm
    
    
    T0 = 500 # K
    # T = T0
    
    P = 2 # atm
    Dab = 0.247 #cm²/s
    P1 = 1 #atm
    T1 = 30 + 273.15 #K
    
    fnh3 = 10 #mol/s
    fo2 = (5/4)*fnh3*1.25 #mol/s com 25% de excesso de o2
    fn2 = fo2*(79/21) #mol/s
    ft = fnh3 + fo2 + fn2
    
    ynh3 = fnh3/ft
    yo2 = fo2/ft
    yn2 = fn2/ft
    
    
    minh3 = 1.7864e-5 #Pa.s
    mio2 = 3.0420e-5 #Pa.s
    min2 = 2.6037e-5 # Pa.s
    
    MMno = 30.0061 #g/mol
    MMh2o = 18.0153 #g/mol
    MMnh3 = 17.0305 #g/mol
    MMo2 = 31.9988 #g/mol
    MMn2 = 28.0134 #g/mol
    
    
    Cpnh3 = 0.59502 #cal/g.K
    Cpo2 = 0.23246 #cal/g.K
    Cpn2 = 0.18148 #cal/g.K
    Cpno_mol = 7.29 #cal/mol.K
    Cpno = Cpno_mol / MMno #cal/g.K
    Cph2o = 0.48021	#cal/g.K
    
    Cpnh3 = 0.59502 * MMnh3 #cal/g.K
    Cpo2 = 0.23246 * MMo2 #cal/g.K
    Cpn2 = 0.18148 * MMn2 #cal/g.K
    Cpno_mol = 7.29 #cal/mol.K
    Cpno = Cpno_mol / MMno  * MMno #cal/g.K
    Cph2o = 0.48021 * MMh2o	#cal/g.K
    
    rhonh3 = 0.83272 #kg/m³    
    rhoo2 = 1.5593 #kg/m³
    rhon2 = 1.3644 # kg/m³
    
    nimix = ynh3*(minh3/rhonh3) + yo2*(mio2/rhoo2) + yn2*(min2/rhon2) #m²/s
    
    A = (pi*D**2)/4 # m²
    L = (1/(N**2) + (d*inch)**2)**(1/2) #in
    e = 1 - (pi*L*d*inch*N**2)/4
    V = A*2*d*(1 - e) #m³
    
    vazao = ft*R*T0/P #m³/s
    
    v0nh3 = fnh3*MMnh3/rhonh3 # m³/s
    v0o2 = fo2*MMo2/rhoo2 #m³/s
    v0n2 = fn2*MMn2/rhon2 #m³/s
    
    v0t = v0nh3 + v0n2 + v0o2
    
    DAB = Dab*(P1/P)*(T/T1)**1.75 *1e-4 #m²/s
    # print(DAB)
    
    U = vazao/A # m/s
    Rey = U*d/(nimix*e)
    
    Jd1_3 = 0.94/(Rey**0.717)
    gama = (1- N*d*inch)**2 #in
    gamaSI = gama/inch # m
    
    Sc = nimix/DAB
    kc = (Jd1_3*DAB*(Sc**(1/3))*(Rey**0.5))/d #m/s
    
    Jd1_5 = 0.664/(((Rey/gama)**0.57)*gama)
    kc2 = (Jd1_5*DAB*(Sc**(1/3))*(Rey**0.5))/d #m/s
    
    Pnh3 = ynh3*P #atm
    Pn2 = yn2*P #atm
    Po2 = yo2*P #atm
    
    Cnh3 = Pnh3/(R*T0) #mol/m³
    Cn2 = Pn2/(R*T0) #mol/m³
    Co2 = Po2/(R*T0) #mol/m³
    
    Ct = P/(R*T0) #mol/m³
    Ct2 = ft/vazao #mol/m³
    
    epsilon = (6/4 + 4/4 - 5/4 - 1)*fnh3/ft
    
    
    
    ag = pi*L*N**2 #1/in
    
    ag_m = ag*inch #1/m
    
    deltaCp = (6/4*Cph2o + 4/4*Cpno - 5/4*Cpo2 - Cpnh3) #cal/g.K
    SomaCPs = Cpnh3 + (fo2*MMo2)/(fnh3*MMnh3)*Cpo2 #cal/g.K (fo2*(MMo2*1e3))/(fnh3*(MMnh3*1e3))   (fo2/fnh3)
    
    
    Fa0 = fnh3*MMnh3 # o Fa0 foi convertido para g por s
    
    deltaHrx = -(54250 - 0.4*(T-298))/MMnh3 #cal/mol de nh3 converti pra cal por g
    
    # X = 0.1
    Ca = ((Cnh3*(1-1*X))/(1 + epsilon*X))*(T0/T) # inverti o T0 e tirei o 4 do x
    
    ra = -kc*ag_m*Ca*MMnh3
    
    dXdv = -ra/Fa0
    
    dTdv = ra*deltaHrx/(Fa0*(SomaCPs + X*deltaCp))
    
    dydv = [dXdv, dTdv]
    return dydv    


pi = np.pi
inch = 39.3701

D = 1 #m peneira
d = 0.044e-3 #m diâmetro do fio
N = 250 # mesh
R = 8.2057e-5 # constante gases atm

A = (pi*D**2)/4 # m²
L = (1/(N**2) + (d*inch)**2)**(1/2) #in
e = 1 - (pi*L*d*inch*N**2)/4
    
Vol = A*2*d*(1 - e) #m³

Vol = d*2*A # m³
print(Vol)



    
    #y = X, T
y0 = [0, 500] 
t = np.linspace(0, 5*Vol, 100)

sol = odeint(palha, y0, t)
X, T = sol[:, 0], sol[:, 1]

t = t/Vol


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(t, T, 'b')
plt.xlabel('Número de Peneiras')
plt.ylabel('Temperatura (K)')
plt.title('Reator T vs V')
plt.grid()
plt.show()


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(t, X, 'b')
plt.xlabel('Número de Peneiras')
plt.ylabel('Conversão (X)')
plt.title('Reator X vs V')
plt.grid()
plt.show()
    
