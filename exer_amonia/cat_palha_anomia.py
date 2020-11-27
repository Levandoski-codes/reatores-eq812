#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 13:56:33 2020

@author: levandoski
"""
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def func_Cp(T, A, B, C, D, E):
    """cal/mol*K"""
    T = T/1000
    return A + B*T + C*T**2 + D*T**3 + E/(T**2)

def func_Cpnh3(T):
    A = 4.779071
    B = 11.89560
    C = -3.67950
    D = 0.45170
    E = 0.045214
    
    return func_Cp(T, A, B, C, D, E)

def func_Cpo2(T):
    A = 7.177904
    B = 2.096791
    C = -0.953187
    D = 0.188411
    E = -0.177246
    
    return func_Cp(T, A, B, C, D, E)

def func_Cpn2(T):
    A = 4.662006
    B = 4.753120
    C =-2.055099	
    D = 0.327386	
    E = 0.126100
    
    return func_Cp(T, A, B, C, D, E)

def func_Cpno(T):
    A = 5.696681
    B = 3.008791
    C =-0.272230	
    D =-0.357901	
    E = 0.051194	
    
    return func_Cp(T, A, B, C, D, E)

def func_Cph2o(T):
    A = 7.192161
    B = 1.633011
    C = 1.623670	
    D = -0.605755	
    E = 0.019632
    
    return func_Cp(T, A, B, C, D, E)


def Cps(T):
    Cpnh3 = func_Cpnh3(T)
    Cpo2 = func_Cpo2(T)
    Cpn2 = func_Cpn2(T)
    Cpno =  func_Cpno(T)
    Cph2o = func_Cph2o(T)
    return Cpnh3, Cpo2, Cpn2, Cpno, Cph2o

def palha(y, v):
    X, T = y
    pi = np.pi
    inch = 39.3701
    
    D = 1 #m peneira
    d = 0.044e-3 #m diametro do fio
    N = 250 # mesh
    R = 8.2057e-5 # constante gases atm
    
    
    T0 = 500 # K
    
    P = 2 # atm
    Dab = 0.247 #cm^2/s
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
    
    rhonh3 = 0.83272 #kg/m^3
    rhoo2 = 1.5593 #kg/m^3
    rhon2 = 1.3644 # kg/m^3
    
    
    Cpnh3, Cpo2, Cpn2, Cpno, Cph2o = Cps(T)# Essa funcao calcula Cps em cada T
    
    Cpnh3 = Cpnh3/MMnh3 #cal/g.K
    Cpo2 = Cpo2/MMo2 #cal/g.K
    Cpn2 = Cpn2/MMn2 #cal/g.K
    Cpno = Cpno/MMno #cal/g.K
    Cph2o = Cph2o/MMh2o #cal/g.K
    
    nimix = ynh3*(minh3/rhonh3) + yo2*(mio2/rhoo2) + yn2*(min2/rhon2) #m^2/s
    
    A = (pi*D**2)/4 # m^2
    L = (1/(N**2) + (d*inch)**2)**(1/2) #in
    e = 1 - (pi*L*d*inch*N**2)/4
    V = A*2*d*(1 - e) #m^3
    
    vazao = ft*R*T0/P #m^3/s
    
    v0nh3 = fnh3*MMnh3/rhonh3 # m^3/s
    v0o2 = fo2*MMo2/rhoo2 # m^3/s
    v0n2 = fn2*MMn2/rhon2 # m^3/s
    
    v0t = v0nh3 + v0n2 + v0o2
    
    DAB = Dab*(P1/P)*(T/T1)**1.75 *1e-4 #m^2/s
    
    U = vazao/A # m/s
    Rey = U*d/(nimix*e)
    
    Jd1_3 = 0.94/(Rey**0.717)
    gama = (1- N*d*inch)**2 #in
    
    Sc = nimix/DAB
    kc3 = (Jd1_3*DAB*(Sc**(1/3))*(Rey**0.5))/d #m/s
    
    Jd1_5 = 0.664/(((Rey/gama)**0.57)*gama)
    kc5 = (Jd1_5*DAB*(Sc**(1/3))*(Rey**0.5))/d #m/s
    
    Pnh3 = ynh3*P #atm
    Pn2 = yn2*P #atm
    Po2 = yo2*P #atm
    
    Cnh3 = Pnh3/(R*T0) #mol/m^3
    Cn2 = Pn2/(R*T0) #mol/m^3
    Co2 = Po2/(R*T0) #mol/m^3
    
    Ct = P/(R*T0) #mol/m^3
    Ct2 = ft/vazao #mol/m^3
    
    epsilon = (6/4 + 4/4 - 5/4 - 1)*fnh3/ft
    
    
    
    ag = pi*L*N**2 #1/in
    
    ag_m = ag*inch #1/m
    
    deltaCp = (6/4*Cph2o + 4/4*Cpno - 5/4*Cpo2 - Cpnh3) #cal/g.K
    SomaCPs = Cpnh3 + (fo2*MMo2)/(fnh3*MMnh3)*Cpo2 #cal/g.K
    
    Fa0 = fnh3*MMnh3 # g/s
    
    deltaHrx = -(54250 - 0.4*(T-298))/MMnh3 #cal/g de NH3
    
    Ca = ((Cnh3*(1-1*X))/(1 + epsilon*X))*(T0/T)
    
    ra = -kc5*ag_m*Ca*MMnh3 # g/m^3s
    
    dXdv = -ra/Fa0 # g/m^3.s / g/s ---> 1/m^3
    
    dTdv = ra*deltaHrx/(Fa0*(SomaCPs + X*deltaCp)) # ---> K/m^3
    
    dydv = [dXdv, dTdv]
    return dydv    


pi = np.pi
inch = 39.3701

D = 1 #m peneira
d = 0.044e-3 #m diametro do fio
N = 250 # mesh
R = 8.2057e-5 # constante gases atm

A = (pi*D**2)/4 # m² Area de um peneira
L = (1/(N**2) + (d*inch)**2)**(1/2) #in
e = 1 - (pi*L*d*inch*N**2)/4 #porosidade
    
Vol = A*2*d*(1 - e) #m³ # Volume de cada peneira

# Vol = d*2*A # m³

numero_peneiras = 5
    
#y = X, T
y0 = [0, 500] 
v = np.linspace(0, numero_peneiras*Vol, 1000)

sol = odeint(palha, y0, v)
X, T = sol[:, 0], sol[:, 1]

v = v/Vol #numero de peneiras


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(v, T, 'b')
plt.xlabel('Numero de Peneiras')
plt.ylabel('Temperatura (K)')
plt.title('Reator T vs n')
plt.grid()
# plt.savefig('temp_n.png')
plt.show()


plt.figure(dpi=100, figsize=(9, 4))
plt.plot(v, X, 'b')
plt.xlabel('Numero de Peneiras')
plt.ylabel('Conversao (X)')
plt.title('Reator X vs n')
plt.grid()
# plt.savefig('convert_n.png')
plt.show()
    
