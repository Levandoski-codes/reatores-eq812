#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np 
from numpy import exp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#isso s[o funcionar[a no colab ou em jupyter notebook

def plots():
    plt.figure(dpi=100, figsize=(9, 4))
    plt.plot(z, T, 'r', label='$T$')
    plt.plot(z, Ta, 'b', label='$T_a$')
    plt.xlabel('z (m)')
    plt.ylabel('Temperatura (K)')
    plt.title('Reator T vs z')
    plt.xlim(0, L)
    plt.ylim(300, 1000)
    plt.legend()
    plt.grid()
    plt.savefig('temp.png', dpi=200)
    plt.show()

    Fn =  Fnf - (1/2)*(Fa - Faf)
    Fh =  Fhf - (3/2)*(Fa - Faf)

    plt.figure(dpi=100, figsize=(9, 4))
    plt.plot(z, Fa, 'r', label='$F_{NH_3}$')
    plt.plot(z, Fn, 'g', label='$F_{N_2}$')
    plt.plot(z, Fh, 'b', label='$F_{H_2}$')
    plt.xlabel('z (m)')
    plt.ylabel('Fs')
    plt.title('Reator Fis vs z')
    plt.xlim(0, L)
    plt.ylim(0, 180)
    plt.legend()
    plt.grid()
    plt.savefig('Fis.png', dpi=200)
    plt.show()

    Ft = Fa + Fn + Fh
    ya = Fa/Ft
    yn = Fn/Ft
    yh = Fh/Ft

    plt.figure(dpi=100, figsize=(9, 4))
    plt.plot(z, ya, 'r', label='$y_{NH_3}$')
    plt.plot(z, yn, 'g', label='$y_{N_2}$')
    plt.plot(z, yh, 'b', label='$y_{H_2}$')
    plt.xlabel('z (m)')
    plt.ylabel('Ys')
    plt.title('Reator yis vs z')
    plt.xlim(0, L)
    plt.ylim(0, 1)
    plt.legend()
    plt.grid()
    plt.savefig('Yis.png', dpi=200)
    plt.show()

def edo(y, z):
    global X
    FA, T, Ta = y
    K298 = exp(-deltaG/Rg/T)
    K = K298*exp((deltaHr/Rg)*(1/Tr - 1/T))

    #Velocidade especifica
    k1 = k10*exp(-ER/T) #Pa/s

    #Determinacao da vazao molar inicial
    FT0 = v0*(P)/Rg/T0 #mol/s

    #Vazão molar inicial de cada composto
    FAf = FT0*x0nh3 #mol/s
    FH2f = FT0*x0h2 #mol/s
    FN2f = FT0*x0n2 #mol/s

    #Vazão molar de nitrogênio e hidrogênio com o decorrer da reação
    FN2 = FN2f - 0.5*(FA - FAf) #mol/s
    FH2 = FH2f - 1.5*(FA - FAf) #mol/s

    #Vazão molar total com o decorrer da reação
    FT = FA + FN2 + FH2 #mol/s


    #Pressão
    yA = FA/FT
    yN2 = FN2/FT
    yH2 = FH2/FT
    Pa = P*yA #Pa
    Ph2 = P*yH2 #Pa
    Pn2 = P*yN2 #Pa

    #Conversão de nitrogênio
    X = 1 - FN2/FN2f

    #Equação da taxa
    r = (k1/(Rg*T))*((K**2)*(((Pn2)*Ph2**1.5)/Pa)-Pa/(Ph2**1.5))

    #Balanço de massa
    dFAdz= 2*Ac*r # 2xAc

    #Balanço de energia 
    dTdz = -beta*r + gama*(Ta - T) #K/m

    #Balanço de energia do fluido refrigerante
    dTadz= gama*(Ta - T) #K/m

    dydz = [dFAdz, dTdz, dTadz]
    return dydz

global X
#Dados fornecidos pelo enunciado
P = 300*101325 #Pa
L = 12 #m
Ac = 1 #m²
v0 = 0.16
x0nh3 = 0.015
x0n2 = 0.985*(1/4)
x0h2 = 0.985*(3/4)
gama = 0.5 #1/m
beta = -2.342 #m²Ks/mol
k10 = (7.794*10**11)*101325 #Pa/s
ER = (2*10**4) #K
deltaHr = (-1.2*10**4)*4.184 #J/mol
deltaG = (4.25*10**3)*4.184 #J/mol

#Constante universal dos gases
Rg = 8.3145 #J/mol/K

#Constante de equilíbrio
Tr = 298 #K

i = 0
Taf = 323 #323

chute = 400
while True and i<200:
    i+=1
    if i > 1:
        Tfinal = Ta[-1]
        dif = (Tfinal - Taf)
        erro = dif/Taf
        print(f'i: {i} x: {X} Tentrada: {Tentrada} Tfinal:{Tfinal} delta: {dif} erro: {erro}')
        if abs(erro) < 1e-4:
            break
        if Tfinal > Taf:
            Tentrada = Tentrada + dif/30
        else:
            Tentrada = Tentrada - dif/30
    else:
        Tentrada = chute
    z = np.linspace(0, L, 100)
    
    T0 = Tentrada #K
    FT0 = v0*(P)/Rg/T0 #mol/s
    Fa0 = FT0*x0nh3
    # Na, T, Ta
    y0 = Fa0, Tentrada, Tentrada
    sol = odeint(edo, y0, z)
    Fa, T, Ta = sol[:, 0], sol[:, 1], sol[:, 2]

plots()
