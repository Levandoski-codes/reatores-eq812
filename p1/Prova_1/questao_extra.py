#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 18:01:22 2020

@author: levandoski
"""
import math



# Dados do problema:
# for = formalina
# F = formaldeido
# nh3 = amonia
# a = agua
# H = HTMA
rho_nh3 = 0.91 #g cm^3
rho_for = 1.1 # g cm ^3
Cp_meio = 4.1868 # J / g K #amonia eh igual

MM_nh3 = 17.0305 # g / mol
MM_F = 30.0260 # g / mol 
MM_H = 140.1863 # g / mol

deltaHpad = -2232.97 # kJ/kg por kg de HTMA

x_F = 0.42 # em massa
V0 = 900 # L
T0 = 25 # C
Ta = 25 # C
U = 480 # W / m^2 K

F_ent = 7.56 # L / min
x_nh3 = 0.25

T = 100 # C

D = 2.54/100 # m

N_F = ((V0 * 1000) * rho_for * x_F) / MM_F # mol de formaldeido
print('N_f:', N_F)

N_nh3 = (4/6) * N_F
print('N_nh3:', N_nh3)

m_nh3 = N_nh3 * MM_nh3
print('m_nh3:', m_nh3)

F_nh3 = F_ent * rho_nh3 * x_nh3 * 1000
print('F_nh3:', F_nh3)

tempo = m_nh3/F_nh3
print('Tempo de reacao:', tempo)

V_nh3 = F_ent*tempo
print('V_nh3:', V_nh3)

deltaH_reac = deltaHpad*(MM_H/(4*MM_nh3))

A = ((F_ent*1000*rho_nh3)/(U*60*(T - Ta))) * \
    (-deltaH_reac*x_nh3 -Cp_meio*(T - T0))#quebra de linha
print('A:', A)

L = A/(math.pi * D)
print('L:', L)

Vs = (L*math.pi*D**2)/4
print('Vs:', Vs*1000)

V = (4/3)*(Vs + V0/1000 + V_nh3/1000)
print('V_t:', V)

d = ((8*V)/(3*math.pi))**(1/3)
print('d:', d)

h = 1.5*d
print('H:', h)
