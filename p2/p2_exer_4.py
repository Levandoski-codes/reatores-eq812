#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 21:47:16 2020

@author: levandoski
"""
from math import log
por = 0.5
dp = 5e-3
ac = 6*(1-por)/dp
print('ac:', ac)

vis = 2e-5
u=10
densi = 0.8

re = densi*dp*u/vis
print('re:', re)

jd = (0.765/(re**(0.82)) + 0.365/(re**0.386))/por

print('jd:', jd)

dab = 6e-6

sc = vis/(densi*dab)

print('sc:', sc)

kc = (jd*dab*(sc**(1/3))*re)/dp

print('kc:', kc)

x = 0.95
L = u*log(1/(1-x))/(kc*ac)

print('L:', L)