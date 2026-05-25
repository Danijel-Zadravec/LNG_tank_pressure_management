# -*- coding: utf-8 -*-
"""
J.E. Hesselgreaves, Richard Law, David Reay - Compact Heat Exchangers. Selection, Design and Operation
str 265
"""
import numpy as np

b = 2/1000 #m
lambd = 15/1000 #m
phi = 60 #degrees
rho = 1000
mass_flow = 0.2
Bp=0.3 #plate width
nu = 1e-6

X =2*np.pi * b / lambd
ench =1/6 *(1+(1+X**2)**(1/2)+4*(1+(X**2)/2)**1/2)
dh =2*b/ench

dh=2*b
u = mass_flow/rho/b/Bp

Re = u*dh/nu


f=0.63*Re**(-0.23)*(1+0.9*(phi-30)/30)

Pr=0.7
Nu=0.205*Pr**(1/3) *1*(f*Re*Re*np.sin(2*np.pi/3))**0.374
