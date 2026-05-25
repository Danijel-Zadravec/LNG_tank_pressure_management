# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import src.properties.vap_props as vap
import src.properties.liq_props as liq
import src.properties.glycol_props as glyc
from heat_exchanger import HEparams, HE_stream, StatesHE, HeatExchanger
area = 10 #m2
k = 1000 #W/m2K



def pi1_calc(pi2, pi3):
    tmp = np.exp(-(1-pi3)*pi2)
    pi1=(1-tmp)/(1-pi3*tmp)
    return pi1

def pi1_cond_calc(pi2):
    pi1=1-np.exp(-pi2)
    return pi1

pi2 = np.linspace(0.001, 8., 200)
pi1_0 = pi1_calc(pi2, 0.0)
pi1_1 = pi1_calc(pi2, 0.9999999)





area = 10 #m2
k = 1000 #W/m2K
he_params = HEparams(area, k)
p = 400000 #Pa
T_lng_in = 293.0 #K
T_lng_out = 293.0 #K
lng_flow = 0.0 #kmol/s
stream = 'lng'
he_lng =HE_stream(T_lng_in, T_lng_out, lng_flow,  p, stream)

T_glycol_in = 318.0
T_glycol_out = 318.0
glycol_flow = 0.0
stream = 'glycol'

he_glycol =HE_stream(T_glycol_in, T_glycol_out, glycol_flow, p, stream)
he_states = StatesHE(he_lng, he_glycol)

he = HeatExchanger(he_params, he_states)
lng_T_in = liq.saturation_temperature(p)+2



flows = np.linspace(0.04,0.16,13)
pi_1_calc =np.zeros_like(flows)
pi_2_calc =np.zeros_like(flows)
pi_3_calc =np.zeros_like(flows)
T_gl_out_lst =np.zeros_like(flows)
flow_gl_lst =np.zeros_like(flows)

for i in range(len(flows)):
    he.update_states(flows[i], lng_T_in, T_lng_out, T_glycol_in, p)
    pi_1_calc[i] = he.states.pi1
    pi_2_calc[i]= he.states.pi2
    pi_3_calc[i] = he.states.pi3
    T_gl_out_lst[i] = he.states.glycol.T_out
    flow_gl_lst[i] =he.states.glycol.mol_flow



plt.figure
plt.plot(pi2, pi1_0, label='pi3=0')
plt.plot(pi2, pi1_1, label='pi3=1')
plt.scatter(pi_2_calc, pi_1_calc, label='rezultati')
plt.legend()
he.update_states(0.2, lng_T_in, T_lng_out, T_glycol_in, p)
pi1_calc(1.4882757567288893, 5.9608609865491405e-06)
