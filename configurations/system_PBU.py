# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:31:20 2022

@author: faksH
"""

from src.tank import Tank_params, Tank_initial, Tank
from src.Evaporator_microchanel import EvapParams, LngEvap, LNGSuperh, GlycolEvap, GlycolSuperh, EvapCommonStates, SuperhCommonStates, EvapStates, Evaporator
from src.Superheater_microchanel import HEparams, LngStream, GlycolStream, StatesCommonHE, StatesHE, HeatExchanger
from src.PBU_microchanel import PBUParams, LngPBU, GlycolPBU, PBUStates, PBU

from src.Pipes import PipeParams, PipeStates, Pipe
from src.Valve import Valve, ValveParams, ValveStates
from src.Prigusenje import ThrottleValveStates, ThrottleValve
from src.System_Pbu_single import SystemLNG_PBU
from src.Mixer_vapor import MixerVstates, MixerV



# SPREMNIK
int_diameter = 7.2 #m
int_length = 18.5 #m
lmbd = 0.0343
delt = 0.1
thickness = 0.05 #m
k_liq = 1/(delt/lmbd+1/2000+1/25) #W/mK #excel
k_vap = 1/(delt/lmbd+1/200+1/25) #W/mK #excek
pressure = 100000. #Pa
temperature_l = 110.0 #K
temperature_v = 113.0 #K
liq_vol_ratio = 0.83 #
T_amb = 293.15 #K








#slosh_times = [20000, 40000, 60000, 80000 , 100000, 120000, 140000, 160000, 180000]


#PBU
T_glycol_in = 313.15
m = 4 #broj prolaza
h=0.01
t=0.00015
n=500
s =0.0012
l=172
W = 0.295
L=0.28 #duljina, m
qm_glyc_pbu = 5.0



#CIJEV PUMPA ISPARIVAC
pipe_pmp_d = 0.1
pipe_pmp_length = 5
pipe_pmp_roughness = 0.00015 #m
pipe_pmp_medium ='lng_liq'
pipe_pmp_k = 2 #W/m2K
T_env = T_amb
pipe_pmp_T_in =290
pipe_pmp_T_out =-1
pipe_pmp_p_in = 600000.0
pipe_pmp_p_out = -1
pipe_pmp_flow = 0.08 #kmol/s


#ISPARIVAC
m_ev = 23 #broj prolaza
h_ev = 0.01
t_ev = 0.00015
n_ev = 500
s_ev = 0.0012
l_ev = 172
W_ev = 0.458
L_ev = 0.405 #duljina, m
T_glycol_out = 318.0
glycol_flow = 0.0




#VENTILI

Kv = 0.2


#PIPE EVAP SUPERHEATER
evap_pipe_d = 0.4
evap_pipe_length = 5
evap_pipe_roughness = 0.00015 #m
evap_pipe_medium ='lng_vap'
evap_pipe_k = 2 #W/m2K


evap_pipe_T_in =290.0
evap_pipe_T_out = -1
evap_pipe_p_in = 600000.0
evap_pipe_p_out = -1
evap_pipe_flow = 0.03 #kmol/s



#SUPERHEATER
m_sup = 32
h_sup = 0.01
t_sup = 0.00015
n_sup = 500
s_sup = 0.0012
l_sup = 172
W_sup = 0.75
L_sup = 0.75 #duljina, m


p = 600000. #Pa
T_lng_in = 293.0 #K
T_lng_out0 = 293.0 #K
lng_flow = 0.0 #kmol/s

T_glycol_in = 313.0
T_glycol_out = 313.0
glycol_flow = 0.0




eng_pipe_d = 0.4
eng_pipe_length = 5
eng_pipe_roughness = 0.00015 #m
eng_pipe_medium ='lng_vap'
eng_pipe_k = 2 #W/m2K
T_env = T_amb


eng_pipe_T_in =-1
eng_pipe_T_out =290.0
eng_pipe_p_in = -1
eng_pipe_p_out = 500000.0
eng_pipe_flow = 0.2 #kmol/s





def create_system_PBU(slosh_times, evap_flow, super_flow, pbu_flow):
    tank_params = Tank_params(int_length, int_diameter, thickness, k_liq, k_vap, T_amb, slosh_times)
    tank_initial = Tank_initial(pressure, temperature_l, temperature_v, liq_vol_ratio)
    tank = Tank(tank_params,tank_initial)



    paramsPBU = PBUParams(m, h, t, n, s, l, W, L)
    evap_lng1 = LngPBU(temperature_l, pressure)
    evap_glyc1 = GlycolPBU(T_glycol_in, T_glycol_in-5, qm_glyc_pbu)
    states1 = PBUStates(evap_lng1, evap_glyc1)
    pbu = PBU(paramsPBU, states1)



    #THROTTLE VALVES
    thrv_st = ThrottleValveStates()
    thr_val = ThrottleValve(thrv_st)
    pipe_pmp_params =PipeParams(pipe_pmp_d, pipe_pmp_length, pipe_pmp_roughness, pipe_pmp_medium, pipe_pmp_k, T_env)
    pipe_pmp_states = PipeStates(pipe_pmp_T_in, pipe_pmp_T_out, pipe_pmp_p_in, pipe_pmp_p_out, pipe_pmp_flow)
    pipe_pmp = Pipe(pipe_pmp_params, pipe_pmp_states)


    evaporator_params = EvapParams(m_ev, h_ev, t_ev, n_ev, s_ev, l_ev, W_ev, L_ev)
    evap_lng = LngEvap(temperature_l, pressure, pipe_pmp_flow)
    evap_glyc = GlycolEvap(T_glycol_in-15, T_glycol_in-30, glycol_flow)
    evap_com = EvapCommonStates()
    superh_lng = LNGSuperh(T_glycol_out, pressure, pipe_pmp_flow)
    superh_glyc = GlycolSuperh(T_glycol_in, T_glycol_in-15, glycol_flow)
    superh_com = SuperhCommonStates()
    evaporator_states = EvapStates(evap_com, evap_lng, evap_glyc, superh_com, superh_lng, superh_glyc)
    evaporator = Evaporator(evaporator_params, evaporator_states)

    valve_par = ValveParams(Kv)
    valve_st = ValveStates()
    valve = Valve(valve_par,valve_st)

    evap_pipe_params =PipeParams(evap_pipe_d, evap_pipe_length, evap_pipe_roughness, evap_pipe_medium, evap_pipe_k, T_env)
    evap_pipe_states = PipeStates(evap_pipe_T_in, evap_pipe_T_out, evap_pipe_p_in, evap_pipe_p_out, evap_pipe_flow)
    evap_pipe = Pipe(evap_pipe_params, evap_pipe_states)


    he_params = HEparams(m_sup, h_sup, t_sup, n_sup, s_sup, l_sup, W_sup, L_sup)
    he_lng =LngStream(T_lng_in, T_lng_out0, lng_flow, p,)

    he_glycol =GlycolStream(T_glycol_in, glycol_flow)
    st_common = StatesCommonHE()
    he_states = StatesHE(he_lng, he_glycol, st_common)

    he_superheater = HeatExchanger(he_params, he_states)

    eng_pipe_params =PipeParams(eng_pipe_d, eng_pipe_length, eng_pipe_roughness, eng_pipe_medium, eng_pipe_k, T_env)

    eng_pipe_states = PipeStates(eng_pipe_T_in, eng_pipe_T_out, eng_pipe_p_in, eng_pipe_p_out, eng_pipe_flow)
    engine_pipe = Pipe(eng_pipe_params, eng_pipe_states)

    systemPBU = SystemLNG_PBU(tank, pbu, thr_val, pipe_pmp, evaporator, valve, evap_pipe, he_superheater, engine_pipe, evap_flow, super_flow, pbu_flow)
    return systemPBU




