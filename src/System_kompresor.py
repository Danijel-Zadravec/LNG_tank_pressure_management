# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 11:05:54 2022

@author: faksH
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
from src.tank import Tank_params, Tank_initial, Tank
from src.Pump import PumpStates, LngPump
from src.Evaporator_microchanel import EvapParams, LngEvap, LNGSuperh, GlycolEvap, GlycolSuperh, EvapCommonStates, SuperhCommonStates, EvapStates, Evaporator
from src.Superheater_microchanel import HEparams, LngStream, GlycolStream, StatesCommonHE, StatesHE, HeatExchanger

from src.Pipes import PipeParams, PipeStates, Pipe
from src.Valve import Valve, ValveParams, ValveStates
from src.Prigusenje import ThrottleValveStates, ThrottleValve
from src.rezimi import Rezimi

import src.properties.liq_props as liq

class Engine:
    def __init__(self, demand, pressures, temperatures, times):
        self.demand = demand #fuel demand, kmol/s
        self.fuel_pressures = pressures #demanded fuel pressure
        self.fuel_temperatures = temperatures #demanded fuel pressure
        self.times = times #calculation times, s


class SystemLNG:
    def __init__(self, tank1, tank2, pump, thrvalve, pipe_pump_evaporator, evaporator, valve1, valve2, pipe_evaporator_superheater, superheater, pipe_superheater_engine, engine):
        self.tank1 = tank1
        self.tank2 = tank2
        self.pump = pump
        self.thrvalve = thrvalve
        self.pipe_pump_evap = pipe_pump_evaporator
        self.evaporator = evaporator
        self.valve1 = valve1 #cijevi
        self.valve2 = valve2 #cijevi
        self.pipe_evap_superh = pipe_evaporator_superheater
        self.superheater = superheater
        self.pipe_superh_engine = pipe_superheater_engine
        self.engine = engine
        self.results = pd.DataFrame({'time': self.engine.times})

    def calculate(self):
        times = self.engine.times
        #print(times)
        demand = self.engine.demand
        temperature_glycol = self.engine.fuel_temperatures
        pressures = self.engine.fuel_pressures

        vap_flow1 = 0.0
        vap_flow2 = 0.0
        for i in range(1,len(times)):
            #print(i)
            pressure_evap = pressures[i-1] + 100 #regulacija!


            if i%1000 == 0:
                print("               " + str(i) + "                    ")
            dt = times[i] - times[i-1]
            p_tank1 = self.tank1.states.common.pressure
            p_tank2 = self.tank2.states.common.pressure


            #print(tank1_vap_flow)
            end_success = True
            if self.tank1.states.liquid.vol_ratio > 0.05:
                self.valve1.calculate(p_tank1, pressure_evap, self.tank1.states.liquid.density)
                self.valve2.calculate(p_tank2, pressure_evap, self.tank2.states.liquid.density)
                tank1_vap_flow = self.valve1.states.gas_mol_flow #kmol/s
                tank2_vap_flow = self.valve2.states.gas_mol_flow #kmol/s

                tank1_liq_flow = demand[i] - tank1_vap_flow - tank2_vap_flow
                tank2_liq_flow = 0.0
                #print(tank1_liq_flow)
                tank1_temperature = self.tank1.states.liquid.temperature # temperatura prije odvođenja LNG K
                self.tank1.update_states(-tank1_liq_flow, -tank1_vap_flow, -1, -1, dt)
                self.tank2.update_states(-tank2_liq_flow, -tank2_vap_flow, -1, -1, dt)
                pump_T_in = tank1_temperature
                tank1_pressure =  tank1.states.common.pressure
                pump_pressure = max(pressure_evap+50000.0, tank1_pressure + 50000.0)
                self.pump.calculate(pump_T_in, tank1_pressure, pump_pressure)
                pump_exit_temperature = self.pump.states.T_out
                self.thrvalve.calculate(pump_exit_temperature, pump_pressure, pressure_evap)
                T_thrvalve =self.thrvalve.states.T_out
                self.pipe_pump_evap.calculate(T_thrvalve, -1, pressure_evap, -1, tank1_liq_flow)
                T_evap_in =self.pipe_pump_evap.states.T_out
                self.evaporator.update_states(T_evap_in, T_evap_in+50, pressure_evap, tank1_liq_flow, temperature_glycol[i], 15.0)
                evap_T_out = self.evaporator.states.superh_lng.T_out
                vapor_flow = tank1_liq_flow + (tank1_vap_flow + tank2_vap_flow)
                self.pipe_evap_superh.calculate(evap_T_out, -1, pressure_evap, -1, vapor_flow)
                T_in_superh = self.pipe_evap_superh.states.T_out
                p_in_superh = self.pipe_evap_superh.states.p_out
                self.superheater.update_states(vapor_flow, T_in_superh, 320.0, p_in_superh, 20.0, temperature_glycol[i])
                T_superh_out = self.superheater.states.lng.T_out
                self.pipe_superh_engine.calculate(T_superh_out, -1, p_in_superh, -1, vapor_flow)
            elif self.tank2.states.liquid.vol_ratio > 0.05:
                self.valve1.calculate(p_tank1, pressure_evap, self.tank1.states.liquid.density)
                self.valve2.calculate(p_tank2, pressure_evap, self.tank2.states.liquid.density)
                tank1_vap_flow = self.valve1.states.gas_mol_flow #kmol/s
                tank2_vap_flow = self.valve2.states.gas_mol_flow #kmol/s

                tank2_liq_flow = demand[i] - tank1_vap_flow - tank2_vap_flow
                tank1_liq_flow = 0.0
                tank2_temperature = self.tank2.states.liquid.temperature # temperatura prije odvođenja LNG K
                self.tank1.update_states(-tank1_liq_flow, -tank1_vap_flow, -1, -1, dt)
                self.tank2.update_states(-tank2_liq_flow, -tank2_vap_flow, -1, -1, dt)
                pump_T_in = tank2_temperature
                tank2_pressure =  tank2.states.common.pressure
                pump_pressure = max(pressure_evap+50000.0, tank2_pressure + 50000.0)
                self.pump.calculate(pump_T_in, tank2_pressure, pump_pressure)
                pump_exit_temperature = self.pump.states.T_out
                self.thrvalve.calculate(pump_exit_temperature, pump_pressure, pressure_evap)
                T_thrvalve =self.thrvalve.states.T_out
                self.pipe_pump_evap.calculate(T_thrvalve, -1, pressure_evap, -1, tank2_liq_flow)
                T_evap_in =self.pipe_pump_evap.states.T_out
                self.evaporator.update_states(T_evap_in, T_evap_in+50, pressure_evap, tank2_liq_flow, temperature_glycol[i], 15.0)
                evap_T_out = self.evaporator.states.superh_lng.T_out
                vapor_flow = tank2_liq_flow + (tank1_vap_flow + tank2_vap_flow)
                self.pipe_evap_superh.calculate(evap_T_out, -1, pressure_evap, -1, vapor_flow)
                T_in_superh = self.pipe_evap_superh.states.T_out
                p_in_superh = self.pipe_evap_superh.states.p_out
                self.superheater.update_states(vapor_flow, T_in_superh, 320.0, p_in_superh, 20.0, temperature_glycol[i])
                T_superh_out = self.superheater.states.lng.T_out
                self.pipe_superh_engine.calculate(T_superh_out, -1, p_in_superh, -1, vapor_flow)
            else:
                end_success = False
                break
        if end_success:
            self.results = self.results[self.results.index < (i+1)]
        else:
            self.results = self.results[self.results.index < (i)]

        self.save_results()



        #     if ((tank1.states.liquid.quantity < 0.1) or (tank2.states.liquid.quantity < 0.1)):
        #         self.results =  pd.DataFrame({'time': times[0:i]})
        #         break
        #     else if
        #     if self.tank1.states.common.pressure > 550000.0:
        #         vap_flow1 = -0.005
        #     elif self.tank1.states.common.pressure < 300000.0:
        #         vap_flow1 = 0.0
        #     if self.tank2.states.common.pressure > 550000.0:
        #         vap_flow2 = -0.005
        #     elif self.tank2.states.common.pressure < 300000.0:
        #         vap_flow2 = 0.0
        #     tank1_liq_flow = -demand[i]/2 - vap_flow1  #dodatna kapljevina za održavanje tlaka
        #     tank2_liq_flow = -demand[i]/2 - vap_flow2  #dodatna kapljevina za održavanje tlaka

        #     tank_vap_temperature = -1
        #     tank1_temperature = self.tank1.states.liquid.temperature # temperatura prije odvođenja LNG K
        #     tank2_temperature = self.tank2.states.liquid.temperature # temperatura prije odvođenja LNG K

        #     tank1_pressure = self.tank1.states.common.pressure
        #     tank2_pressure = self.tank2.states.common.pressure

        #     #print('tnk')
        #     self.tank1.update_states(tank1_liq_flow, vap_flow1, -1, tank_vap_temperature, dt)
        #     self.tank2.update_states(tank2_liq_flow, vap_flow2, -1, tank_vap_temperature, dt)
        #     pump_T_in = (tank1_temperature + tank2_temperature) / 2.0
        #     pump_p_in = (tank1_pressure + tank2_pressure) / 2
        #     self.pump.calculate(pump_T_in, pump_p_in, pump_pressure)
        #     pump_exit_temperature = self.pump.states.T_out
        #     #print('pp1')
        #     total_liq_flow =- (tank1_liq_flow + tank2_liq_flow)
        #     self.pipe_pump_evap.calculate(pump_exit_temperature, -1, pump_pressure, -1, total_liq_flow)
        #     evap_in_T = self.pipe_pump_evap.states.T_out
        #     evap_p_in = self.pipe_pump_evap.states.p_out
        #     #print('evap')
        #     self.evaporator.update_states(total_liq_flow, evap_in_T, evap_p_in, 323.15)
        #     evap_T_out = self.evaporator.states.lng.T_sat
        #     #print('pp2')
        #     vapor_flow = total_liq_flow - (vap_flow1 + vap_flow2)
        #     self.pipe_evap_superh.calculate(evap_T_out, -1, evap_p_in, -1, vapor_flow)
        #     T_in_superh = self.pipe_evap_superh.states.T_out
        #     p_in_superh = self.pipe_evap_superh.states.p_out
        #     T_superh_out = 0.1 + temperature[i]
        #     #print('superh')
        #     self.superheater.update_states(vapor_flow, T_in_superh, T_superh_out, 323.15, p_in_superh)
        #     #vap_return = vapor_flow - demand[i]
        #     #print('pp3')
        #     self.pipe_superh_engine.calculate(T_superh_out, -1, p_in_superh, -1, vapor_flow)

        # self.save_results()

    def save_results(self):
        '''
        https://stackoverflow.com/questions/21443963/pandas-multilevel-column-names

        '''
        col = 1
        cols = [(' ', 'time')]
        for k in self.tank1.save.liquid:
            self.results[col] = self.tank1.save.liquid[k]
            cols.append(('tank1 liq', k))
            col = col + 1
        for k in self.tank1.save.vapor:
            self.results[col] = self.tank1.save.vapor[k]
            cols.append(('tank1 vap', k))
            col = col + 1
        for k in self.tank1.save.common:
            self.results[col] = self.tank1.save.common[k]
            cols.append(('tank1 com', k))
            col = col + 1
        for k in self.tank2.save.liquid:
            self.results[col] = self.tank2.save.liquid[k]
            cols.append(('tank2 liq', k))
            col = col + 1
        for k in self.tank2.save.vapor:
            self.results[col] = self.tank2.save.vapor[k]
            cols.append(('tank2 vap', k))
            col = col + 1
        for k in self.tank2.save.common:
            self.results[col] = self.tank2.save.common[k]
            cols.append(('tank2 com', k))
            col = col + 1
        for k in self.pump.save.states:
            self.results[col] = self.pump.save.states[k]
            cols.append(('pump', k))
            col = col + 1
        for k in self.thrvalve.save.states:
            #print(self.thrvalve.save.states[k])
            self.results[col] = pd.Series(self.thrvalve.save.states[k])
            cols.append(('throttle valve', k))
            col = col + 1
        for k in self.pipe_pump_evap.save.states:
            self.results[col] = pd.Series(self.pipe_pump_evap.save.states[k])
            cols.append(('p. pmp evap', k))
            col = col + 1

        for k in self.evaporator.save.evap_lng:
            self.results[col] =  pd.Series(self.evaporator.save.evap_lng[k])
            cols.append(('evap lng evap', k))
            col = col + 1
        for k in self.evaporator.save.evap_glyc:
            self.results[col] =  pd.Series(self.evaporator.save.evap_glyc[k])
            cols.append(('evap glycol evap', k))
            col = col + 1
        for k in self.evaporator.save.evap_com:
            self.results[col] =  pd.Series(self.evaporator.save.evap_com[k])
            cols.append(('evap comm evap', k))
            col = col + 1
        for k in self.evaporator.save.superh_lng:
            self.results[col] =  pd.Series(self.evaporator.save.superh_lng[k])
            cols.append(('evap lng superh', k))
            col = col + 1
        for k in self.evaporator.save.superh_glyc:
            self.results[col] =  pd.Series(self.evaporator.save.superh_glyc[k])
            cols.append(('evap glycol superh', k))
            col = col + 1
        for k in self.evaporator.save.superh_com:
            self.results[col] = self.evaporator.save.superh_com[k]
            cols.append(('evap comm superh', k))
            col = col + 1
        for k in self.valve1.save.states:
            self.results[col] = self.valve1.save.states[k]
            cols.append(('BOG valve 1', k))
            col = col + 1
        for k in self.valve2.save.states:
            self.results[col] = self.valve2.save.states[k]
            cols.append(('BOG valve 2', k))
            col = col + 1
        for k in self.pipe_evap_superh.save.states:
            self.results[col] = self.pipe_evap_superh.save.states[k]
            cols.append(('p. evp suph', k))
            col = col + 1
        for k in self.superheater.save.lng:
            self.results[col] = self.superheater.save.lng[k]
            cols.append(('super lng', k))
            col = col + 1
        for k in self.superheater.save.glycol:
            self.results[col] = self.superheater.save.glycol[k]
            cols.append(('super glycol', k))
            col = col + 1
        for k in self.superheater.save.common:
            self.results[col] = self.superheater.save.common[k]
            cols.append(('super comm', k))
            col = col + 1
        for k in self.pipe_superh_engine.save.states:
            self.results[col] = self.pipe_superh_engine.save.states[k]
            cols.append(('p. suph eng', k))
            col = col + 1
        #print(cols)
        self.results.columns=pd.MultiIndex.from_tuples(cols)
        self.results = self.results.copy()





def create_engine_inputs(times, demands, pressures, temperatures):
    times_inp = []
    demands_inp = []
    pressures_inp = []
    temperatures_inp = []
    for i in range(1,len(times)):
        times_tmp = np.arange(times[i-1], times[i], 10.0)
        len_lst = len(times_tmp)
        times_inp = times_inp + times_tmp.tolist()
        demands_inp = demands_inp + [demands[i-1]]*len_lst
        pressures_inp = pressures_inp + [pressures[i-1]]*len_lst
        temperatures_inp = temperatures_inp + [temperatures[i-1]]*len_lst
    return (times_inp, demands_inp, pressures_inp, temperatures_inp)





if __name__ == "__main__":
#     #primjer koristenja:
    rezim = Rezimi()

    # SPREMNIK
    int_length = 14.5 #m
    int_diameter = 5.9 #m
    thickness = 0.05 #m
    k_liq = 1/(0.08/0.015+1/2000+1/20) #W/mK #excel
    k_vap = 1/(0.08/0.015+1/200+1/20) #W/mK #excel
    pressure = 100000. #Pa
    temperature_l = 110 #K
    temperature_v = 113 #K
    liq_vol_ratio = 0.83 #
    T_amb = 293.15 #K
    tank_params1 = Tank_params(int_length, int_diameter, thickness, k_liq, k_vap, T_amb)
    tank_initial1 = Tank_initial(pressure, temperature_l, temperature_v, liq_vol_ratio)
    tank1 = Tank(tank_params1,tank_initial1)

    tank_params2 = Tank_params(int_length, int_diameter, thickness, k_liq, k_vap, T_amb)
    tank_initial2 = Tank_initial(pressure, temperature_l, temperature_v, liq_vol_ratio)
    tank2 = Tank(tank_params2,tank_initial2)


    #PUMPA
    eta = 0.8
    pmp_sts = PumpStates(temperature_l, temperature_v, pressure, 600000)
    lng_pump = LngPump(eta, pmp_sts)


    #THROTTLE VALVES
    thrv_st = ThrottleValveStates()
    thr_val = ThrottleValve(thrv_st)


    #CIJEV PUMPA ISPARIVAC
    pipe_pmp_d = 0.1
    pipe_pmp_length = 1
    pipe_pmp_roughness = 0.00015 #m
    pipe_pmp_medium ='lng_liq'
    pipe_pmp_k = 2 #W/m2K
    T_env = T_amb
    pipe_pmp_params =PipeParams(pipe_pmp_d, pipe_pmp_length, pipe_pmp_roughness, pipe_pmp_medium, pipe_pmp_k, T_env)
    pipe_pmp_T_in =290
    pipe_pmp_T_out =-1
    pipe_pmp_p_in = 600000.0
    pipe_pmp_p_out = -1
    pipe_pmp_flow = 0.08 #kmol/s
    pipe_pmp_states = PipeStates(pipe_pmp_T_in, pipe_pmp_T_out, pipe_pmp_p_in, pipe_pmp_p_out, pipe_pmp_flow)
    pipe_pmp = Pipe(pipe_pmp_params, pipe_pmp_states)

    #ISPARIVAC
    #ev_areas_lst = [9.0, 8.0, 7.0, 6.0, 5.5 ,5., 4.5, 4., 3.5, 3., 2.5, 2.25, 2.0, 1.75, 1.5, 1.3, 1.25, 1.125, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2] #m2
    #k_ev = 1200. #W/m2K
    #evaporator_params = EVparams(ev_areas_lst, k_ev)
    #evaporator_lng =EV_lng(temperature_l, pressure, pipe_pmp_flow)

    #evaporator_glycol =EV_glycol(T_glycol_in, T_glycol_out, glycol_flow)
    #st_common = StatesCommonEV()
    #evaporator_states = EV_States(evaporator_lng, evaporator_glycol, st_common)
    #lng_evaporator = Evaporator(evaporator_params, evaporator_states)

    m_ev = 22 #broj prolaza
    h_ev = 0.01
    t_ev = 0.00015
    n_ev = 500
    s_ev = 0.0012
    l_ev = 172
    W_ev = 0.43
    L_ev = 0.405 #duljina, m
    T_glycol_in = 318.0
    T_glycol_out = 318.0
    glycol_flow = 0.0

    evaporator_params = EvapParams(m_ev, h_ev, t_ev, n_ev, s_ev, l_ev, W_ev, L_ev)
    evap_lng = LngEvap(temperature_l, pressure, pipe_pmp_flow)
    evap_glyc = GlycolEvap(T_glycol_in-15, T_glycol_in-30, glycol_flow)
    evap_com = EvapCommonStates()
    superh_lng = LNGSuperh(T_glycol_out, pressure, pipe_pmp_flow)
    superh_glyc = GlycolSuperh(T_glycol_in, T_glycol_in-15, glycol_flow)
    superh_com = SuperhCommonStates()
    evaporator_states = EvapStates(evap_com, evap_lng, evap_glyc, superh_com, superh_lng, superh_glyc)
    evaporator = Evaporator(evaporator_params, evaporator_states)




    #VENTILI

    Kv1 = 1
    Kv2 = 1
    valve_par1 = ValveParams(Kv1)
    valve_par2 = ValveParams(Kv2)
    valve_st1 = ValveStates()
    valve_st2 = ValveStates()
    valve1 = Valve(valve_par1,valve_st1)
    valve2 = Valve(valve_par2,valve_st2)

    #PIPE EVAP SUPERHEATER
    evap_pipe_d = 0.4
    evap_pipe_length = 2
    evap_pipe_roughness = 0.00015 #m
    evap_pipe_medium ='lng_vap'
    evap_pipe_k = 2 #W/m2K
    evap_pipe_params =PipeParams(evap_pipe_d, evap_pipe_length, evap_pipe_roughness, evap_pipe_medium, evap_pipe_k, T_env)


    evap_pipe_T_in =290.0
    evap_pipe_T_out = -1
    evap_pipe_p_in = 600000.0
    evap_pipe_p_out = -1
    evap_pipe_flow = 0.03 #kmol/s
    evap_pipe_states = PipeStates(evap_pipe_T_in, evap_pipe_T_out, evap_pipe_p_in, evap_pipe_p_out, evap_pipe_flow)
    evap_pipe = Pipe(evap_pipe_params, evap_pipe_states)

    #SUPERHEATER
    m_sup = 32
    h_sup = 0.01
    t_sup = 0.00015
    n_sup = 500
    s_sup = 0.0012
    l_sup = 172
    W_sup = 0.75
    L_sup = 0.75 #duljina, m

    he_params = HEparams(m_sup, h_sup, t_sup, n_sup, s_sup, l_sup, W_sup, L_sup)
    p = 600000. #Pa
    T_lng_in = 293.0 #K
    T_lng_out0 = 293.0 #K
    lng_flow = 0.0 #kmol/s
    stream = 'lng'
    he_lng =LngStream(T_lng_in, T_lng_out0, lng_flow, p,)

    T_glycol_in = 318.0
    T_glycol_out = 318.0
    glycol_flow = 0.0
    stream = 'glycol'

    he_glycol =GlycolStream(T_glycol_in, glycol_flow)
    st_common = StatesCommonHE()
    he_states = StatesHE(he_lng, he_glycol, st_common)

    he_superheater = HeatExchanger(he_params, he_states)



    eng_pipe_d = 0.4
    eng_pipe_length = 5
    eng_pipe_roughness = 0.00015 #m
    eng_pipe_medium ='lng_vap'
    eng_pipe_k = 2 #W/m2K
    T_env = T_amb
    eng_pipe_params =PipeParams(eng_pipe_d, eng_pipe_length, eng_pipe_roughness, eng_pipe_medium, eng_pipe_k, T_env)


    eng_pipe_T_in =-1
    eng_pipe_T_out =290.0
    eng_pipe_p_in = -1
    eng_pipe_p_out = 500000.0
    eng_pipe_flow = 0.2 #kmol/s
    eng_pipe_states = PipeStates(eng_pipe_T_in, eng_pipe_T_out, eng_pipe_p_in, eng_pipe_p_out, eng_pipe_flow)
    engine_pipe = Pipe(eng_pipe_params, eng_pipe_states)


    times = [0.01, 1*24*60*60-60*30, 1*24*60*60+30*60, 2*24*60*60-60*30, 2*24*60*60+60*30, 3*24*60*60-60*30, 3*24*60*60+60*30, 4*24*60*60-60*30, 4*24*60*60+60*30, 5*24*60*60-60*30]
    #times = [0.01, 200, 500, 2*60*60*2, 2*60*60*3]
    demands = [rezim.port_winter, rezim.manouvering_winter, rezim.service_winter, rezim.manouvering_winter, rezim.port_winter, rezim.manouvering_winter, rezim.service_winter, rezim.manouvering_winter, rezim.port_winter,] #2. 0.03

    pres = 600000.0 #6 bar -prema dokumentu  P1381 Inquiry No.02-707-000 LNG System-Final.pdf
    pressures = [pres]*len(demands)
    temp = liq.saturation_temperature(pres) + 5.0
    temperatures_glyc = [273.15+40.0]*len(demands) #  -prema dokumentu  P1381 Inquiry No.02-707-000 LNG System-Final.pdf
    (times_inp, demands_inp, pressures_inp, temperatures_inp) = create_engine_inputs(times, demands, pressures, temperatures_glyc)

    time_end = 24*60*6*3-1
    endd = 10
    engine = Engine(demands_inp[0:], pressures_inp[0:], temperatures_inp[0:], times_inp[0:])
    system = SystemLNG(tank1, tank2, lng_pump, thr_val, pipe_pmp, evaporator, valve1, valve2, evap_pipe, he_superheater, engine_pipe, engine)
    system.calculate()
    rezultati = system.results

    plt.figure()
    plt.plot(tank1.save.common['pressure'])
    #plt.figure()
    plt.plot(tank2.save.common['pressure'])
    plt.plot(tank2.save.liquid['temperature'])
    plt.plot(tank2.save.common['T_sat'])


    plt.figure()
    plt.plot(tank2.save.common['evaporation'])
    plt.figure()
    plt.plot(tank1.save.vapor['flow_quantity'])
    plt.plot(tank1.save.liquid['flow_quantity'])

    plt.figure()
    plt.plot(tank1.save.liquid['quantity'])
    plt.plot(tank2.save.liquid['quantity'])
    plt.plot(tank2.save.liquid['flow_quantity'])

    plt.plot(tank1.save.liquid['temperature'])
    plt.plot(tank1.save.common['T_sat'])
    thr_val.states.p_in
    thr_val.states.p_out
    lng_pump.states.T_in
    lng_pump.states.p_in
    lng_pump.states.p_out

    from src.properties.liq_props_lowlev import MetLiq
    met_liq = MetLiq()

    met_liq.set_T(lng_pump.states.T_in)
    met_liq.get_p()
    lng_pump.states.p_in
    lng_pump.states.p_out
    lng_pump.states.T_in

    res_all = system.results
    res_all_head = res_all.head()
    times_used = res_all[' ']['time'].to_numpy()
    max_time = times_used[-1]
    #%%
    plt.figure(figsize=(8,6))
    tms = times_used/60/60
    pressures1 = np.array(system.tank1.save.common['pressure']) / 100000.0
    pressures2 = np.array(system.tank2.save.common['pressure']) / 100000.0

    plt.plot(tms, pressures1, label='spremnik 1')
    plt.plot(tms, pressures2, label='spremnik 2')
    plt.legend()
    plt.ylabel('tlak u spremniku, bar')
    plt.xlabel('vrijeme, h')
    plt.grid()
    plt.savefig('1dan_tlak.jpg')


    #%%
    plt.figure(figsize=(8,6))
    tms = times_used/60/60
    T_liq1 = np.array(system.tank1.save.liquid['temperature'])
    T_liq2 = np.array(system.tank2.save.liquid['temperature'])
    T_vap1 = np.array(system.tank1.save.vapor['temperature'])
    T_vap2 = np.array(system.tank2.save.vapor['temperature'])

    plt.plot(tms, T_liq1, label='spremnik 1 kapljevina')
    plt.plot(tms, T_liq2, label='spremnik 2 kapljevina')
    plt.plot(tms, T_vap1, label='spremnik 1 para')
    plt.plot(tms, T_vap2, label='spremnik 2 para')

    plt.legend()
    plt.ylabel('temperatura, K')
    plt.xlabel('vrijeme, h')
    plt.grid()
    plt.savefig('1dan_spremnik_temperatura.jpg')


    #%%


    plt.figure(figsize=(8,6))
    tms = times_used/60/60
    vol_rat1 = np.array(system.tank1.save.liquid['vol_ratio']) *100
    vol_rat2 = np.array(system.tank2.save.liquid['vol_ratio']) *100

    plt.plot(tms, vol_rat1, label='spremnik 1')
    plt.plot(tms, vol_rat2, label='spremnik 2')
    plt.legend()
    plt.ylabel('udio volumena kapljevine, %')
    plt.xlabel('vrijeme, h')
    plt.grid()
    plt.savefig('1dan_vol_kapljevine.jpg')


    #%%


    plt.figure(figsize=(8,6))
    tms = times_used/60/60
    flow_liq1 = np.array(system.tank1.save.liquid['flow'])* 16.04*3600
    flow_liq2= np.array(system.tank2.save.liquid['flow']) * 16.04*3600
    flow_vap1 = np.array(system.tank1.save.vapor['flow']) * 16.04*3600
    flow_vap2 = np.array(system.tank2.save.vapor['flow']) * 16.04*3600

    plt.plot(tms, flow_liq1, label='kapljevina spremnik 1')
    plt.plot(tms, flow_liq2, label='kapljevina spremnik 2')
    plt.plot(tms, flow_vap1, label='para spremnik 1')
    plt.plot(tms, flow_vap2, label='para spremnik 2')


    plt.legend()
    plt.ylabel('protok kg/h')
    plt.xlabel('vrijeme, h')
    plt.grid()
    plt.savefig('1dan_protok.jpg')




    #%%


    plt.figure(figsize=(8,6))
    tms = times_used/60/60
    T_tcs = np.array(system.pipe_superh_engine.save.states['T_out']) - 273.15

    plt.plot(tms, T_tcs, label='TCS')
    plt.legend()
    plt.ylabel('temperatura, °C')
    plt.xlabel('vrijeme, h')
    plt.grid()
    plt.savefig('1dan_temperatura_.jpg')



#%%



    #%%
    plt.figure(figsize=(8,6))
    times_plt = times_used

    T_l = np.array(system.tank1.save.liquid['temperature']) -273.15
    T_l = T_l[1:]
    T_v = np.array(system.tank1.save.vapor['temperature']) -273.15
    T_v = T_v[1:]
    T_evap =  np.array(system.evaporator.save.lng['T_sat']) -273.15
    T_evap = T_evap[1:]
    T_superh = np.array(system.superheater.save.lng['T_out']) -273.15
    T_superh = T_superh[1:]
    plt.plot(times_plt, T_l, label='spremnik kapljevina')
    plt.plot(times_plt, T_v, label='spremnik para')
    #plt.plot(times, T_evap, label='isparivač')
    #plt.plot(times, T_superh, label='pregrijač')

    plt.legend()
    plt.ylabel('temperatura, °C')
    plt.xlabel('vrijeme, h')
    plt.grid()

    plt.savefig('trial_temperature.jpg')


    #%%
    plt.figure(figsize=(8,6))
    T_l = np.array(system.tank1.save.liquid['quantity'])
    T_l = T_l[1:]
    T_v = np.array(system.tank1.save.vapor['quantity'])
    T_v = T_v[1:]

    plt.plot(times_plt, T_l, label='kapljevina')
    plt.plot(times_plt, T_v, label='para')
    #plt.plot(times, T_evap, label='isparivač')
    #plt.plot(times, T_superh, label='pregrijač')

    plt.legend()
    plt.ylabel('količina u spremniku, kmol')
    plt.xlabel('vrijeme, h')
    plt.grid()

    plt.savefig('trial_quantity.jpg')



    #%%


    plt.figure(figsize=(8,6))
    T_l = np.array(system.evaporator.save.glycol['T_out'])-273.15
    T_l = T_l[1:]
    T_v = np.array(system.superheater.save.glycol['T_out'])-273.15
    T_v = T_v[1:]

    plt.plot(times_plt, T_l, label='isparivač')
    plt.plot(times_plt, T_v, label='pregrijač')
    #plt.plot(times, T_evap, label='isparivač')
    #plt.plot(times, T_superh, label='pregrijač')

    plt.legend()
    plt.ylabel('temperatura glikola na izlazu, °C')
    plt.xlabel('vrijeme, h')
    plt.grid()

    plt.savefig('trial_T_glycol.jpg')


    #%%
    plt.figure(figsize=(8,6))
    T_l = np.array(system.evaporator.save.glycol['mas_flow'])
    T_l = T_l[1:]
    T_v = np.array(system.superheater.save.glycol['mas_flow'])
    T_v = T_v[1:]

    plt.plot(times_plt, T_l, label='isparivač')
    plt.plot(times_plt, T_v, label='pregrijač')
    #plt.plot(times, T_evap, label='isparivač')
    #plt.plot(times, T_superh, label='pregrijač')

    plt.legend()
    plt.ylabel('protok glikola na izlazu, kg/s')
    plt.xlabel('vrijeme, h')
    plt.grid()

    plt.savefig('trial_flow_glycol.jpg')


    #%%
    plt.figure(figsize=(8,6))
    T_l = np.array(system.engine.demand)
    T_l = T_l[1:]

    plt.plot(times_plt, T_l)
    plt.legend()
    plt.ylabel('protok LNG, kmol/s')
    plt.xlabel('vrijeme, h')
    plt.grid()

    #%%


    plt.figure()
    plt.plot(tank1.save.common['pressure'])
    plt.figure()
    plt.plot(tank2.save.common['pressure'])

    plt.figure()
    plt.plot(tank2.save.common['evaporation'])
    plt.figure()
    plt.plot(tank1.save.vapor['flow_quantity'])
    plt.figure()
    plt.plot(tank1.save.liquid['quantity'])


aaa=valve1.save.states
aaaa=evaporator.save.superh_com


