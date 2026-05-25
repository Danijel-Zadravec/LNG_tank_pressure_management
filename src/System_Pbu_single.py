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
from src.Evaporator_microchanel import EvapParams, LngEvap, LNGSuperh, GlycolEvap, GlycolSuperh, EvapCommonStates, SuperhCommonStates, EvapStates, Evaporator
from src.Superheater_microchanel import HEparams, LngStream, GlycolStream, StatesCommonHE, StatesHE, HeatExchanger
from src.PBU_microchanel import PBUParams, LngPBU, GlycolPBU, PBUStates, PBU
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

class SystemLNG_PBU:
    def __init__(self, tank, pbu, thrvalve, pipe_pump_evaporator, evaporator, valve, pipe_evaporator_superheater, superheater, pipe_superheater_engine, evap_flow, super_flow, pbu_flow):
        self.pbu = pbu
        self.tank = tank
        self.thrvalve = thrvalve
        self.pipe_pump_evap = pipe_pump_evaporator
        self.evaporator = evaporator
        self.valve = valve #cijevi
        self.pipe_evap_superh = pipe_evaporator_superheater
        self.superheater = superheater
        self.pipe_superh_engine = pipe_superheater_engine
        self.engine = None
        self.results = None
        self.time = 0
        self.protok_evap = evap_flow
        self.protok_super = super_flow
        self.protok_pbu = pbu_flow
        self.not_working_times = []
        self.holding_time = None
        self.BOG_excess = 0.0


    def calculate(self, engine):

        time_holding = 0.0

        self.engine = engine
        times = self.engine.times
        #print(times)
        demand = self.engine.demand
        temperature_glycol = self.engine.fuel_temperatures
        pressures = self.engine.fuel_pressures
        pbu_counter = 0
        pbu_state = False
        pbu2_state = False

        h24s = 24*60*60 #sekunde u 24 dana
        h48s = 48*60*60 #sekunde u 48 dana

        bog_state = False
        histereza = 0.0


        not_w = False
        time_notw = 0.0 # ne radi zbog preniskog tlaka

        for i in range(1,len(times)):
            #print(i)
            pressure_evap = pressures[i-1] + 100 #regulacija!
            p_min_noflow = pressure_evap + 50000.0-100.0 # nema oduzimanja kapljevine - ne moze se ppostici protok
            p_PBU_ON = pressure_evap + 50000.0 + 10000.0
            p_pbu_max = pressure_evap + 150000.0 - 100.0 #do ovog tlaka radi PBU
            p_BOG_on = pressure_evap + 160000.0 - 100.0 #odvođenje BOG iz spremnika - previsok tlak
            p_BOG_off = pressure_evap + 100000.0 - 100.0
            if i%1000 == 0:
                print("               " + str(i) + "                    ")
            dt = times[i] - times[i-1]
            p_tank = self.tank.states.common.pressure


            #print(tank_vap_flow)
            end_success = True

            self.time = self.time + dt

            tank_day = (((self.time % h48s) - h24s) < 0)


            tank_vol_ratio = self.tank.states.liquid.vol_ratio
            #while (p_tank < p_min_noflow) and (p_tank2 < p_min_noflow): #prenizak tlak u oba


            while not((p_tank > p_min_noflow)): #prenizak tlak u oba

                if not not_w:
                    time_notw = 0.0
                    not_w = True

                T_tankV = self.tank.states.vapor.temperature
                pbu_counter = pbu_counter+1
                print(pbu_counter)
                # pustanje BOG - nema jer je tlak jako nizak
                self.valve.calculate(p_tank, T_tankV,  p_tank, self.tank.states.liquid.density) #stavljen isti tlak tako da nema protoka

                #PBU spremnik 1
                T_liq_t = self.tank.states.liquid.temperature
                self.pbu.update_states(T_liq_t, p_tank, temperature_glycol[i], self.protok_pbu)
                if tank_vol_ratio > 0.04:
                    pbu_flow = self.pbu.states.evap_lng.mol_flow
                else:
                    pbu_flow = 0.0
                T_vap = self.pbu.states.evap_lng.T_sat
                self.tank.update_states(-pbu_flow, pbu_flow, -1, T_vap, dt)
                p_tank = self.tank.states.common.pressure

                tank_T = T_liq_t
                tank_p = p_tank
                self.thrvalve.no_flow(tank_T, tank_p)
                self.pipe_pump_evap.no_flow(tank_T, tank_p)
                self.evaporator.no_flow(tank_T, tank_p, temperature_glycol[i])
                self.pipe_evap_superh.no_flow(tank_T, tank_p)
                self.superheater.no_flow(tank_T, tank_p, temperature_glycol[i])
                self.pipe_superh_engine.no_flow(tank_T, tank_p)

                #print(self.pipe_pump_evap.states.flow)

                time_notw = time_notw + dt
                if (tank_vol_ratio < 0.04): #premalo lng u oba spr
                    break

            if not_w:
                self.not_working_times = self.not_working_times + [time_notw]
                not_w = False

            tank_vol_ratio = self.tank.states.liquid.vol_ratio



            if (tank_vol_ratio > (0.04+histereza)) and (p_tank >= p_min_noflow):
                T_tankV = self.tank.states.vapor.temperature
                histereza = 0.0
                if p_tank<p_PBU_ON:
                    pbu_state = True
                elif p_tank>p_pbu_max:
                    pbu_state = False


                if p_tank > p_BOG_on:
                    bog_state = True
                elif p_tank < p_BOG_off:
                    bog_state = False



                if bog_state:
                    self.valve.calculate(p_tank, T_tankV, pressure_evap, self.tank.states.liquid.density)


                else:
                    self.valve.calculate(p_tank, T_tankV, p_tank, self.tank.states.liquid.density)

                tank_vap_flow = self.valve.states.gas_mol_flow #kmol/s

                # tank_vap_flow = min(demand[i]-0.00002, tank_vap_flow)
                # tank_vap_flow = max(tank_vap_flow, 0.0)

                # self.valve.states.gas_mol_flow = tank_vap_flow


                tank_liq_flow = demand[i] - tank_vap_flow


                if tank_liq_flow < 0.0:
                    self.BOG_excess = self.BOG_excess - tank_liq_flow*dt
                    tank_liq_flow = 0.0


                #print(tank_liq_flow)
                #assert tank_liq_flow >= 0.000005

                tank_temperature = self.tank.states.liquid.temperature # temperatura prije odvođenja LNG K
                if pbu_state == False:
                    self.pbu.ne_radi()
                    self.tank.update_states(-tank_liq_flow, -tank_vap_flow, -1, -1, dt)

                else:
                    self.pbu.update_states(T_liq_t, p_tank, temperature_glycol[i], self.protok_pbu)
                    pbu_flow = self.pbu.states.evap_lng.mol_flow
                    tank_vap_flow_tank = pbu_flow
                    tank_liq_flow_tank = tank_liq_flow + pbu_flow
                    T_vap = self.pbu.states.evap_lng.T_sat
                    self.tank.update_states(-tank_liq_flow_tank, tank_vap_flow_tank, -1, T_vap, dt)



                self.thrvalve.calculate(tank_temperature, p_tank, pressure_evap)
                T_thrvalve =self.thrvalve.states.T_out
                self.pipe_pump_evap.calculate(T_thrvalve, -1, pressure_evap, -1, tank_liq_flow)
                T_evap_in =self.pipe_pump_evap.states.T_out
                self.evaporator.update_states(T_evap_in, T_evap_in+50, pressure_evap, tank_liq_flow, temperature_glycol[i], self.protok_evap)
                evap_T_out = self.evaporator.states.superh_lng.T_out
                valve_T = self.valve.states.T_out
                vapor_flow = demand[i]
                T_mixer_out = evap_T_out

                self.pipe_evap_superh.calculate(T_mixer_out, -1, pressure_evap, -1, vapor_flow)
                T_in_superh = self.pipe_evap_superh.states.T_out
                p_in_superh = self.pipe_evap_superh.states.p_out
                self.superheater.update_states(vapor_flow, T_in_superh, 320.0, p_in_superh, self.protok_super, temperature_glycol[i])
                T_superh_out = self.superheater.states.lng.T_out
                self.pipe_superh_engine.calculate(T_superh_out, -1, p_in_superh, -1, vapor_flow)
                #print(self.pipe_pump_evap.states.flow)


            else:
                end_success = False
                break
        it_ugas = 0
        p = self.tank.states.common.pressure

        if end_success:
            time_sv = np.arange(i+pbu_counter+1+it_ugas)*10
        else:
            time_sv = np.arange(i+pbu_counter+it_ugas)*10
        self.results  = pd.DataFrame({'time': time_sv})
        self.save_results()



    def save_results(self):
        '''
        https://stackoverflow.com/questions/21443963/pandas-multilevel-column-names

        '''
        col = 1
        cols = [(' ', 'time')]
        for k in self.tank.save.liquid:
            self.results[col] = self.tank.save.liquid[k]
            cols.append(('tank liq', k))
            col = col + 1
        for k in self.tank.save.vapor:
            self.results[col] = self.tank.save.vapor[k]
            cols.append(('tank vap', k))
            col = col + 1
        for k in self.tank.save.common:
            self.results[col] = self.tank.save.common[k]
            cols.append(('tank com', k))
            col = col + 1
        for k in self.pbu.save.evap_lng:
            self.results[col] = self.pbu.save.evap_lng[k]
            cols.append(('pbu lng', k))
            col = col + 1
        for k in self.pbu.save.evap_glyc:
            self.results[col] = self.pbu.save.evap_glyc[k]
            cols.append(('pbu glycol', k))
            col = col + 1
        for k in self.pbu.save.evap_com:
            self.results[col] = self.pbu.save.evap_com[k]
            cols.append(('pbu common', k))
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
        for k in self.valve.save.states:
            self.results[col] = self.valve.save.states[k]
            cols.append(('BOG valve 1', k))
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
