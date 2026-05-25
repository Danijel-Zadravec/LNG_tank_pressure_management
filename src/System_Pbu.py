# -*- coding: utf-8 -*-
"""
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
    def __init__(self, tank1, tank2, pbu1, pbu2, thrvalve, pipe_pump_evaporator, evaporator, valve1, valve2, mixer, pipe_evaporator_superheater, superheater, pipe_superheater_engine,  evap_flow, super_flow, pbu_flow):
        self.pbu1 = pbu1
        self.pbu2 = pbu2
        self.tank1 = tank1
        self.tank2 = tank2
        self.thrvalve = thrvalve
        self.pipe_pump_evap = pipe_pump_evaporator
        self.evaporator = evaporator
        self.valve1 = valve1 #cijevi
        self.valve2 = valve2 #cijevi
        self.mixer = mixer
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


    def calculate(self, engine):

        time_holding = 0.0

        self.engine = engine
        times = self.engine.times
        #print(times)
        demand = self.engine.demand
        temperature_glycol = self.engine.fuel_temperatures
        pressures = self.engine.fuel_pressures
        pbu_counter = 0
        vap_flow1 = 0.0
        vap_flow2 = 0.0
        pbu1_state = False
        pbu2_state = False

        h24s = 24*60*60 #sekunde u 24 dana
        h48s = 48*60*60 #sekunde u 48 dana

        bog1_state = False
        bog2_state = False
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
            p_tank1 = self.tank1.states.common.pressure
            p_tank2 = self.tank2.states.common.pressure


            #print(tank1_vap_flow)
            end_success = True

            self.time = self.time + dt

            tank1_day = (((self.time % h48s) - h24s) < 0)


            tank1_vol_ratio = self.tank1.states.liquid.vol_ratio
            tank2_vol_ratio = self.tank2.states.liquid.vol_ratio
            #while (p_tank1 < p_min_noflow) and (p_tank2 < p_min_noflow): #prenizak tlak u oba


            while not((p_tank1 > p_min_noflow) and (tank1_vol_ratio>0.05)) and not((p_tank2 > p_min_noflow) and (tank2_vol_ratio>0.05)): #prenizak tlak u oba

                if not not_w:
                    time_notw = 0.0
                    not_w = True

                T_tankV_1 = self.tank1.states.vapor.temperature
                T_tankV_2 = self.tank2.states.vapor.temperature
                pbu_counter = pbu_counter+1
                print(pbu_counter)
                # pustanje BOG - nema jer je tlak jako nizak
                self.valve1.calculate(p_tank1, T_tankV_1,  p_tank1, self.tank1.states.liquid.density) #stavljen isti tlak tako da nema protoka
                self.valve2.calculate(p_tank2, T_tankV_2, p_tank2, self.tank2.states.liquid.density)

                #PBU spremnik 1
                T_liq_t1 = self.tank1.states.liquid.temperature
                self.pbu1.update_states(T_liq_t1, p_tank1, temperature_glycol[i], self.protok_pbu)
                if tank1_vol_ratio > 0.05:
                    pbu1_flow = self.pbu1.states.evap_lng.mol_flow
                else:
                    pbu1_flow = 0.0
                T_vap1 = self.pbu1.states.evap_lng.T_sat
                self.tank1.update_states(-pbu1_flow, pbu1_flow, -1, T_vap1, dt)
                p_tank1 = self.tank1.states.common.pressure

                # PBU spremnik 2
                T_liq_t2 = self.tank2.states.liquid.temperature
                self.pbu2.update_states(T_liq_t2, p_tank2, temperature_glycol[i], self.protok_pbu)
                if tank2_vol_ratio > 0.05:
                    pbu2_flow = self.pbu2.states.evap_lng.mol_flow
                else:
                    pbu2_flow = 0.0
                T_vap2 = self.pbu2.states.evap_lng.T_sat
                self.tank2.update_states(-pbu2_flow, pbu2_flow, -1, T_vap2, dt)
                p_tank2 = self.tank2.states.common.pressure
                if p_tank1 > p_tank2:
                    tank_T = T_liq_t1
                    tank_p = p_tank1
                else:
                    tank_T = T_liq_t2
                    tank_p = p_tank2
                self.thrvalve.no_flow(tank_T, tank_p)
                self.pipe_pump_evap.no_flow(tank_T, tank_p)
                self.evaporator.no_flow(tank_T, tank_p, temperature_glycol[i])
                self.pipe_evap_superh.no_flow(tank_T, tank_p)
                self.superheater.no_flow(tank_T, tank_p, temperature_glycol[i])
                self.pipe_superh_engine.no_flow(tank_T, tank_p)
                self.mixer.ne_radi()

                #print(self.pipe_pump_evap.states.flow)

                time_notw = time_notw + dt
                if (tank1_vol_ratio < 0.05) and (tank2_vol_ratio < 0.05): #premalo lng u oba spr
                    break

            if not_w:
                self.not_working_times = self.not_working_times + [time_notw]
                not_w = False

            tank1_vol_ratio = self.tank1.states.liquid.vol_ratio
            tank2_vol_ratio = self.tank2.states.liquid.vol_ratio



            if (tank1_vol_ratio > (0.05+histereza)) and (p_tank1 >= p_min_noflow) and (tank1_day or self.tank1.states.liquid.vol_ratio < 0.05 or p_tank2 <= p_min_noflow):
                T_tankV_1 = self.tank1.states.vapor.temperature
                T_tankV_2 = self.tank2.states.vapor.temperature
                histereza = 0.0
                if p_tank1<p_PBU_ON:
                    pbu1_state = True
                elif p_tank1>p_pbu_max:
                    pbu1_state = False
                if p_tank2<p_PBU_ON:
                    pbu2_state = True
                elif p_tank2>p_pbu_max:
                    pbu2_state = False

                if p_tank1 > p_BOG_on:
                    bog1_state = True
                elif p_tank1 < p_BOG_off:
                    bog1_state = False
                if  p_tank2 > p_BOG_on:
                    bog2_state = True
                elif p_tank2 < p_BOG_off:
                    bog2_state = False


                if bog1_state:
                    self.valve1.calculate(p_tank1, T_tankV_1, pressure_evap, self.tank1.states.liquid.density)
                else:
                    self.valve1.calculate(p_tank1, T_tankV_1, p_tank1, self.tank1.states.liquid.density)

                if bog2_state:
                    self.valve2.calculate(p_tank2, T_tankV_2, pressure_evap, self.tank2.states.liquid.density)
                else:
                    self.valve2.calculate(p_tank2, T_tankV_2, p_tank2, self.tank2.states.liquid.density)

                tank1_vap_flow = self.valve1.states.gas_mol_flow #kmol/s
                tank2_vap_flow = self.valve2.states.gas_mol_flow #kmol/s

                tank1_vap_flow = min(demand[i]-0.00002, tank1_vap_flow)
                tank1_vap_flow = max(tank1_vap_flow, 0.0)
                tank2_vap_flow = min(demand[i]-tank1_vap_flow-0.00001, tank2_vap_flow)
                tank2_vap_flow = max(tank2_vap_flow, 0.0)

                self.valve1.states.gas_mol_flow = tank1_vap_flow
                self.valve2.states.gas_mol_flow = tank2_vap_flow


                tank1_liq_flow = demand[i] - tank1_vap_flow - tank2_vap_flow
                tank2_liq_flow = 0.0
                #print(tank1_liq_flow)
                assert tank1_liq_flow >= 0.000005

                tank1_temperature = self.tank1.states.liquid.temperature # temperatura prije odvođenja LNG K
                if pbu1_state == False:
                    self.pbu1.ne_radi()
                    self.tank1.update_states(-tank1_liq_flow, -tank1_vap_flow, -1, -1, dt)

                else:
                    self.pbu1.update_states(T_liq_t1, p_tank1, temperature_glycol[i], self.protok_pbu)
                    pbu1_flow = self.pbu1.states.evap_lng.mol_flow
                    tank1_vap_flow_tank = pbu1_flow
                    tank1_liq_flow_tank = tank1_liq_flow + pbu1_flow
                    T_vap1 = self.pbu1.states.evap_lng.T_sat
                    self.tank1.update_states(-tank1_liq_flow_tank, tank1_vap_flow_tank, -1, T_vap1, dt)

                if (pbu2_state == False) or (tank2_vol_ratio<0.05):
                    self.pbu2.ne_radi()
                    self.tank2.update_states(-tank2_liq_flow, -tank2_vap_flow, -1, -1, dt)
                else:
                    T_liq_t2 = self.tank2.states.liquid.temperature
                    self.pbu2.update_states(T_liq_t2, p_tank2, temperature_glycol[i], self.protok_pbu)
                    pbu2_flow = self.pbu2.states.evap_lng.mol_flow
                    T_vap2 = self.pbu2.states.evap_lng.T_sat
                    self.tank2.update_states(-pbu2_flow, pbu2_flow, -1, T_vap2, dt)


                self.thrvalve.calculate(tank1_temperature, p_tank1, pressure_evap)
                T_thrvalve =self.thrvalve.states.T_out
                self.pipe_pump_evap.calculate(T_thrvalve, -1, pressure_evap, -1, tank1_liq_flow)
                T_evap_in =self.pipe_pump_evap.states.T_out
                self.evaporator.update_states(T_evap_in, T_evap_in+50, pressure_evap, tank1_liq_flow, temperature_glycol[i], self.protok_evap)
                evap_T_out = self.evaporator.states.superh_lng.T_out
                valve1_T = self.valve1.states.T_out
                valve2_T = self.valve2.states.T_out
                self.mixer.calculate(evap_T_out, valve1_T, valve2_T, tank1_liq_flow, tank1_vap_flow, tank2_vap_flow, pressure_evap)
                vapor_flow = self.mixer.states.qn_tot
                T_mixer_out = self.mixer.states.T_out

                self.pipe_evap_superh.calculate(T_mixer_out, -1, pressure_evap, -1, vapor_flow)
                T_in_superh = self.pipe_evap_superh.states.T_out
                p_in_superh = self.pipe_evap_superh.states.p_out
                self.superheater.update_states(vapor_flow, T_in_superh, 320.0, p_in_superh, self.protok_super, temperature_glycol[i])
                T_superh_out = self.superheater.states.lng.T_out
                self.pipe_superh_engine.calculate(T_superh_out, -1, p_in_superh, -1, vapor_flow)
                #print(self.pipe_pump_evap.states.flow)

            elif (tank2_vol_ratio > 0.05) and (p_tank2 >= p_min_noflow):
                T_tankV_1 = self.tank1.states.vapor.temperature
                T_tankV_2 = self.tank2.states.vapor.temperature
                histereza = 0.02
                if p_tank1<p_PBU_ON:
                    pbu1_state = True
                elif p_tank1>p_pbu_max:
                    pbu1_state = False
                if p_tank2<p_PBU_ON:
                    pbu2_state = True
                elif p_tank2>p_pbu_max:
                    pbu2_state = False


                if p_tank1 > p_BOG_on:
                    bog1_state = True
                elif p_tank1 < p_BOG_off:
                    bog1_state = False
                if  p_tank2 > p_BOG_on:
                    bog2_state = True
                elif p_tank2 < p_BOG_off:
                    bog2_state = False


                if bog1_state:
                    self.valve1.calculate(p_tank1, T_tankV_1, pressure_evap, self.tank1.states.liquid.density)
                else:
                    self.valve1.calculate(p_tank1, T_tankV_1, p_tank1, self.tank1.states.liquid.density)

                if bog2_state:
                    self.valve2.calculate(p_tank2, T_tankV_2, pressure_evap, self.tank2.states.liquid.density)
                else:
                    self.valve2.calculate(p_tank2, T_tankV_2, p_tank2, self.tank2.states.liquid.density)


                tank1_vap_flow = self.valve1.states.gas_mol_flow #kmol/s
                tank2_vap_flow = self.valve2.states.gas_mol_flow #kmol/s






                tank1_vap_flow = min(demand[i]-0.00002, tank1_vap_flow)
                tank1_vap_flow = max(tank1_vap_flow, 0.0)

                tank2_vap_flow = min(demand[i]-tank1_vap_flow-0.00001, tank2_vap_flow)
                tank2_vap_flow = max(tank2_vap_flow, 0.0)

                self.valve1.states.gas_mol_flow = tank1_vap_flow
                self.valve2.states.gas_mol_flow = tank2_vap_flow

                tank2_liq_flow = demand[i] - tank1_vap_flow - tank2_vap_flow
                tank1_liq_flow = 0.0
                assert tank2_liq_flow >= 0.000005



                tank2_temperature = self.tank2.states.liquid.temperature # temperatura prije odvođenja LNG K

                if (pbu1_state == False) or (tank1_vol_ratio<0.05):
                    self.pbu1.ne_radi()
                    self.tank1.update_states(-tank1_liq_flow, -tank1_vap_flow, -1, -1, dt)
                else:
                    T_liq_t1 = self.tank1.states.liquid.temperature
                    self.pbu1.update_states(T_liq_t1, p_tank1, temperature_glycol[i], self.protok_pbu)
                    pbu1_flow = self.pbu1.states.evap_lng.mol_flow
                    T_vap1 = self.pbu1.states.evap_lng.T_sat
                    self.tank1.update_states(-pbu1_flow, pbu1_flow, -1, T_vap1, dt)

                if pbu2_state == False:
                    self.pbu2.ne_radi()
                    self.tank2.update_states(-tank2_liq_flow, -tank2_vap_flow, -1, -1, dt)
                else:
                    self.pbu2.update_states(T_liq_t2, p_tank2, temperature_glycol[i], self.protok_pbu)
                    pbu2_flow = self.pbu2.states.evap_lng.mol_flow
                    tank2_vap_flow_tank = pbu2_flow
                    tank2_liq_flow_tank = tank2_liq_flow + pbu2_flow
                    T_vap2 = self.pbu2.states.evap_lng.T_sat
                    self.tank2.update_states(-tank2_liq_flow_tank, tank2_vap_flow_tank, -1, T_vap2, dt)

                self.thrvalve.calculate(tank2_temperature, p_tank2, pressure_evap)
                T_thrvalve =self.thrvalve.states.T_out
                self.pipe_pump_evap.calculate(T_thrvalve, -1, pressure_evap, -1, tank2_liq_flow)
                T_evap_in =self.pipe_pump_evap.states.T_out
                self.evaporator.update_states(T_evap_in, T_evap_in+50, pressure_evap, tank2_liq_flow, temperature_glycol[i], self.protok_evap)
                evap_T_out = self.evaporator.states.superh_lng.T_out

                valve1_T = self.valve1.states.T_out
                valve2_T = self.valve2.states.T_out
                self.mixer.calculate(evap_T_out, valve1_T, valve2_T, tank2_liq_flow, tank1_vap_flow, tank2_vap_flow, pressure_evap)
                vapor_flow = self.mixer.states.qn_tot
                T_mixer_out = self.mixer.states.T_out
                #print('T_mixer_out 2: ' + str(T_mixer_out))

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
        p1 = self.tank1.states.common.pressure
        p2 = self.tank2.states.common.pressure


        if True:
            while (p1<810000.0) and (p2<810000.0):
                T_tankV_1 = self.tank1.states.vapor.temperature
                T_tankV_2 = self.tank2.states.vapor.temperature
                it_ugas = it_ugas+1
                self.valve1.calculate(p_tank1, T_tankV_1, p_tank1, self.tank1.states.liquid.density)
                self.valve2.calculate(p_tank2, T_tankV_2, p_tank2, self.tank2.states.liquid.density)
                self.pbu1.ne_radi()
                self.pbu2.ne_radi()
                self.mixer.ne_radi()
                self.tank1.update_states(-0.0, 0.0, -1, -1, dt)
                self.tank2.update_states(-0.0, 0.0, -1, -1, dt)
                self.thrvalve.no_flow(0.0, 0.0)
                self.pipe_pump_evap.no_flow(0.0, 0.0)
                self.evaporator.no_flow(0, 0, 0)
                self.pipe_evap_superh.no_flow(0, 0)
                self.superheater.no_flow(0, 0, 0)
                self.pipe_superh_engine.no_flow(0, 0)
                if it_ugas%100 == 0:
                    print("               " + str(it_ugas) + "          " + str(p1) + ' ' + str(p2))
                p1 = self.tank1.states.common.pressure
                p2 = self.tank2.states.common.pressure
                time_holding = time_holding+dt
            self.holding_time = time_holding
        if end_success:
            time_sv = np.arange(i+pbu_counter+1+it_ugas)*10
        else:
            time_sv = np.arange(i+pbu_counter+it_ugas)*10
        self.results  = pd.DataFrame({'time': time_sv})
        self.save_results()



    def save_results(self):
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
        for k in self.pbu1.save.evap_lng:
            self.results[col] = self.pbu1.save.evap_lng[k]
            cols.append(('pbu1 lng', k))
            col = col + 1
        for k in self.pbu1.save.evap_glyc:
            self.results[col] = self.pbu1.save.evap_glyc[k]
            cols.append(('pbu1 glycol', k))
            col = col + 1
        for k in self.pbu1.save.evap_com:
            self.results[col] = self.pbu1.save.evap_com[k]
            cols.append(('pbu1 common', k))
            col = col + 1
        for k in self.pbu2.save.evap_lng:
            self.results[col] = self.pbu2.save.evap_lng[k]
            cols.append(('pbu2 lng', k))
            col = col + 1
        for k in self.pbu2.save.evap_glyc:
            self.results[col] = self.pbu2.save.evap_glyc[k]
            cols.append(('pbu2 glycol', k))
            col = col + 1
        for k in self.pbu2.save.evap_com:
            self.results[col] = self.pbu2.save.evap_com[k]
            cols.append(('pbu2 common', k))
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

