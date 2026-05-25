# -*- coding: utf-8 -*-
"""
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


class SystemLNG_Pump:
    def __init__(self, tank, pump, thrvalve, pipe_pump_evaporator, evaporator, valve, pipe_evaporator_superheater, superheater, pipe_superheater_engine, evap_flow, super_flow):
        self.tank = tank
        self.pump = pump
        self.thrvalve = thrvalve
        self.pipe_pump_evap = pipe_pump_evaporator
        self.evaporator = evaporator
        self.valve = valve #
        self.pipe_evap_superh = pipe_evaporator_superheater
        self.superheater = superheater

        self.pipe_superh_engine = pipe_superheater_engine
        self.engine = None
        self.results = None
        self.time = 0
        self.protok_evap = evap_flow
        self.protok_super = super_flow
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

        vap_flow = 0.0

        h24s = 24*60*60 #sekunde u 24 dana
        h48s = 48*60*60 #sekunde u 48 dana


        for i in range(1,len(times)):
            #print(i)
            pressure_evap = pressures[i-1] + 100 #regulacija!
            pressure_bog_valve = pressure_evap + 50000.0

            if i%1000 == 0:
                print("               " + str(i) + "                    ")
            dt = times[i] - times[i-1]
            p_tank = self.tank.states.common.pressure
            T_tankV = self.tank.states.vapor.temperature

            #print(tank_vap_flow)
            end_success = True
            self.time = self.time + dt

            tank_day = (((self.time % h48s) - h24s) < 0)


            if self.tank.states.liquid.vol_ratio > 0.05:
                self.valve.calculate(p_tank, T_tankV, pressure_bog_valve, self.tank.states.liquid.density)
                tank_vap_flow = self.valve.states.gas_mol_flow #kmol/s

                tank_liq_flow = demand[i] - tank_vap_flow


                if tank_liq_flow < 0.0:
                    self.BOG_excess = self.BOG_excess - tank_liq_flow*dt
                    tank_liq_flow = 0.0


                #print(tank_liq_flow)
                tank_temperature = self.tank.states.liquid.temperature # temperatura prije odvođenja LNG K
                self.tank.update_states(-tank_liq_flow, -tank_vap_flow, -1, -1, dt)
                #if i == 1999 or i==2000 or i==2001:
                #    print([self.tank.states.common.pressure, self.tank2.states.common.pressure])
                pump_T_in = tank_temperature
                tank_pressure =  self.tank.states.common.pressure
                pump_pressure = max(pressure_evap+50000.0, tank_pressure + 50000.0)
                self.pump.calculate(pump_T_in, tank_pressure, pump_pressure)
                pump_exit_temperature = self.pump.states.T_out
                self.thrvalve.calculate(pump_exit_temperature, pump_pressure, pressure_evap)
                T_thrvalve =self.thrvalve.states.T_out
                self.pipe_pump_evap.calculate(T_thrvalve, -1, pressure_evap, -1, tank_liq_flow)
                T_evap_in =self.pipe_pump_evap.states.T_out
                self.evaporator.update_states(T_evap_in, T_evap_in+50, pressure_evap, tank_liq_flow, temperature_glycol[i],  self.protok_evap)

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
            else:
                end_success = False
                break


        it_ugas = 0
        p = self.tank.states.common.pressure


        if end_success:
            time_sv = np.arange(i+1+it_ugas)*10
        else:
            time_sv = np.arange(i+it_ugas)*10

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
