# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy

from src.tank import Tank_params, Tank_initial, Tank
from src.Pump import PumpStates, LngPump
from src.Evaporator import EVparams, EV_lng, EV_glycol, StatesCommonEV, EV_States, Evaporator
from src.Superheater import HEparams, LngStream, GlycolStream, StatesHE, HeatExchanger, StatesCommonHE
from src.Pipes import PipeParams, PipeStates, Pipe
import src.properties.liq_props as liq

class Engine:
    def __init__(self, demand, pressures, temperatures, times):
        self.demand = demand #fuel demand, kmol/s
        self.fuel_pressures = pressures #demanded fuel pressure
        self.fuel_temperatures = temperatures #demanded fuel pressure
        self.times = times #calculation times, s


class SystemLNG:
    def __init__(self, tank1, tank2, pump, pipe_pump_evaporator, evaporator, pipe_evaporator_superheater, superheater, pipe_superheater_engine, engine):
        self.tank1 = tank1
        self.tank2 = tank2
        self.pump = pump
        self.pipe_pump_evap = pipe_pump_evaporator
        self.evaporator = evaporator
        self.pipe_evap_superh = pipe_evaporator_superheater
        self.superheater = superheater
        self.pipe_superh_engine = pipe_superheater_engine
        self.engine = engine
        self.results = pd.DataFrame({'time': self.engine.times})

    def calculate(self):
        times = self.engine.times
        #print(times)
        demand = self.engine.demand
        temperature = self.engine.fuel_temperatures
        pressures = self.engine.fuel_pressures
        pump_pressure = pressures[0] + 100 #regulacija!
        vap_flow1 = 0.0
        vap_flow2 = 0.0

        for i in range(1,len(times)):
            #print(i)
            if i%1000 == 0:
                print("               " + str(i) + "                    ")
            dt = times[i] - times[i-1]
            if ((tank1.states.liquid.quantity < 0.1) or (tank2.states.liquid.quantity < 0.1)):
                self.results =  pd.DataFrame({'time': times[0:i]})
                break
            if self.tank1.states.common.pressure > 550000.0:
                vap_flow1 = -0.005
            elif self.tank1.states.common.pressure < 300000.0:
                vap_flow1 = 0.0
            if self.tank2.states.common.pressure > 550000.0:
                vap_flow2 = -0.005
            elif self.tank2.states.common.pressure < 300000.0:
                vap_flow2 = 0.0
            tank1_liq_flow = -demand[i]/2 - vap_flow1  #dodatna kapljevina za održavanje tlaka
            tank2_liq_flow = -demand[i]/2 - vap_flow2  #dodatna kapljevina za održavanje tlaka

            tank_vap_temperature = -1
            tank1_temperature = self.tank1.states.liquid.temperature # temperatura prije odvođenja LNG K
            tank2_temperature = self.tank2.states.liquid.temperature # temperatura prije odvođenja LNG K

            tank1_pressure = self.tank1.states.common.pressure
            tank2_pressure = self.tank2.states.common.pressure

            #print('tnk')
            self.tank1.update_states(tank1_liq_flow, vap_flow1, -1, tank_vap_temperature, dt)
            self.tank2.update_states(tank2_liq_flow, vap_flow2, -1, tank_vap_temperature, dt)
            pump_T_in = (tank1_temperature + tank2_temperature) / 2.0
            pump_p_in = (tank1_pressure + tank2_pressure) / 2
            self.pump.calculate(pump_T_in, pump_p_in, pump_pressure)
            pump_exit_temperature = self.pump.states.T_out
            #print('pp1')
            total_liq_flow =- (tank1_liq_flow + tank2_liq_flow)
            self.pipe_pump_evap.calculate(pump_exit_temperature, -1, pump_pressure, -1, total_liq_flow)
            evap_in_T = self.pipe_pump_evap.states.T_out
            evap_p_in = self.pipe_pump_evap.states.p_out
            #print('evap')
            self.evaporator.update_states(total_liq_flow, evap_in_T, evap_p_in, 323.15)
            evap_T_out = self.evaporator.states.lng.T_sat
            #print('pp2')
            vapor_flow = total_liq_flow - (vap_flow1 + vap_flow2)
            self.pipe_evap_superh.calculate(evap_T_out, -1, evap_p_in, -1, vapor_flow)
            T_in_superh = self.pipe_evap_superh.states.T_out
            p_in_superh = self.pipe_evap_superh.states.p_out
            T_superh_out = 0.1 + temperature[i]
            #print('superh')
            self.superheater.update_states(vapor_flow, T_in_superh, T_superh_out, 323.15, p_in_superh)
            #vap_return = vapor_flow - demand[i]
            #print('pp3')
            self.pipe_superh_engine.calculate(T_superh_out, -1, p_in_superh, -1, vapor_flow)
        self.results = self.results[self.results.index < (i+1)]
        self.save_results()

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
        for k in self.pipe_pump_evap.save.states:
            self.results[col] = self.pipe_pump_evap.save.states[k]
            cols.append(('p. pmp evap', k))
            col = col + 1
        for k in self.evaporator.save.lng:
            self.results[col] = self.evaporator.save.lng[k]
            cols.append(('evap lng', k))
            col = col + 1
        for k in self.evaporator.save.glycol:
            self.results[col] = self.evaporator.save.glycol[k]
            cols.append(('evap glycol', k))
            col = col + 1
        for k in self.evaporator.save.common:
            self.results[col] = self.evaporator.save.common[k]
            cols.append(('evap comm', k))
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
