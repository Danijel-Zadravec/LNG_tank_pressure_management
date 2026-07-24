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

MIN_HEEL_FRACTION = 0.05


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
        self.startup_time_s = 0.0
        self.pressure_recovery_time_s = 0.0
        self.holding_time = None
        self.BOG_excess = 0.0
        self._reset_operation_log()

    def _reset_operation_log(self):
        """Initialise row-aligned operating data for every saved system state."""
        self.result_time_s = [0.0]
        self.operation_log = {
            'engine_demand_kmol_s': [0.0],
            'fuel_supplied_kmol_s': [0.0],
            'liquid_fuel_kmol_s': [0.0],
            'fuel_mode': ['initial'],
            'system_mode': ['initial'],
            'pbu_active': [False],
            'bog_valve_active': [False],
            'bog_removed_kmol_s': [0.0],
            'bog_used_kmol_s': [0.0],
            'bog_excess_kmol_s': [0.0],
            'cumulative_bog_excess_kmol': [0.0],
        }

    def _save_operation_row(
        self,
        dt,
        engine_demand,
        fuel_supplied,
        liquid_fuel,
        fuel_mode,
        system_mode,
        pbu_active,
        bog_valve_active,
        bog_removed,
        bog_used,
        bog_excess,
    ):
        """Append one operating-data row for one physical model time step."""
        self.result_time_s.append(self.result_time_s[-1] + float(dt))
        self.operation_log['engine_demand_kmol_s'].append(
            float(engine_demand)
        )
        self.operation_log['fuel_supplied_kmol_s'].append(
            float(fuel_supplied)
        )
        self.operation_log['liquid_fuel_kmol_s'].append(
            float(liquid_fuel)
        )
        self.operation_log['fuel_mode'].append(str(fuel_mode))
        self.operation_log['system_mode'].append(str(system_mode))
        self.operation_log['pbu_active'].append(bool(pbu_active))
        self.operation_log['bog_valve_active'].append(
            bool(bog_valve_active)
        )
        self.operation_log['bog_removed_kmol_s'].append(
            float(bog_removed)
        )
        self.operation_log['bog_used_kmol_s'].append(
            float(bog_used)
        )
        self.operation_log['bog_excess_kmol_s'].append(
            float(bog_excess)
        )
        self.operation_log['cumulative_bog_excess_kmol'].append(
            float(self.BOG_excess)
        )

    def _run_pbu_no_fuel_step(
        self,
        dt,
        glycol_temperature,
        engine_demand,
        fuel_mode,
        system_mode,
    ):
        """Advance one PBU-only step without supplying LNG to the engine."""
        p_tank = self.tank.states.common.pressure
        T_tank_vap = self.tank.states.vapor.temperature
        liquid_density = self.tank.states.liquid.density

        # Equal valve pressures impose zero BOG flow during pressure build-up.
        self.valve.calculate(
            p_tank,
            T_tank_vap,
            p_tank,
            liquid_density,
        )

        T_liq = self.tank.states.liquid.temperature
        self.pbu.update_states(
            T_liq,
            p_tank,
            glycol_temperature,
            self.protok_pbu,
        )

        if self.tank.states.liquid.vol_ratio > MIN_HEEL_FRACTION:
            pbu_flow = max(
                float(self.pbu.states.evap_lng.mol_flow),
                0.0,
            )
        else:
            pbu_flow = 0.0

        T_vap_return = self.pbu.states.evap_lng.T_sat
        self.tank.update_states(
            -pbu_flow,
            pbu_flow,
            -1,
            T_vap_return,
            dt,
        )

        p_updated = self.tank.states.common.pressure
        self.thrvalve.no_flow(T_liq, p_updated)
        self.pipe_pump_evap.no_flow(T_liq, p_updated)
        self.evaporator.no_flow(
            T_liq,
            p_updated,
            glycol_temperature,
        )
        self.pipe_evap_superh.no_flow(T_liq, p_updated)
        self.superheater.no_flow(
            T_liq,
            p_updated,
            glycol_temperature,
        )
        self.pipe_superh_engine.no_flow(T_liq, p_updated)

        self._save_operation_row(
            dt=dt,
            engine_demand=max(float(engine_demand), 0.0),
            fuel_supplied=0.0,
            liquid_fuel=0.0,
            fuel_mode=fuel_mode,
            system_mode=system_mode,
            pbu_active=True,
            bog_valve_active=False,
            bog_removed=0.0,
            bog_used=0.0,
            bog_excess=0.0,
        )

    def calculate(self, engine):

        self.BOG_excess = 0.0
        self.not_working_times = []
        self.startup_time_s = 0.0
        self.pressure_recovery_time_s = 0.0
        self._reset_operation_log()

        self.engine = engine
        times = np.asarray(self.engine.times, dtype=float)
        demand = np.asarray(self.engine.demand, dtype=float)
        temperature_glycol = np.asarray(
            self.engine.fuel_temperatures,
            dtype=float,
        )
        pressures = np.asarray(
            self.engine.fuel_pressures,
            dtype=float,
        )

        if len(times) < 2:
            raise ValueError(
                'At least two engine-profile time points are required.'
            )
        if not (
            len(times) == len(demand)
            == len(temperature_glycol)
            == len(pressures)
        ):
            raise ValueError(
                'Engine time, demand, pressure, and temperature arrays '
                'must have equal lengths.'
            )

        pbu_state = False
        bog_state = False
        hysteresis = 0.0

        # ---------------------------------------------------------------
        # One-time pre-voyage pressure build-up
        # ---------------------------------------------------------------
        positive_demand_indices = np.flatnonzero(demand > 0.0)
        if positive_demand_indices.size:
            startup_idx = int(positive_demand_indices[0])
            dt_startup = float(times[1] - times[0])
            pressure_evap_start = pressures[startup_idx] + 100.0
            startup_pressure = (
                pressure_evap_start + 150000.0 - 100.0
            )

            previous_pressure = self.tank.states.common.pressure
            stagnant_steps = 0

            while (
                self.tank.states.common.pressure < startup_pressure
                and self.tank.states.liquid.vol_ratio
                > MIN_HEEL_FRACTION
            ):
                self._run_pbu_no_fuel_step(
                    dt=dt_startup,
                    glycol_temperature=(
                        temperature_glycol[startup_idx]
                    ),
                    engine_demand=0.0,
                    fuel_mode='pre_voyage',
                    system_mode='PBU_prepressurization',
                )
                self.startup_time_s += dt_startup

                current_pressure = self.tank.states.common.pressure
                if current_pressure <= previous_pressure + 1.0e-6:
                    stagnant_steps += 1
                else:
                    stagnant_steps = 0
                previous_pressure = current_pressure

                if stagnant_steps >= 100:
                    raise RuntimeError(
                        'PBU pre-pressurization did not increase tank '
                        'pressure for 100 consecutive time steps.'
                    )

        # ---------------------------------------------------------------
        # Voyage simulation
        # ---------------------------------------------------------------
        for i in range(1, len(times)):
            pressure_evap = pressures[i - 1] + 100.0
            p_min_noflow = pressure_evap + 50000.0 - 100.0
            p_pbu_on = pressure_evap + 50000.0 + 10000.0
            p_pbu_max = pressure_evap + 150000.0 - 100.0
            p_bog_on = pressure_evap + 160000.0 - 100.0
            p_bog_off = pressure_evap + 100000.0 - 100.0

            if i % 1000 == 0:
                print(f'               {i}                    ')

            dt = float(times[i] - times[i - 1])
            p_tank = self.tank.states.common.pressure
            tank_vol_ratio = self.tank.states.liquid.vol_ratio
            engine_demand = max(float(demand[i]), 0.0)

            if tank_vol_ratio <= MIN_HEEL_FRACTION:
                break

            # If pressure becomes insufficient during the voyage, the
            # voyage clock continues. The PBU recovers pressure while the
            # corresponding LNG demand is recorded as a supply shortfall.
            if p_tank <= p_min_noflow:
                self._run_pbu_no_fuel_step(
                    dt=dt,
                    glycol_temperature=temperature_glycol[i],
                    engine_demand=engine_demand,
                    fuel_mode=(
                        'LNG'
                        if engine_demand > 0.0
                        else 'conventional'
                    ),
                    system_mode='PBU_pressure_recovery',
                )
                self.pressure_recovery_time_s += dt
                self.not_working_times.append(dt)
                pbu_state = True
                continue

            if not (
                tank_vol_ratio > MIN_HEEL_FRACTION + hysteresis
                and p_tank >= p_min_noflow
            ):
                break

            T_tank_vap = self.tank.states.vapor.temperature
            hysteresis = 0.0

            if p_tank < p_pbu_on:
                pbu_state = True
            elif p_tank > p_pbu_max:
                pbu_state = False

            if p_tank > p_bog_on:
                bog_state = True
            elif p_tank < p_bog_off:
                bog_state = False

            if bog_state:
                self.valve.calculate(
                    p_tank,
                    T_tank_vap,
                    pressure_evap,
                    self.tank.states.liquid.density,
                )
            else:
                self.valve.calculate(
                    p_tank,
                    T_tank_vap,
                    p_tank,
                    self.tank.states.liquid.density,
                )

            bog_removed_flow = max(
                float(self.valve.states.gas_mol_flow),
                0.0,
            )
            bog_used_flow = min(
                bog_removed_flow,
                engine_demand,
            )
            bog_excess_flow = max(
                bog_removed_flow - engine_demand,
                0.0,
            )

            tank_vap_flow = bog_removed_flow
            tank_liq_flow = engine_demand - bog_used_flow
            self.BOG_excess += bog_excess_flow * dt

            tank_temperature = self.tank.states.liquid.temperature
            if not pbu_state:
                self.pbu.ne_radi()
                self.tank.update_states(
                    -tank_liq_flow,
                    -tank_vap_flow,
                    -1,
                    -1,
                    dt,
                )
            else:
                self.pbu.update_states(
                    tank_temperature,
                    p_tank,
                    temperature_glycol[i],
                    self.protok_pbu,
                )
                pbu_flow = self.pbu.states.evap_lng.mol_flow
                tank_vap_flow_tank = pbu_flow - tank_vap_flow
                tank_liq_flow_tank = tank_liq_flow + pbu_flow
                T_vap = self.pbu.states.evap_lng.T_sat
                self.tank.update_states(
                    -tank_liq_flow_tank,
                    tank_vap_flow_tank,
                    -1,
                    T_vap,
                    dt,
                )

            self.thrvalve.calculate(
                tank_temperature,
                p_tank,
                pressure_evap,
            )
            T_thrvalve = self.thrvalve.states.T_out
            self.pipe_pump_evap.calculate(
                T_thrvalve,
                -1,
                pressure_evap,
                -1,
                tank_liq_flow,
            )
            T_evap_in = self.pipe_pump_evap.states.T_out
            self.evaporator.update_states(
                T_evap_in,
                T_evap_in + 50.0,
                pressure_evap,
                tank_liq_flow,
                temperature_glycol[i],
                self.protok_evap,
            )
            evap_T_out = self.evaporator.states.superh_lng.T_out
            vapor_flow = engine_demand
            T_mixer_out = evap_T_out

            self.pipe_evap_superh.calculate(
                T_mixer_out,
                -1,
                pressure_evap,
                -1,
                vapor_flow,
            )
            T_in_superh = self.pipe_evap_superh.states.T_out
            p_in_superh = self.pipe_evap_superh.states.p_out
            self.superheater.update_states(
                vapor_flow,
                T_in_superh,
                320.0,
                p_in_superh,
                self.protok_super,
                temperature_glycol[i],
            )
            T_superh_out = self.superheater.states.lng.T_out
            self.pipe_superh_engine.calculate(
                T_superh_out,
                -1,
                p_in_superh,
                -1,
                vapor_flow,
            )

            self._save_operation_row(
                dt=dt,
                engine_demand=engine_demand,
                fuel_supplied=engine_demand,
                liquid_fuel=tank_liq_flow,
                fuel_mode=(
                    'LNG'
                    if engine_demand > 0.0
                    else 'conventional'
                ),
                system_mode=(
                    'LNG_supply'
                    if engine_demand > 0.0
                    else 'conventional_fuel'
                ),
                pbu_active=pbu_state,
                bog_valve_active=(
                    bog_state and bog_removed_flow > 0.0
                ),
                bog_removed=bog_removed_flow,
                bog_used=bog_used_flow,
                bog_excess=bog_excess_flow,
            )

        self.results = pd.DataFrame({
            'time': np.asarray(
                self.result_time_s,
                dtype=float,
            )
        })
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

        for key, values in self.operation_log.items():
            if len(values) != len(self.results):
                raise RuntimeError(
                    f'PBU operation-log length mismatch for {key}: '
                    f'{len(values)} values versus '
                    f'{len(self.results)} result rows.'
                )
            self.results[col] = values
            cols.append(('operation', key))
            col = col + 1

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
