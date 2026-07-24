# -*- coding: utf-8 -*-
"""
Shared utilities for the one-at-a-time (OAT) sensitivity analysis of the
three FGSS (Fuel Gas Supply System) dual-fuel configurations:

    - PBU        (configurations/system_PBU.py)
    - Pump       (configurations/system_Pump.py)
    - Compressor (configurations/system_Compressor.py)

The swept parameters are:
    tank_ins_scale   : tinsulation thermal conductivity scale, λins​/λins,0​  -> {0.5, 1.0, 1.5}
    initial_fill     : initial liquid volume ratio, -                        -> {0.75, 0.83, 0.90}
    T_amb            : ambient temperature, K                                -> {278.15, 293.15, 308.15}
    fuel_flow_scale  : engine fuel demand scale, -                           -> {0.75, 1.00, 1.25}

For every case the following outputs are reported:
    - termination_reason (minimum_heel / MAWP / profile_end) and simulation duration
    - time to minimum heel (days), time to MAWP (days)
    - maximum tank pressure (bar), pressure at termination (bar), margin to MAWP (bar)
    - total BOG removed / BOG used as fuel / excess BOG requiring management (kg)
    - pump / compressor / total and specific mechanical energy
    - PBU, evaporator and superheater thermal energy (kWh)
    - mode-separated BOG removal for all three FGSS configurations

All of the above (except `bog_excess_kg` for the PBU system, and the final
rows themselves) are truncated to the first-occurring "termination" event
(MAWP exceedance or the 5% minimum-heel threshold): the underlying
System_*_single.calculate() loops have no pressure-relief cutoff, so they
keep running past a MAWP exceedance; results beyond that point are not
physically admissible and are excluded from the aggregated metrics.

In addition, the FULL (untruncated) time profile of every case -- time_s,
day, pressure_bar, fill_pct, liquid_flow_kg_per_h, vapor_flow_kg_per_h,
bog_excess_kg (cumulative) -- is saved to its own CSV under
results/sensitivity/timeprofiles/, for raw-data inspection / diagnosing any
individual run (see `save_time_profile`).

The configuration factories in `configurations.system_PBU`,
`configurations.system_Pump`, and `configurations.system_Compressor` are reused.
The tank insulation, initial filling, and ambient-temperature parameters are
varied through the corresponding module-level settings immediately before each
factory call. The insulation parameter scales the insulation thermal
conductivity `lmbd`, rather than the overall tank heat-transfer coefficients
`k_liq` and `k_vap`; those coefficients are re-derived using the same equations
as the configuration modules.
"""
import importlib
import shutil
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.rezimi import Rezimi  # noqa: E402

M = 16.04  # kg/kmol, molar mass of methane (LNG modeled as pure methane)
DT = 10.0  # s, fixed simulation time step (see create_engine_inputs)

HEEL_FRACTION = 0.05  # minimum liquid volume ratio (heel) before operation stops
HEEL_TOL = 1.0e-8
MAWP_BAR = 8.1  # bar(a), tank maximum allowable working pressure

RESULTS_DIR = ROOT / 'results' / 'sensitivity'
TIMEPROFILES_DIR = RESULTS_DIR / 'timeprofiles'

BOG_BALANCE_TOL_KG = 1.0


rezim = Rezimi()

SLOSH_TIMES = []


# ---------------------------------------------------------------------------
# Engine / demand profile (identical to PBU_dual.py, Pump_dual.py, Compressor_dual.py)
# ---------------------------------------------------------------------------
class Engine:
    def __init__(self, demand, pressures, temperatures, times):
        self.demand = demand  # fuel demand, kmol/s
        self.fuel_pressures = pressures  # demanded fuel pressure
        self.fuel_temperatures = temperatures  # demanded fuel temperature
        self.times = times  # calculation times, s


def create_engine_inputs(times, demands, pressures, temperatures):
    times_inp = []
    demands_inp = []
    pressures_inp = []
    temperatures_inp = []
    for i in range(1, len(times)):
        times_tmp = np.arange(times[i - 1], times[i], DT)
        len_lst = len(times_tmp)
        times_inp = times_inp + times_tmp.tolist()
        demands_inp = demands_inp + [demands[i - 1]] * len_lst
        pressures_inp = pressures_inp + [pressures[i - 1]] * len_lst
        temperatures_inp = temperatures_inp + [temperatures[i - 1]] * len_lst
    return (times_inp, demands_inp, pressures_inp, temperatures_inp)


def build_engine(fuel_flow_scale: float, max_days: float = 30.0) -> Engine:
    """Generate an extended periodic dual-fuel operating profile."""
    hour = 3600.0
    day = 24.0 * hour

    times = np.arange(0.0, max_days * day + DT, DT)

    # Repeated 48 h operating cycle:
    # 23 h port, 1 h manoeuvring, 23 h sailing, 1 h manouvering.
    cycle_phase = np.mod(times, 48.0 * hour)

    base_demand = np.select(
        [
            cycle_phase < 23.0 * hour,
            cycle_phase < 24.0 * hour,
            cycle_phase < 47.0 * hour,
        ],
        [
            rezim.port_winter,
            rezim.manouvering_winter,
            rezim.service_winter,
        ],
        default=rezim.manouvering_winter,
    )

    # Two days of conventional-fuel operation in each six-day block:
    # days 2–4, 8–10, 14–16, ...
    day_in_block = np.mod(times / day, 6.0)
    conventional_fuel = (day_in_block >= 2.0) & (day_in_block < 4.0)

    demands = np.where(conventional_fuel, 0.0, base_demand)
    demands = demands * 0.59 * fuel_flow_scale

    pressures = np.full_like(times, 600000.0)
    heating_medium_temperatures = np.full_like(times, 273.15 + 30.0)  # matches base cases' temperatures_glyc

    return Engine(
        demands,
        pressures,
        heating_medium_temperatures,
        times,
    )
# ---------------------------------------------------------------------------
# Parameter grid (OAT sensitivity, all other parameters held at base value)
# ---------------------------------------------------------------------------
PARAM_META = OrderedDict([
    ('tank_ins_scale', {'label': 'Insulation thermal conductivity, \u03bb/\u03bb$_0$ (\u2013)', 'values': [0.5, 1.0, 1.5], 'base': 1.0}),
    ('initial_fill', {'label': 'Initial filling (\u2013)', 'values': [0.75, 0.83, 0.90], 'base': 0.83}),
    ('T_amb', {'label': 'Ambient temperature (K)', 'values': [278.15, 293.15, 308.15], 'base': 293.15}),
    ('fuel_flow_scale', {'label': 'Fuel flow scaling (\u2013)', 'values': [0.75, 1.00, 1.25], 'base': 1.00}),
])

METRICS_META = OrderedDict([
    ('voyage_time_to_min_heel_days', 'Voyage time to minimum heel (days)'),
    ('max_pressure_bar', 'Maximum tank pressure (bar)'),
    ('pressure_margin_to_mawp_bar', 'Minimum pressure margin to MAWP (bar)'),
    ('bog_excess_kg', 'Excess BOG requiring management (kg)'),
    ('bog_used_as_fuel_kg', 'BOG used as fuel (kg)'),
    ('mechanical_energy_kWh', 'Mechanical (pump + compressor) energy (kWh)'),
])


# ---------------------------------------------------------------------------
# System specifications: maps a system key to its configuration module /
# factory function / fixed flow parameters (identical to the base scripts).
# ---------------------------------------------------------------------------
SYSTEM_SPECS = {
    'PBU': {
        'config_module': 'configurations.system_PBU',
        'builder_name': 'create_system_PBU',
        'extra_args': (5.0,),  # pbu_flow
        'evap_flow': 11.0,
        'super_flow': 12.0,
        'label': 'PBU',
    },
    'Pump': {
        'config_module': 'configurations.system_Pump',
        'builder_name': 'create_system_Pump',
        'extra_args': (),
        'evap_flow': 11.0,
        'super_flow': 12.0,
        'label': 'Pump',
    },
    'Compressor': {
        'config_module': 'configurations.system_Compressor',
        'builder_name': 'create_system_Compressor',
        'extra_args': (0.003,),  # comp_flow
        'evap_flow': 11.0,
        'super_flow': 12.0,
        'label': 'Compressor',
    },
}


def _lock_baseline_params(spec):
    '''Reads the (unmodified) baseline insulation conductivity and wall thickness
    once, before any monkey-patching happens, so that tank_ins_scale always
    scales the *original* conductivity regardless of how many cases already ran.'''
    if '_lambda_ins_0' not in spec:
        cfg = importlib.import_module(spec['config_module'])
        spec['_lambda_ins_0'] = cfg.lmbd
        spec['_delt'] = cfg.delt
    return spec['_lambda_ins_0'], spec['_delt']


for _spec in SYSTEM_SPECS.values():
    _lock_baseline_params(_spec)


# ---------------------------------------------------------------------------
# Simulation + metric extraction
# ---------------------------------------------------------------------------
def _get_positive_power_W(res, power_col, active_mask=None):
    """Return non-negative thermal power, zeroing inactive rows."""
    if power_col not in res.columns:
        raise KeyError(f'Missing thermal-power column: {power_col}')

    power_W = np.clip(
        res[power_col].to_numpy(dtype=float),
        0.0,
        None,
    )

    if active_mask is not None:
        active_mask = np.asarray(active_mask, dtype=bool)
        if len(active_mask) != len(power_W):
            raise RuntimeError(
                f'Thermal-power mask length mismatch for {power_col}: '
                f'{len(active_mask)} versus {len(power_W)}.'
            )
        power_W = np.where(active_mask, power_W, 0.0)

    return power_W


def _thermal_power_arrays(
    res,
    evaporator_flow_kmol_s,
    superheater_flow_kmol_s,
):
    """Extract PBU, evaporator and superheater thermal powers in W."""
    n = len(res)

    evaporator_active = (
        np.asarray(evaporator_flow_kmol_s, dtype=float) > 0.0
    )
    superheater_active = (
        np.asarray(superheater_flow_kmol_s, dtype=float) > 0.0
    )

    if ('pbu common', 'heat_flow') in res.columns:
        if ('pbu lng', 'mol_flow') not in res.columns:
            raise KeyError(
                "Missing PBU flow column: ('pbu lng', 'mol_flow')"
            )
        pbu_active = (
            res[
                ('pbu lng', 'mol_flow')
            ].to_numpy(dtype=float)
            > 0.0
        )
        pbu_power_W = _get_positive_power_W(
            res,
            ('pbu common', 'heat_flow'),
            pbu_active,
        )
    else:
        pbu_power_W = np.zeros(n)

    evaporator_evap_power_W = _get_positive_power_W(
        res,
        ('evap comm evap', 'heat_flow'),
        evaporator_active,
    )
    evaporator_superheat_power_W = _get_positive_power_W(
        res,
        ('evap comm superh', 'heat_flow'),
        evaporator_active,
    )
    evaporator_power_W = (
        evaporator_evap_power_W
        + evaporator_superheat_power_W
    )

    superheater_power_W = _get_positive_power_W(
        res,
        ('super comm', 'heat_flow'),
        superheater_active,
    )

    return {
        'pbu_power_W': pbu_power_W,
        'evaporator_evap_power_W': evaporator_evap_power_W,
        'evaporator_superheat_power_W':
            evaporator_superheat_power_W,
        'evaporator_power_W': evaporator_power_W,
        'superheater_power_W': superheater_power_W,
    }


def _align_engine_demand_to_results(
    engine,
    time_s,
    system_key,
):
    """
    Align a Pump/Compressor engine-demand profile with the saved result rows.

    Pump and Compressor simulations can stop when the minimum heel is reached,
    so their result arrays may contain only the initial prefix of the generated
    30-day engine profile. The saved result rows remain one-to-one aligned with
    that prefix.
    """
    demand_arr = np.asarray(engine.demand, dtype=float)
    engine_time_s = np.asarray(engine.times, dtype=float)
    result_time_s = np.asarray(time_s, dtype=float)
    n = len(result_time_s)

    if n > len(demand_arr) or n > len(engine_time_s):
        raise RuntimeError(
            f'{system_key}: result profile is longer than the available '
            f'engine profile: {n} result rows, '
            f'{len(demand_arr)} demand values and '
            f'{len(engine_time_s)} engine-time values.'
        )

    expected_time_s = engine_time_s[:n]
    if not np.allclose(
        result_time_s,
        expected_time_s,
        rtol=0.0,
        atol=1.0e-9,
    ):
        max_error_s = float(
            np.max(np.abs(result_time_s - expected_time_s))
        )
        raise RuntimeError(
            f'{system_key}: result times are not aligned with the '
            f'initial engine-profile segment; maximum time difference '
            f'is {max_error_s:.6g} s.'
        )

    demand_aligned = np.zeros(n)
    if n > 1:
        demand_aligned[1:] = demand_arr[1:n]

    return demand_aligned

def save_time_profile(system, engine, system_key, out_csv):
    """Save the complete, untruncated time profile of one simulation."""
    res = system.results
    time_s = res[(' ', 'time')].to_numpy(dtype=float)
    day = time_s / 86400.0
    pressure_bar = (
        res[('tank com', 'pressure')].to_numpy(dtype=float)
        / 1.0e5
    )
    fill_pct = (
        res[('tank liq', 'vol_ratio')].to_numpy(dtype=float)
        * 100.0
    )
    liquid_flow_kg_per_h = (
        -res[('tank liq', 'flow')].to_numpy(dtype=float)
        * M * 3600.0
    )
    net_tank_vapor_flow_kg_per_h = (
        -res[('tank vap', 'flow')].to_numpy(dtype=float)
        * M * 3600.0
    )
    liquid_temperature_K = (
        res[('tank liq', 'temperature')].to_numpy(dtype=float)
    )
    vapor_temperature_K = (
        res[('tank vap', 'temperature')].to_numpy(dtype=float)
    )
    saturation_temperature_K = (
        res[('tank com', 'T_sat')].to_numpy(dtype=float)
    )
    vapor_quantity_kmol = (
        res[('tank vap', 'quantity')].to_numpy(dtype=float)
    )
    liquid_quantity_kmol = (
        res[('tank liq', 'quantity')].to_numpy(dtype=float)
    )
    net_phase_change_kmol_per_step = (
        res[('tank com', 'evaporation')].to_numpy(dtype=float)
    )
    liquid_heat_gain_kJ_per_step = (
        res[('tank liq', 'surf_heat_flow')].to_numpy(dtype=float)
    )
    vapor_heat_gain_kJ_per_step = (
        res[('tank vap', 'surf_heat_flow')].to_numpy(dtype=float)
    )
    liquid_superheat_K = (
        liquid_temperature_K - saturation_temperature_K
    )

    n = len(time_s)
    time_step_s = np.diff(time_s, prepend=time_s[0])

    if system_key == 'PBU':
        required_operation_columns = [
            ('operation', 'engine_demand_kmol_s'),
            ('operation', 'fuel_supplied_kmol_s'),
            ('operation', 'liquid_fuel_kmol_s'),
            ('operation', 'fuel_mode'),
            ('operation', 'system_mode'),
            ('operation', 'pbu_active'),
            ('operation', 'bog_valve_active'),
            ('operation', 'bog_removed_kmol_s'),
            ('operation', 'bog_used_kmol_s'),
            ('operation', 'bog_excess_kmol_s'),
            ('operation', 'cumulative_bog_excess_kmol'),
        ]
        missing = [
            col for col in required_operation_columns
            if col not in res.columns
        ]
        if missing:
            raise KeyError(
                f'PBU operating-log columns are missing: {missing}'
            )

        engine_demand = res[
            ('operation', 'engine_demand_kmol_s')
        ].to_numpy(dtype=float)
        fuel_supplied = res[
            ('operation', 'fuel_supplied_kmol_s')
        ].to_numpy(dtype=float)
        evaporator_flow = res[
            ('operation', 'liquid_fuel_kmol_s')
        ].to_numpy(dtype=float)
        fuel_mode = res[
            ('operation', 'fuel_mode')
        ].astype(str).to_numpy()
        system_mode = res[
            ('operation', 'system_mode')
        ].astype(str).to_numpy()
        pbu_active = res[
            ('operation', 'pbu_active')
        ].astype(bool).to_numpy()
        bog_valve_active = res[
            ('operation', 'bog_valve_active')
        ].astype(bool).to_numpy()
        bog_removed_flow = res[
            ('operation', 'bog_removed_kmol_s')
        ].to_numpy(dtype=float)
        bog_used_flow = res[
            ('operation', 'bog_used_kmol_s')
        ].to_numpy(dtype=float)
        bog_excess_flow = res[
            ('operation', 'bog_excess_kmol_s')
        ].to_numpy(dtype=float)
        cumulative_bog_excess_kg = (
            res[
                ('operation', 'cumulative_bog_excess_kmol')
            ].to_numpy(dtype=float)
            * M
        )

    else:
        engine_demand = _align_engine_demand_to_results(
            engine,
            time_s,
            system_key,
        )
        fuel_supplied = engine_demand.copy()
        evaporator_flow = np.clip(
            -res[
                ('tank liq', 'flow')
            ].to_numpy(dtype=float),
            0.0,
            None,
        )
        fuel_mode = np.where(
            engine_demand > 0.0,
            'LNG',
            'conventional',
        )
        fuel_mode[0] = 'initial'
        system_mode = fuel_mode.copy()
        pbu_active = np.zeros(n, dtype=bool)

        if ('BOG valve 1', 'gas_mol_flow') in res.columns:
            bog_removed_flow = np.clip(
                res[
                    ('BOG valve 1', 'gas_mol_flow')
                ].to_numpy(dtype=float),
                0.0,
                None,
            )
        else:
            bog_removed_flow = np.clip(
                -res[
                    ('tank vap', 'flow')
                ].to_numpy(dtype=float),
                0.0,
                None,
            )

        bog_used_flow = np.minimum(
            bog_removed_flow,
            np.clip(engine_demand, 0.0, None),
        )
        bog_excess_flow = np.clip(
            bog_removed_flow - engine_demand,
            0.0,
            None,
        )
        bog_valve_active = bog_removed_flow > 0.0
        cumulative_bog_excess_kg = (
            np.cumsum(bog_excess_flow * time_step_s)
            * M
        )

    engine_demand_kg_per_h = engine_demand * M * 3600.0
    fuel_supplied_kg_per_h = fuel_supplied * M * 3600.0
    lng_supply_shortfall_flow = np.clip(
        engine_demand - fuel_supplied,
        0.0,
        None,
    )
    lng_supply_shortfall_kg_per_h = (
        lng_supply_shortfall_flow * M * 3600.0
    )
    liquid_fuel_kg_per_h = evaporator_flow * M * 3600.0
    conventional_fuel = fuel_mode == 'conventional'
    startup = system_mode == 'PBU_prepressurization'
    pressure_recovery = system_mode == 'PBU_pressure_recovery'
    bog_removed_kg_per_h = bog_removed_flow * M * 3600.0
    bog_used_kg_per_h = bog_used_flow * M * 3600.0
    bog_excess_kg_per_h = bog_excess_flow * M * 3600.0
    cumulative_bog_removed_kg = (
        np.cumsum(bog_removed_flow * time_step_s) * M
    )
    cumulative_bog_used_kg = (
        np.cumsum(bog_used_flow * time_step_s) * M
    )
    cumulative_requested_fuel_kg = (
        np.cumsum(engine_demand * time_step_s) * M
    )
    cumulative_supplied_fuel_kg = (
        np.cumsum(fuel_supplied * time_step_s) * M
    )
    cumulative_lng_supply_shortfall_kg = (
        np.cumsum(
            lng_supply_shortfall_flow * time_step_s
        )
        * M
    )
    cumulative_voyage_time_s = np.cumsum(
        np.where(startup, 0.0, time_step_s)
    )

    powers = _thermal_power_arrays(
        res,
        evaporator_flow_kmol_s=evaporator_flow,
        superheater_flow_kmol_s=fuel_supplied,
    )
    cumulative_pbu_thermal_energy_kWh = np.cumsum(
        powers['pbu_power_W'] * time_step_s
    ) / 3.6e6
    cumulative_evaporator_thermal_energy_kWh = np.cumsum(
        powers['evaporator_power_W'] * time_step_s
    ) / 3.6e6
    cumulative_superheater_thermal_energy_kWh = np.cumsum(
        powers['superheater_power_W'] * time_step_s
    ) / 3.6e6
    cumulative_startup_pbu_thermal_energy_kWh = np.cumsum(
        powers['pbu_power_W'] * time_step_s * startup
    ) / 3.6e6
    cumulative_operational_pbu_thermal_energy_kWh = np.cumsum(
        powers['pbu_power_W'] * time_step_s * ~startup
    ) / 3.6e6

    profile = pd.DataFrame({
        'time_s': time_s,
        'day': day,
        'time_step_s': time_step_s,
        'voyage_time_s': cumulative_voyage_time_s,
        'pressure_bar': pressure_bar,
        'fill_pct': fill_pct,
        'liquid_flow_kg_per_h': liquid_flow_kg_per_h,
        'net_tank_vapor_flow_kg_per_h':
            net_tank_vapor_flow_kg_per_h,
        'liquid_temperature_K': liquid_temperature_K,
        'vapor_temperature_K': vapor_temperature_K,
        'saturation_temperature_K': saturation_temperature_K,
        'liquid_superheat_K': liquid_superheat_K,
        'vapor_quantity_kmol': vapor_quantity_kmol,
        'liquid_quantity_kmol': liquid_quantity_kmol,
        'net_phase_change_kmol_per_step':
            net_phase_change_kmol_per_step,
        'liquid_heat_gain_kJ_per_step':
            liquid_heat_gain_kJ_per_step,
        'vapor_heat_gain_kJ_per_step':
            vapor_heat_gain_kJ_per_step,
        'engine_demand_kg_per_h': engine_demand_kg_per_h,
        'fuel_supplied_kg_per_h': fuel_supplied_kg_per_h,
        'lng_supply_shortfall_kg_per_h':
            lng_supply_shortfall_kg_per_h,
        'cumulative_requested_fuel_kg':
            cumulative_requested_fuel_kg,
        'cumulative_supplied_fuel_kg':
            cumulative_supplied_fuel_kg,
        'cumulative_lng_supply_shortfall_kg':
            cumulative_lng_supply_shortfall_kg,
        'liquid_fuel_to_evaporator_kg_per_h':
            liquid_fuel_kg_per_h,
        'fuel_mode': fuel_mode,
        'system_mode': system_mode,
        'conventional_fuel': conventional_fuel,
        'pre_voyage_startup': startup,
        'pressure_recovery': pressure_recovery,
        'pbu_active': pbu_active,
        'bog_valve_or_compressor_active': bog_valve_active,
        'bog_removed_kg_per_h': bog_removed_kg_per_h,
        'bog_used_kg_per_h': bog_used_kg_per_h,
        'bog_excess_kg_per_h': bog_excess_kg_per_h,
        'cumulative_bog_removed_kg':
            cumulative_bog_removed_kg,
        'cumulative_bog_used_kg': cumulative_bog_used_kg,
        'cumulative_bog_excess_kg':
            cumulative_bog_excess_kg,
        'pbu_thermal_power_W': powers['pbu_power_W'],
        'evaporator_evaporation_power_W':
            powers['evaporator_evap_power_W'],
        'evaporator_internal_superheat_power_W':
            powers['evaporator_superheat_power_W'],
        'evaporator_total_thermal_power_W':
            powers['evaporator_power_W'],
        'superheater_thermal_power_W':
            powers['superheater_power_W'],
        'cumulative_pbu_thermal_energy_kWh':
            cumulative_pbu_thermal_energy_kWh,
        'cumulative_startup_pbu_thermal_energy_kWh':
            cumulative_startup_pbu_thermal_energy_kWh,
        'cumulative_operational_pbu_thermal_energy_kWh':
            cumulative_operational_pbu_thermal_energy_kWh,
        'cumulative_evaporator_thermal_energy_kWh':
            cumulative_evaporator_thermal_energy_kWh,
        'cumulative_superheater_thermal_energy_kWh':
            cumulative_superheater_thermal_energy_kWh,
        'pbu_liquid_reynolds': (
            res[('pbu lng', 'Re')].to_numpy(dtype=float)
            if ('pbu lng', 'Re') in res.columns
            else np.full(n, np.nan)
        ),
    })

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    profile.to_csv(out_csv, index=False)

def run_case(system_key, params, profile_path=None):
    '''
    Builds and runs one simulation of `system_key` ('PBU', 'Pump' or 'Compressor')
    with the given parameter dict {tank_ins_scale, initial_fill, T_amb, fuel_flow_scale}
    and returns a dict of the reported metrics.

    If `profile_path` is given, the full (untruncated) time profile of the
    run is also saved there via `save_time_profile`.
    '''
    spec = SYSTEM_SPECS[system_key]
    lambda0, delt = spec['_lambda_ins_0'], spec['_delt']
    cfg = importlib.import_module(spec['config_module'])

    # monkey-patch the tank / ambient parameters read by create_system_* at call
    # time. `tank_ins_scale` scales the insulation thermal conductivity
    # (lambda), and k_liq/k_vap are re-derived with the same formula used by
    # the configuration modules (delt/lambda + convective resistances).
    lambda_ins = lambda0 * params['tank_ins_scale']
    cfg.k_liq = 1.0 / (delt / lambda_ins + 1.0 / 2000.0 + 1.0 / 25.0)
    cfg.k_vap = 1.0 / (delt / lambda_ins + 1.0 / 200.0 + 1.0 / 25.0)
    cfg.T_amb = params['T_amb']
    cfg.T_env = params['T_amb']
    cfg.liq_vol_ratio = params['initial_fill']

    builder = getattr(cfg, spec['builder_name'])
    system = builder(list(SLOSH_TIMES), spec['evap_flow'], spec['super_flow'], *spec['extra_args'])

    engine = build_engine(params['fuel_flow_scale'])
    system.calculate(engine)

    if profile_path is not None:
        save_time_profile(system, engine, system_key, Path(profile_path))

    return extract_metrics(system, engine, system_key)


def extract_metrics(system, engine, system_key):
    """Derive admissible pressure, BOG and energy metrics."""
    res = system.results
    time_s = res[(' ', 'time')].to_numpy(dtype=float)
    days = time_s / 86400.0
    time_step_s = np.diff(time_s, prepend=time_s[0])
    n = len(days)

    fill_pct = (
        res[('tank liq', 'vol_ratio')].to_numpy(dtype=float)
        * 100.0
    )
    fill_fraction = fill_pct / 100.0
    pressure_bar = (
        res[('tank com', 'pressure')].to_numpy(dtype=float)
        / 1.0e5
    )

    heel_indices = np.flatnonzero(
        fill_fraction <= HEEL_FRACTION + HEEL_TOL
    )
    heel_detected = bool(heel_indices.size)
    heel_index = (
        int(heel_indices[0])
        if heel_detected
        else None
    )

    mawp_indices = np.flatnonzero(pressure_bar >= MAWP_BAR)
    mawp_detected = bool(mawp_indices.size)
    mawp_index = (
        int(mawp_indices[0])
        if mawp_detected
        else None
    )

    if mawp_detected and (
        not heel_detected or mawp_index <= heel_index
    ):
        termination_reason = 'MAWP'
        term_idx = mawp_index
    elif heel_detected:
        termination_reason = 'minimum_heel'
        term_idx = heel_index
    else:
        termination_reason = 'profile_end'
        term_idx = n - 1

    heel_reached = termination_reason == 'minimum_heel'
    mawp_exceeded = termination_reason == 'MAWP'

    time_to_min_heel_days = (
        float(days[term_idx])
        if heel_reached
        else float('nan')
    )
    time_to_mawp_days = (
        float(days[term_idx])
        if mawp_exceeded
        else float('nan')
    )

    sl = slice(0, term_idx + 1)
    simulation_duration_days = float(days[term_idx])
    final_fill_pct = float(fill_pct[term_idx])

    max_pressure_bar = float(np.max(pressure_bar[sl]))
    pressure_at_termination_bar = float(pressure_bar[term_idx])
    pressure_margin_to_mawp_bar = float(
        MAWP_BAR - max_pressure_bar
    )

    if system_key == 'PBU':
        required_operation_columns = [
            ('operation', 'engine_demand_kmol_s'),
            ('operation', 'fuel_supplied_kmol_s'),
            ('operation', 'liquid_fuel_kmol_s'),
            ('operation', 'fuel_mode'),
            ('operation', 'system_mode'),
            ('operation', 'bog_removed_kmol_s'),
            ('operation', 'bog_used_kmol_s'),
            ('operation', 'bog_excess_kmol_s'),
        ]
        missing = [
            col for col in required_operation_columns
            if col not in res.columns
        ]
        if missing:
            raise KeyError(
                f'PBU operating-log columns are missing: {missing}'
            )

        demand_aligned = res[
            ('operation', 'engine_demand_kmol_s')
        ].to_numpy(dtype=float)
        fuel_supplied = res[
            ('operation', 'fuel_supplied_kmol_s')
        ].to_numpy(dtype=float)
        evaporator_flow = res[
            ('operation', 'liquid_fuel_kmol_s')
        ].to_numpy(dtype=float)
        fuel_mode = res[
            ('operation', 'fuel_mode')
        ].astype(str).to_numpy()
        system_mode = res[
            ('operation', 'system_mode')
        ].astype(str).to_numpy()
        bog_removed_flow = np.clip(
            res[
                ('operation', 'bog_removed_kmol_s')
            ].to_numpy(dtype=float),
            0.0,
            None,
        )
        bog_used_flow = np.clip(
            res[
                ('operation', 'bog_used_kmol_s')
            ].to_numpy(dtype=float),
            0.0,
            None,
        )
        bog_excess_flow = np.clip(
            res[
                ('operation', 'bog_excess_kmol_s')
            ].to_numpy(dtype=float),
            0.0,
            None,
        )

    else:
        demand_aligned = _align_engine_demand_to_results(
            engine,
            time_s,
            system_key,
        )
        fuel_supplied = demand_aligned.copy()
        evaporator_flow = np.clip(
            -res[
                ('tank liq', 'flow')
            ].to_numpy(dtype=float),
            0.0,
            None,
        )
        fuel_mode = np.where(
            demand_aligned > 0.0,
            'LNG',
            'conventional',
        )
        fuel_mode[0] = 'initial'
        system_mode = fuel_mode.copy()

        if ('BOG valve 1', 'gas_mol_flow') in res.columns:
            bog_removed_flow = np.clip(
                res[
                    ('BOG valve 1', 'gas_mol_flow')
                ].to_numpy(dtype=float),
                0.0,
                None,
            )
        else:
            bog_removed_flow = np.clip(
                -res[
                    ('tank vap', 'flow')
                ].to_numpy(dtype=float),
                0.0,
                None,
            )

        positive_demand = np.clip(
            demand_aligned,
            0.0,
            None,
        )
        bog_used_flow = np.minimum(
            bog_removed_flow,
            positive_demand,
        )
        bog_excess_flow = np.clip(
            bog_removed_flow - positive_demand,
            0.0,
            None,
        )

    total_bog_removed_kg = float(
        np.sum(bog_removed_flow[sl] * time_step_s[sl])
        * M
    )
    bog_used_as_fuel_kg = float(
        np.sum(bog_used_flow[sl] * time_step_s[sl])
        * M
    )
    bog_excess_kg = float(
        np.sum(bog_excess_flow[sl] * time_step_s[sl])
        * M
    )
    bog_balance_error_kg = (
        total_bog_removed_kg
        - bog_used_as_fuel_kg
        - bog_excess_kg
    )

    if abs(bog_balance_error_kg) > BOG_BALANCE_TOL_KG:
        raise RuntimeError(
            f'{system_key}: inconsistent BOG balance: '
            f'error = {bog_balance_error_kg:.2f} kg.'
        )

    lng_mask = fuel_mode == 'LNG'
    conventional_mask = fuel_mode == 'conventional'

    bog_removed_lng_kg = float(
        np.sum(
            bog_removed_flow[sl]
            * time_step_s[sl]
            * lng_mask[sl]
        )
        * M
    )
    bog_removed_conventional_kg = float(
        np.sum(
            bog_removed_flow[sl]
            * time_step_s[sl]
            * conventional_mask[sl]
        )
        * M
    )
    conventional_operation_days = float(
        np.sum(
            time_step_s[sl] * conventional_mask[sl]
        )
        / 86400.0
    )
    requested_fuel_kg = float(
        np.sum(demand_aligned[sl] * time_step_s[sl]) * M
    )

    fuel_supplied_kg = float(
        np.sum(fuel_supplied[sl] * time_step_s[sl]) * M
    )

    lng_supply_shortfall_flow = np.clip(
        demand_aligned - fuel_supplied,
        0.0,
        None,
    )
    lng_supply_shortfall_kg = float(
        np.sum(
            lng_supply_shortfall_flow[sl]
            * time_step_s[sl]
        )
        * M
    )

    if ('pump', 'H_mol_in') in res.columns:
        liq_flow = np.clip(
            -res[
                ('tank liq', 'flow')
            ].to_numpy(dtype=float),
            0.0,
            None,
        )
        dH_pump = (
            res[
                ('pump', 'H_mol_out')
            ].to_numpy(dtype=float)
            - res[
                ('pump', 'H_mol_in')
            ].to_numpy(dtype=float)
        )
        pump_energy_kWh = float(
            np.sum(
                liq_flow[sl]
                * dH_pump[sl]
                * time_step_s[sl]
            )
            / 3600.0
        )
    else:
        pump_energy_kWh = 0.0

    if ('compressor 1', 'H_mol_in') in res.columns:
        vap_flow = np.clip(
            -res[
                ('tank vap', 'flow')
            ].to_numpy(dtype=float),
            0.0,
            None,
        )
        dH_comp = (
            res[
                ('compressor 1', 'H_mol_out')
            ].to_numpy(dtype=float)
            - res[
                ('compressor 1', 'H_mol_in')
            ].to_numpy(dtype=float)
        )
        compressor_energy_kWh = float(
            np.sum(
                vap_flow[sl]
                * dH_comp[sl]
                * time_step_s[sl]
            )
            / 3600.0
        )
    else:
        compressor_energy_kWh = 0.0

    mechanical_energy_kWh = (
        pump_energy_kWh + compressor_energy_kWh
    )

    if system_key != 'PBU' and fuel_supplied_kg > 0.0:
        specific_mechanical_energy_kWh_per_t = (
            mechanical_energy_kWh
            / (fuel_supplied_kg / 1000.0)
        )
    else:
        specific_mechanical_energy_kWh_per_t = float('nan')

    powers = _thermal_power_arrays(
        res,
        evaporator_flow_kmol_s=evaporator_flow,
        superheater_flow_kmol_s=fuel_supplied,
    )
    pbu_thermal_energy_kWh = float(
        np.sum(
            powers['pbu_power_W'][sl]
            * time_step_s[sl]
        )
        / 3.6e6
    )
    evaporator_evaporation_energy_kWh = float(
        np.sum(
            powers['evaporator_evap_power_W'][sl]
            * time_step_s[sl]
        )
        / 3.6e6
    )
    evaporator_internal_superheat_energy_kWh = float(
        np.sum(
            powers['evaporator_superheat_power_W'][sl]
            * time_step_s[sl]
        )
        / 3.6e6
    )
    evaporator_thermal_energy_kWh = float(
        np.sum(
            powers['evaporator_power_W'][sl]
            * time_step_s[sl]
        )
        / 3.6e6
    )
    superheater_thermal_energy_kWh = float(
        np.sum(
            powers['superheater_power_W'][sl]
            * time_step_s[sl]
        )
        / 3.6e6
    )
    total_process_thermal_energy_kWh = (
        pbu_thermal_energy_kWh
        + evaporator_thermal_energy_kWh
        + superheater_thermal_energy_kWh
    )

    startup_mask = system_mode == 'PBU_prepressurization'
    pressure_recovery_mask = (
        system_mode == 'PBU_pressure_recovery'
    )

    startup_time_days = float(
        np.sum(time_step_s[sl] * startup_mask[sl])
        / 86400.0
    )
    pressure_recovery_time_days = float(
        np.sum(
            time_step_s[sl]
            * pressure_recovery_mask[sl]
        )
        / 86400.0
    )
    voyage_simulation_duration_days = float(
        np.sum(time_step_s[sl] * ~startup_mask[sl])
        / 86400.0
    )

    startup_pbu_thermal_energy_kWh = float(
        np.sum(
            powers['pbu_power_W'][sl]
            * time_step_s[sl]
            * startup_mask[sl]
        )
        / 3.6e6
    )
    operational_pbu_thermal_energy_kWh = (
        pbu_thermal_energy_kWh
        - startup_pbu_thermal_energy_kWh
    )

    elapsed_time_to_min_heel_days = time_to_min_heel_days
    if heel_reached:
        voyage_time_to_min_heel_days = float(
            np.sum(
                time_step_s[sl]
                * ~startup_mask[sl]
            )
            / 86400.0
        )
    else:
        voyage_time_to_min_heel_days = float('nan')

    # Retained as an alias for compatibility with existing result files.
    pbu_pressurization_time_days = startup_time_days

    return {
        'termination_reason': termination_reason,
        'heel_reached': heel_reached,
        'mawp_exceeded': mawp_exceeded,
        'results_truncated_at_mawp': mawp_exceeded,
        'simulation_duration_days': simulation_duration_days,
        'voyage_simulation_duration_days':
            voyage_simulation_duration_days,
        'time_to_min_heel_days': time_to_min_heel_days,
        'elapsed_time_to_min_heel_days':
            elapsed_time_to_min_heel_days,
        'voyage_time_to_min_heel_days':
            voyage_time_to_min_heel_days,
        'time_to_mawp_days': time_to_mawp_days,
        'final_fill_pct': final_fill_pct,
        'max_pressure_bar': max_pressure_bar,
        'pressure_at_termination_bar':
            pressure_at_termination_bar,
        'pressure_margin_to_mawp_bar':
            pressure_margin_to_mawp_bar,
        'total_bog_removed_kg': total_bog_removed_kg,
        'bog_used_as_fuel_kg': bog_used_as_fuel_kg,
        'bog_excess_kg': bog_excess_kg,
        'bog_removed_lng_kg': bog_removed_lng_kg,
        'bog_removed_conventional_kg':
            bog_removed_conventional_kg,
        'conventional_operation_days':
            conventional_operation_days,
        'requested_fuel_kg': requested_fuel_kg,
        'fuel_supplied_kg': fuel_supplied_kg,
        'lng_supply_shortfall_kg':
            lng_supply_shortfall_kg,
        'unserved_fuel_kg': lng_supply_shortfall_kg,
        'pump_energy_kWh': pump_energy_kWh,
        'compressor_energy_kWh': compressor_energy_kWh,
        'mechanical_energy_kWh': mechanical_energy_kWh,
        'specific_mechanical_energy_kWh_per_t':
            specific_mechanical_energy_kWh_per_t,
        'pbu_thermal_energy_kWh': pbu_thermal_energy_kWh,
        'startup_pbu_thermal_energy_kWh':
            startup_pbu_thermal_energy_kWh,
        'operational_pbu_thermal_energy_kWh':
            operational_pbu_thermal_energy_kWh,
        'evaporator_evaporation_energy_kWh':
            evaporator_evaporation_energy_kWh,
        'evaporator_internal_superheat_energy_kWh':
            evaporator_internal_superheat_energy_kWh,
        'evaporator_thermal_energy_kWh':
            evaporator_thermal_energy_kWh,
        'superheater_thermal_energy_kWh':
            superheater_thermal_energy_kWh,
        'total_process_thermal_energy_kWh':
            total_process_thermal_energy_kWh,
        'startup_time_days': startup_time_days,
        'pressure_recovery_time_days':
            pressure_recovery_time_days,
        'pbu_pressurization_time_days':
            pbu_pressurization_time_days,
        'bog_balance_error_kg': bog_balance_error_kg,
    }


# ---------------------------------------------------------------------------
# OAT sweep driver
# ---------------------------------------------------------------------------
def run_oat_sensitivity(system_key, out_dir=None, save_time_profiles=True):
    '''
    Runs the baseline case plus the low/high level of every parameter in
    PARAM_META (holding the other parameters at their base value), and saves
    a long-form CSV with one row per case.

    If `save_time_profiles` is True (default), the full (untruncated) time
    profile of every case (time_s, day, pressure_bar, fill_pct,
    liquid_flow_kg_per_h, vapor_flow_kg_per_h, bog_excess_kg) is also saved,
    one CSV per case, under `out_dir/timeprofiles/`. Cases whose parameter
    value equals the baseline reuse (copy) the baseline's time profile file
    instead of recomputing it, mirroring the metric-reuse optimization below.

    Returns (df, base_metrics).
    '''
    out_dir = Path(out_dir) if out_dir is not None else RESULTS_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    profiles_dir = out_dir / 'timeprofiles'
    if save_time_profiles:
        profiles_dir.mkdir(parents=True, exist_ok=True)

    base_params = {name: meta['base'] for name, meta in PARAM_META.items()}
    print(f"[{system_key}] running baseline case ...")
    base_profile_path = profiles_dir / f'{system_key}_dual_timeprofile_baseline.csv' if save_time_profiles else None
    base_metrics = run_case(system_key, base_params, profile_path=base_profile_path)

    rows = [{'parameter': 'baseline', 'level': 'base', 'value': np.nan, **base_metrics}]

    for pname, meta in PARAM_META.items():
        for val in meta['values']:
            if np.isclose(val, meta['base']):
                level = 'base'
                metrics = base_metrics  # avoid recomputing the identical baseline case
                if save_time_profiles:
                    case_profile_path = profiles_dir / f'{system_key}_dual_timeprofile_{pname}_base.csv'
                    shutil.copyfile(base_profile_path, case_profile_path)
            else:
                level = 'low' if val < meta['base'] else 'high'
                params = dict(base_params)
                params[pname] = val
                print(f"[{system_key}] running {pname} = {val} ({level}) ...")
                case_profile_path = profiles_dir / f'{system_key}_dual_timeprofile_{pname}_{level}.csv' if save_time_profiles else None
                metrics = run_case(system_key, params, profile_path=case_profile_path)

            rows.append({'parameter': pname, 'level': level, 'value': val, **metrics})

    df = pd.DataFrame(rows)
    profile_end_cases = df.loc[
        df['termination_reason'] == 'profile_end',
        ['parameter', 'level', 'value'],
    ]

    if not profile_end_cases.empty:
        print(
            f'[{system_key}] WARNING: the following cases reached '
            'the end of the generated profile before reaching the '
            'minimum heel or MAWP:'
        )
        print(profile_end_cases.to_string(index=False))

    out_csv = out_dir / f'{system_key}_dual_sensitivity_results.csv'
    df.to_csv(out_csv, index=False)
    print(f"[{system_key}] sensitivity results saved to {out_csv}")
    return df, base_metrics


# ---------------------------------------------------------------------------
# Publication-style tornado diagram (one figure per FGSS, 6 metric panels)
# ---------------------------------------------------------------------------
def plot_tornado(
    df,
    base_metrics,
    system_label,
    out_path,
    na_metrics=None,
    na_note='Not applicable for this FGSS',
):
    """Generate a six-panel OAT tornado diagram."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    na_metrics = set(na_metrics) if na_metrics else set()

    metric_meta = OrderedDict(METRICS_META)
    if system_label == 'PBU':
        metric_meta.pop('mechanical_energy_kWh', None)
        metric_meta['pbu_thermal_energy_kWh'] = (
            'PBU thermal energy (kWh)'
        )

    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 10,
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'figure.dpi': 150,
    })

    color_low = '#4C72B0'
    color_high = '#DD8452'
    mawp_hatch = '///'

    metric_keys = list(metric_meta.keys())
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(16, 8),
        constrained_layout=True,
    )
    axes = axes.ravel()

    any_mawp = bool(
        (df['termination_reason'] == 'MAWP').any()
    )

    for ax, metric_key in zip(axes, metric_keys):
        metric_label = metric_meta[metric_key]
        ax.set_title(metric_label)
        ax.set_xlabel(metric_label)

        if metric_key in na_metrics:
            ax.text(
                0.5,
                0.5,
                na_note,
                transform=ax.transAxes,
                ha='center',
                va='center',
                fontsize=10,
                color='gray',
                style='italic',
                wrap=True,
            )
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            continue

        base_val = base_metrics[metric_key]

        bars = []
        for pname, meta in PARAM_META.items():
            sub = df[df['parameter'] == pname]
            low_row = sub.loc[sub['level'] == 'low'].iloc[0]
            high_row = sub.loc[sub['level'] == 'high'].iloc[0]
            bars.append({
                'label': meta['label'],
                'low': low_row[metric_key],
                'high': high_row[metric_key],
                'low_reason': low_row['termination_reason'],
                'high_reason': high_row['termination_reason'],
                'low_mawp_time': low_row['time_to_mawp_days'],
                'high_mawp_time': high_row['time_to_mawp_days'],
            })

        def sensitivity_range(row):
            difference = row['high'] - row['low']
            return (
                abs(difference)
                if np.isfinite(difference)
                else -1.0
            )

        bars.sort(key=sensitivity_range)
        labels = [row['label'] for row in bars]
        y = np.arange(len(labels))

        for yi, row in zip(y, bars):
            notes = []

            for level, value, color, reason, mawp_time in [
                (
                    'Low',
                    row['low'],
                    color_low,
                    row['low_reason'],
                    row['low_mawp_time'],
                ),
                (
                    'High',
                    row['high'],
                    color_high,
                    row['high_reason'],
                    row['high_mawp_time'],
                ),
            ]:
                if np.isfinite(base_val) and np.isfinite(value):
                    ax.barh(
                        yi,
                        value - base_val,
                        left=base_val,
                        color=color,
                        edgecolor='black',
                        height=0.6,
                        zorder=3,
                        hatch=(
                            mawp_hatch
                            if reason == 'MAWP'
                            else None
                        ),
                    )
                elif reason == 'MAWP':
                    if np.isfinite(mawp_time):
                        notes.append(
                            f'{level}: MAWP at {mawp_time:.2f} d'
                        )
                    else:
                        notes.append(f'{level}: MAWP')
                else:
                    notes.append(
                        f'{level}: not reached before termination'
                    )

            if notes:
                ax.text(
                    0.02,
                    yi,
                    '; '.join(notes),
                    transform=ax.get_yaxis_transform(),
                    ha='left',
                    va='center',
                    fontsize=7.5,
                    color='gray',
                    style='italic',
                    zorder=5,
                )

        if np.isfinite(base_val):
            ax.axvline(
                base_val,
                color='black',
                linewidth=1.0,
                zorder=4,
            )

        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.grid(
            axis='x',
            linestyle='--',
            alpha=0.5,
            zorder=0,
        )
        ax.margins(y=0.15)

    handles = [
        Patch(
            facecolor=color_low,
            edgecolor='black',
            label='Low level',
        ),
        Patch(
            facecolor=color_high,
            edgecolor='black',
            label='High level',
        ),
        Line2D(
            [0],
            [0],
            color='black',
            linewidth=1.0,
            label='Baseline',
        ),
    ]

    if any_mawp:
        handles.append(
            Patch(
                facecolor='white',
                edgecolor='black',
                hatch=mawp_hatch,
                label='Value evaluated up to MAWP termination',
            )
        )

    fig.legend(
        handles=handles,
        loc='lower center',
        ncol=len(handles),
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle(
        f'One-at-a-time sensitivity — '
        f'{system_label} FGSS (dual fuel)',
        fontsize=13,
    )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out_path,
        dpi=300,
        bbox_inches='tight',
    )
    plt.close(fig)
    print(
        f'[{system_label}] tornado diagram saved to {out_path}'
    )

