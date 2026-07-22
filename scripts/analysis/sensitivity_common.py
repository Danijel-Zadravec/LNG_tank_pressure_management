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
    - pump / compressor / total mechanical energy (kWh)

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

The three base cases (PBU_dual.py, Pump_dual.py, Compressor_dual.py) are
reused unmodified: the tank insulation / initial filling / ambient
temperature parameters are varied by monkey-patching the corresponding
module-level globals in the `configurations.system_*` modules right before
calling their `create_system_*` factory function (these functions read the
parameters from their module's globals at call time), so the underlying
physical model stays identical to the base cases. The insulation parameter
scales the insulation thermal conductivity `lmbd`, not the overall tank
heat-transfer coefficients `k_liq`/`k_vap` directly (those are re-derived
from the scaled conductivity using the same formula as the configuration
modules).
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
    ('time_to_min_heel_days', 'Time to minimum heel (days)'),
    ('max_pressure_bar', 'Maximum tank pressure (bar)'),
    ('pressure_at_termination_bar', 'Pressure at termination (bar)'),
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
def save_time_profile(system, engine, system_key, out_csv):
    '''
    Saves the FULL (untruncated) time profile of a finished system.calculate()
    run to `out_csv`, with columns:
        time_s, day, pressure_bar, fill_pct,
        liquid_flow_kg_per_h, vapor_flow_kg_per_h, bog_excess_kg

    Sign convention: liquid_flow_kg_per_h / vapor_flow_kg_per_h are positive
    when LNG/BOG is being withdrawn from the tank (i.e. the negative of the
    raw ('tank liq'|'tank vap','flow') columns, which are positive when mass
    is added to the tank).

    `bog_excess_kg` is the CUMULATIVE excess BOG (kg) that could not be
    absorbed as fuel, reconstructed at every row from BOG-removed vs. engine
    demand (same criterion the model itself uses for system.BOG_excess). This
    exact reconstruction relies on results rows being 1:1 aligned with
    engine.demand, which holds for Pump/Compressor but NOT for PBU (its
    internal low-pressure pressurization-stall loop inserts extra rows -- see
    `not_working_times`); for PBU this column is therefore reported as NaN.
    '''
    res = system.results
    time_s = res[(' ', 'time')].to_numpy(dtype=float)
    day = time_s / 86400.0
    pressure_bar = res[('tank com', 'pressure')].to_numpy(dtype=float) / 1.0e5
    fill_pct = res[('tank liq', 'vol_ratio')].to_numpy(dtype=float) * 100.0
    liquid_flow_kg_per_h = -res[('tank liq', 'flow')].to_numpy(dtype=float) * M * 3600.0
    vapor_flow_kg_per_h = -res[('tank vap', 'flow')].to_numpy(dtype=float) * M * 3600.0
    liquid_temperature_K = res[('tank liq', 'temperature')].to_numpy(dtype=float)
    vapor_temperature_K = res[('tank vap', 'temperature')].to_numpy(dtype=float)
    saturation_temperature_K = res[('tank com', 'T_sat')].to_numpy(dtype=float)
    vapor_quantity_kmol = res[('tank vap', 'quantity')].to_numpy(dtype=float)
    liquid_quantity_kmol = res[('tank liq', 'quantity')].to_numpy(dtype=float)
    evaporation_per_step_kmol = res[('tank com', 'evaporation')].to_numpy(dtype=float)
    liquid_heat_gain_kJ = res[('tank liq', 'surf_heat_flow')].to_numpy(dtype=float)
    vapor_heat_gain_kJ = res[('tank vap', 'surf_heat_flow')].to_numpy(dtype=float)
    liquid_superheat_K = (liquid_temperature_K - saturation_temperature_K)
    n = len(time_s)
    if system_key == 'PBU':
        bog_excess_kg = np.full(n, np.nan)
    else:
        if ('BOG valve 1', 'gas_mol_flow') in res.columns:
            bog_removed_flow = np.clip(res[('BOG valve 1', 'gas_mol_flow')].to_numpy(dtype=float), 0.0, None)
        else:
            bog_removed_flow = np.clip(-res[('tank vap', 'flow')].to_numpy(dtype=float), 0.0, None)
        demand_arr = np.asarray(engine.demand, dtype=float)
        demand_aligned = np.zeros(n)
        demand_aligned[1:n] = demand_arr[1:n]
        if len(demand_arr) != n:
            raise RuntimeError(
                f'{system_key}: results and engine demand are not aligned: '
                f'{n} result rows versus {len(demand_arr)} demand values.'
            )
        engine_demand_kg_per_h = demand_aligned * M * 3600.0
        conventional_fuel = demand_aligned <= 0.0
        bog_removed_kg_per_h = bog_removed_flow * M * 3600.0
        cumulative_bog_removed_kg = (
            np.cumsum(bog_removed_flow) * DT * M
        )
        excess_flow = np.clip(bog_removed_flow - demand_aligned, 0.0, None)
        bog_excess_kg = np.cumsum(excess_flow) * DT * M

    profile = pd.DataFrame({
        'time_s': time_s,
        'day': day,
        'pressure_bar': pressure_bar,
        'fill_pct': fill_pct,
        'liquid_flow_kg_per_h': liquid_flow_kg_per_h,
        'vapor_flow_kg_per_h': vapor_flow_kg_per_h,
        'bog_excess_kg': bog_excess_kg,
        'liquid_temperature_K': liquid_temperature_K,
        'vapor_temperature_K': vapor_temperature_K,
        'saturation_temperature_K': saturation_temperature_K,
        'vapor_quantity_kmol': vapor_quantity_kmol,
        'liquid_quantity_kmol': liquid_quantity_kmol,
        'evaporation_per_step_kmol': evaporation_per_step_kmol,
        'liquid_heat_gain_kJ': liquid_heat_gain_kJ,
        'vapor_heat_gain_kJ': vapor_heat_gain_kJ,
        'liquid_superheat_K': liquid_superheat_K,
        'engine_demand_kg_per_h': engine_demand_kg_per_h,
        'conventional_fuel': conventional_fuel,
        'bog_removed_kg_per_h': bog_removed_kg_per_h,
        'cumulative_bog_removed_kg': cumulative_bog_removed_kg,
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
    '''
    Derives the reported metrics from a finished system.calculate() run.

    Metrics that are time-indexed series (pressure, BOG flow, pump/compressor
    enthalpy rise) are truncated to the first-occurring "termination" event
    (MAWP exceedance or the 5% minimum-heel threshold): the underlying
    System_*_single.calculate() loops have no pressure-relief cutoff, so the
    raw arrays keep going past a MAWP exceedance / do not stop exactly at the
    heel threshold row; results beyond the termination point are not
    physically admissible and are excluded here.

    Exception: `bog_excess_kg` for the PBU system is `system.BOG_excess`, a
    *scalar* accumulated by the model over the entire simulated profile. PBU
    inserts extra pressurization sub-steps (see `not_working_times`) that
    break the 1:1 row <-> engine.demand index alignment that would be needed
    to reconstruct a truncated, time-resolved value, so for PBU it is
    reported over the full simulated profile, not just up to the reported
    termination point (for Pump/Compressor, where that alignment holds, the
    truncated value is reconstructed exactly from BOG-removed vs. demand).
    '''
    res = system.results
    time_s = res[(' ', 'time')].to_numpy(dtype=float)
    days = time_s / 86400.0
    n = len(days)

    fill_pct = res[('tank liq', 'vol_ratio')].to_numpy(dtype=float) * 100.0
    fill_fraction = fill_pct / 100.0
    pressure_bar = res[('tank com', 'pressure')].to_numpy(dtype=float) / 1.0e5

    # --- time to minimum heel: first threshold crossing, NOT argmin over the
    #     record (fill ratio can plateau/creep up during conventional-fuel
    #     operation, so the smallest in-record value is not necessarily the
    #     heel limit being reached) ---
    heel_indices = np.flatnonzero(fill_fraction <= HEEL_FRACTION + HEEL_TOL)
    heel_reached = bool(heel_indices.size)
    heel_detected = bool(heel_indices.size)
    heel_index = int(heel_indices[0]) if heel_detected else None
    time_to_min_heel_days = float(days[heel_index]) if heel_reached else float('nan')

    # --- MAWP threshold crossing ---
    mawp_indices = np.flatnonzero(pressure_bar >= MAWP_BAR)
    mawp_exceeded = bool(mawp_indices.size)
    mawp_detected = bool(mawp_indices.size)
    mawp_index = int(mawp_indices[0]) if mawp_detected else None
    time_to_mawp_days = float(days[mawp_index]) if mawp_exceeded else float('nan')

    # --- termination reason: whichever of {MAWP, minimum_heel} occurs first;
    #     if neither occurs within the simulated profile -> profile_end ---
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

    sl = slice(0, term_idx + 1)  # admissible window: up to (and incl.) termination
    simulation_duration_days = float(days[term_idx])
    final_fill_pct = float(fill_pct[term_idx])

    max_pressure_bar = float(np.max(pressure_bar[sl]))
    pressure_at_termination_bar = float(pressure_bar[term_idx])
    pressure_margin_to_mawp_bar = float(MAWP_BAR - max_pressure_bar)

    # --- BOG removed from the tank (admissible window only) ---
    if ('BOG valve 1', 'gas_mol_flow') in res.columns:
        # PBU / Pump: BOG offtake goes through a dedicated BOG relief valve
        bog_removed_flow = np.clip(res[('BOG valve 1', 'gas_mol_flow')].to_numpy(dtype=float), 0.0, None)
    else:
        # Compressor: BOG is drawn directly from the tank vapor space
        bog_removed_flow = np.clip(-res[('tank vap', 'flow')].to_numpy(dtype=float), 0.0, None)
    total_bog_removed_kg = float(np.sum(bog_removed_flow[sl])) * DT * M

    if system_key == 'PBU':
        # PBU excess BOG is available only as a full-profile scalar.
        # Do not report it if the admissible simulation terminates at MAWP.
        if termination_reason == 'MAWP':
            bog_excess_kg = float('nan')
        else:
            bog_excess_kg = float(system.BOG_excess) * M
    else:
        # exact reconstruction: rows are 1:1 aligned with engine.demand for
        # Pump/Compressor (no extra pressurization sub-steps), so we can
        # reproduce system.BOG_excess's own criterion (demand < BOG removed)
        # on the truncated window.
        demand_arr = np.asarray(engine.demand, dtype=float)
        demand_aligned = np.zeros(n)
        demand_aligned[1:n] = demand_arr[1:n]
        excess_flow = np.clip(bog_removed_flow - demand_aligned, 0.0, None)
        bog_excess_kg = float(np.sum(excess_flow[sl])) * DT * M

    if np.isfinite(bog_excess_kg):
        raw_bog_used_kg = total_bog_removed_kg - bog_excess_kg

        if raw_bog_used_kg < -BOG_BALANCE_TOL_KG:
            raise RuntimeError(
                f'{system_key}: inconsistent BOG balance: '
                f'total removed = {total_bog_removed_kg:.2f} kg, '
                f'excess = {bog_excess_kg:.2f} kg.'
            )

        bog_used_as_fuel_kg = max(raw_bog_used_kg, 0.0)
        bog_balance_error_kg = (
            total_bog_removed_kg
            - bog_used_as_fuel_kg
            - bog_excess_kg
        )
    else:
        bog_used_as_fuel_kg = float('nan')
        bog_balance_error_kg = float('nan')
            
    # --- Pump energy (kWh); zero for PBU (no mechanical liquid pump) ---
    if ('pump', 'H_mol_in') in res.columns:
        liq_flow = np.clip(-res[('tank liq', 'flow')].to_numpy(dtype=float), 0.0, None)  # kmol/s
        dH_pump = (res[('pump', 'H_mol_out')].to_numpy(dtype=float) - res[('pump', 'H_mol_in')].to_numpy(dtype=float))  # kJ/kmol
        pump_energy_kWh = float(np.sum((liq_flow * dH_pump)[sl])) * DT / 3600.0
    else:
        pump_energy_kWh = 0.0

    # --- Compressor energy (kWh); zero for PBU and Pump (no compressor) ---
    if ('compressor 1', 'H_mol_in') in res.columns:
        vap_flow = np.clip(-res[('tank vap', 'flow')].to_numpy(dtype=float), 0.0, None)  # kmol/s
        dH_comp = (res[('compressor 1', 'H_mol_out')].to_numpy(dtype=float) - res[('compressor 1', 'H_mol_in')].to_numpy(dtype=float))  # kJ/kmol
        compressor_energy_kWh = float(np.sum((vap_flow * dH_comp)[sl])) * DT / 3600.0
    else:
        compressor_energy_kWh = 0.0

    mechanical_energy_kWh = pump_energy_kWh + compressor_energy_kWh

    # PBU-only: total time spent in the low-pressure pressurization stall (the
    # engine receives no fuel during this stall; see the `while` loop in
    # System_Pbu_single.calculate). NaN for Pump/Compressor (no such stall).
    if hasattr(system, 'not_working_times'):
        pbu_pressurization_time_days = float(np.sum(system.not_working_times)) / 86400.0
    else:
        pbu_pressurization_time_days = float('nan')
    if np.isfinite(time_to_min_heel_days):
        voyage_time_to_min_heel_days = time_to_min_heel_days

        if (
            system_key == 'PBU'
            and np.isfinite(pbu_pressurization_time_days)
        ):
            voyage_time_to_min_heel_days = max(
                time_to_min_heel_days
                - pbu_pressurization_time_days,
                0.0,
            )
    else:
        voyage_time_to_min_heel_days = float('nan')

    return {
        'termination_reason': termination_reason,
        'heel_reached': heel_reached,
        'mawp_exceeded': mawp_exceeded,
        'simulation_duration_days': simulation_duration_days,
        'time_to_min_heel_days': time_to_min_heel_days,
        'voyage_time_to_min_heel_days': voyage_time_to_min_heel_days,
        'time_to_mawp_days': time_to_mawp_days,
        'final_fill_pct': final_fill_pct,
        'max_pressure_bar': max_pressure_bar,
        'pressure_at_termination_bar': pressure_at_termination_bar,
        'pressure_margin_to_mawp_bar': pressure_margin_to_mawp_bar,
        'total_bog_removed_kg': total_bog_removed_kg,
        'bog_used_as_fuel_kg': bog_used_as_fuel_kg,
        'bog_excess_kg': bog_excess_kg,
        'pump_energy_kWh': pump_energy_kWh,
        'compressor_energy_kWh': compressor_energy_kWh,
        'mechanical_energy_kWh': mechanical_energy_kWh,
        'pbu_pressurization_time_days': pbu_pressurization_time_days,
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
def plot_tornado(df, base_metrics, system_label, out_path, na_metrics=None, na_note='Not applicable for this FGSS'):
    '''
    `na_metrics`: optional set of METRICS_META keys that are structurally not
    applicable for this system (e.g. `mechanical_energy_kWh` for PBU, which
    has no rotating machinery and instead consumes thermal energy from the
    glycol heating medium): the panel is drawn empty with an explanatory note
    instead of a row of zero-width bars.
    '''
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    na_metrics = set(na_metrics) if na_metrics else set()

    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 10,
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'figure.dpi': 150,
    })

    color_low = '#4C72B0'
    color_high = '#DD8452'

    metric_keys = list(METRICS_META.keys())
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.ravel()

    for ax, metric_key in zip(axes, metric_keys):
        ax.set_title(METRICS_META[metric_key])
        ax.set_xlabel(METRICS_META[metric_key])

        if metric_key in na_metrics:
            ax.text(0.5, 0.5, na_note, transform=ax.transAxes, ha='center', va='center',
                    fontsize=10, color='gray', style='italic', wrap=True)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            continue

        base_val = base_metrics[metric_key]

        bars = []
        for pname, meta in PARAM_META.items():
            sub = df[df['parameter'] == pname]
            low_val = sub.loc[sub['level'] == 'low', metric_key].iloc[0]
            high_val = sub.loc[sub['level'] == 'high', metric_key].iloc[0]
            bars.append((meta['label'], low_val, high_val))

        # sort by sensitivity range (largest at top); NaN ranges (e.g. heel
        # not reached within the simulated profile) sink to the bottom
        bars.sort(key=lambda r: abs(r[2] - r[1]) if np.isfinite(r[2] - r[1]) else -1.0)
        labels = [b[0] for b in bars]
        lows = [b[1] for b in bars]
        highs = [b[2] for b in bars]

        y = np.arange(len(labels))
        for yi, lo, hi in zip(y, lows, highs):
            if not (np.isfinite(base_val) and np.isfinite(lo) and np.isfinite(hi)):
                ax.text(0.02, yi, 'not reached before termination', transform=ax.get_yaxis_transform(),
                        ha='left', va='center', fontsize=8, color='gray', style='italic', zorder=3)
                continue
            ax.barh(yi, lo - base_val, left=base_val, color=color_low, edgecolor='black', height=0.6, zorder=3)
            ax.barh(yi, hi - base_val, left=base_val, color=color_high, edgecolor='black', height=0.6, zorder=3)

        if np.isfinite(base_val):
            ax.axvline(base_val, color='black', linewidth=1.0, zorder=4)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.grid(axis='x', linestyle='--', alpha=0.5, zorder=0)
        ax.margins(y=0.15)

    handles = [
        Patch(facecolor=color_low, edgecolor='black', label='Low level'),
        Patch(facecolor=color_high, edgecolor='black', label='High level'),
        Line2D([0], [0], color='black', linewidth=1.0, label='Baseline'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(f'One-at-a-time sensitivity \u2014 {system_label} FGSS (dual fuel)', fontsize=13)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"[{system_label}] tornado diagram saved to {out_path}")
