#!/usr/bin/env python3
"""
Run the final FGSS base-case analyses used in Section 5.3.

The script performs the following tasks:

1. Re-runs the LNG-only base case for:
       - FGSS 1 (PBU)
       - FGSS 2 (Pump)
       - FGSS 3 (Compressor)

2. Saves, for every re-run:
       - the comprehensive time profile produced by sensitivity_common.py;
       - the complete raw system.results table with flattened column names;
       - one summary CSV containing all extracted performance indicators.

3. Re-runs the three dual-fuel base cases over a fixed 10-day voyage window.

4. Creates combined CSV files for manuscript revision:
       - LNG_base_case_summary.csv
       - dual_base_case_summary.csv
       - all_base_case_summary.csv
       - manuscript_base_case_table.csv
       - base_case_run_manifest.json

The LNG simulations use a fixed six-day voyage window and the dual-fuel
simulations use a fixed ten-day voyage window. For FGSS 1, pre-voyage PBU
pressurization is performed before the voyage clock starts; therefore, elapsed
simulation time may exceed the prescribed voyage duration by the startup time.
FGSS 1 terminates earlier if the tank reaches the 5% minimum heel.

Expected project structure
--------------------------
project/
├── configurations/
│   ├── system_PBU.py
│   ├── system_Pump.py
│   └── system_Compressor.py
├── src/
│   ├── System_Pbu_single.py
│   ├── System_Pump_single.py
│   ├── System_Compressor_single.py
│   └── ...
├── scripts/
│   └── analysis/
│       ├── sensitivity_common.py
│       └── run_base_case_analysis.py
└── results/
    ├── sensitivity/
    └── base_case/

The single-tank PBU configuration is created through
configurations.system_PBU, which imports src.System_Pbu_single.
"""

from __future__ import annotations

import importlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


# =============================================================================
# PROJECT PATHS
# =============================================================================

ANALYSIS_DIR = Path(__file__).resolve().parent
ROOT = ANALYSIS_DIR.parent.parent

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

import sensitivity_common as sc  # noqa: E402


# =============================================================================
# USER SETTINGS
# =============================================================================

RUN_LNG_BASE_CASES = True

LNG_PROFILE_DAYS = 6.0
DUAL_FUEL_PROFILE_DAYS = 10.0

BASE_PARAMS = {
    "tank_ins_scale": 1.0,
    "initial_fill": 0.83,
    "T_amb": 293.15,
    "fuel_flow_scale": 1.0,
}

RESULTS_DIR = ROOT / "results" / "base_case"
TIMEPROFILE_DIR = RESULTS_DIR / "timeprofiles"
RAW_RESULTS_DIR = RESULTS_DIR / "raw_model_results"
SUMMARY_DIR = RESULTS_DIR / "summaries"


SYSTEM_LABELS = {
    "PBU": "FGSS 1 (PBU)",
    "Pump": "FGSS 2 (Pump)",
    "Compressor": "FGSS 3 (Compressor)",
}

SYSTEM_FILE_LABELS = {
    "PBU": "FGSS1_PBU",
    "Pump": "FGSS2_Pump",
    "Compressor": "FGSS3_Compressor",
}


# =============================================================================
# MANUSCRIPT TABLE DEFINITIONS
# =============================================================================

# These retain all quantities already reported in the manuscript and add the
# most decision-relevant pressure, BOG, and energy indicators.
MANUSCRIPT_METRICS = [
    (
        "planned_voyage_duration_days",
        "Planned voyage duration",
        "d",
        1.0,
    ),
    (
        "voyage_simulation_duration_days",
        "Analysed voyage duration",
        "d",
        1.0,
    ),
    (
        "simulation_duration_days",
        "Elapsed simulation duration",
        "d",
        1.0,
    ),
    (
        "voyage_time_to_min_heel_days",
        "Voyage time to minimum heel",
        "d",
        1.0,
    ),
    (
        "final_fill_pct",
        "Final tank filling",
        "%",
        1.0,
    ),
    (
        "max_pressure_bar",
        "Maximum tank pressure",
        "bar",
        1.0,
    ),
    (
        "pressure_at_termination_bar",
        "Final tank pressure",
        "bar",
        1.0,
    ),
    (
        "pressure_margin_to_mawp_bar",
        "Pressure margin to MAWP",
        "bar",
        1.0,
    ),
    (
        "total_bog_removed_kg",
        "Total BOG removed",
        "kg",
        1.0,
    ),
    (
        "bog_used_as_fuel_kg",
        "BOG used as fuel",
        "kg",
        1.0,
    ),
    (
        "bog_excess_kg",
        "Excess BOG requiring management",
        "kg",
        1.0,
    ),
    (
        "bog_removed_conventional_kg",
        "BOG removed during conventional-fuel operation",
        "kg",
        1.0,
    ),
    (
        "pbu_thermal_energy_kWh",
        "PBU thermal energy",
        "kWh",
        1.0,
    ),
    (
        "startup_pbu_thermal_energy_kWh",
        "Startup PBU thermal energy",
        "kWh",
        1.0,
    ),
    (
        "pump_energy_kWh",
        "Pump mechanical energy",
        "kWh",
        1.0,
    ),
    (
        "compressor_energy_kWh",
        "Compressor mechanical energy",
        "kWh",
        1.0,
    ),
    (
        "mechanical_energy_kWh",
        "Total mechanical energy",
        "kWh",
        1.0,
    ),
    (
        "specific_mechanical_energy_kWh_per_t",
        "Specific mechanical energy",
        "kWh t-1",
        1.0,
    ),
    (
        "requested_fuel_kg",
        "Requested LNG fuel",
        "kg",
        1.0,
    ),
    (
        "fuel_supplied_kg",
        "Supplied LNG fuel",
        "kg",
        1.0,
    ),
    (
        "lng_supply_shortfall_kg",
        "LNG-supply shortfall",
        "kg",
        1.0,
    ),
]


# =============================================================================
# ENGINE PROFILES
# =============================================================================

def build_lng_engine(
    fuel_flow_scale: float = 1.0,
    max_days: float = LNG_PROFILE_DAYS,
) -> sc.Engine:
    """
    Generate the repeated LNG-only operating profile.

    One 48 h cycle consists of:
        23 h port
         1 h manoeuvring
        23 h sailing
         1 h manoeuvring

    Unlike the dual-fuel profile, LNG demand is never intentionally set to zero.
    """
    hour = 3600.0
    day = 24.0 * hour

    times = np.arange(
        0.0,
        max_days * day + sc.DT,
        sc.DT,
    )

    cycle_phase = np.mod(times, 48.0 * hour)

    base_demand = np.select(
        [
            cycle_phase < 23.0 * hour,
            cycle_phase < 24.0 * hour,
            cycle_phase < 47.0 * hour,
        ],
        [
            sc.rezim.port_winter,
            sc.rezim.manouvering_winter,
            sc.rezim.service_winter,
        ],
        default=sc.rezim.manouvering_winter,
    )

    demands = (
        base_demand
        * 0.59
        * float(fuel_flow_scale)
    )

    pressures = np.full_like(
        times,
        600000.0,
        dtype=float,
    )
    heating_medium_temperatures = np.full_like(
        times,
        273.15 + 30.0,
        dtype=float,
    )

    return sc.Engine(
        demands,
        pressures,
        heating_medium_temperatures,
        times,
    )


# =============================================================================
# SYSTEM CONSTRUCTION
# =============================================================================

def build_system(
    system_key: str,
    params: dict[str, float],
):
    """
    Build one FGSS using the same baseline-parameter handling as the
    sensitivity analysis.
    """
    if system_key not in sc.SYSTEM_SPECS:
        raise KeyError(f"Unknown system key: {system_key}")

    spec = sc.SYSTEM_SPECS[system_key]
    lambda0, delt = sc._lock_baseline_params(spec)
    cfg = importlib.import_module(spec["config_module"])

    lambda_ins = lambda0 * params["tank_ins_scale"]

    cfg.k_liq = 1.0 / (
        delt / lambda_ins
        + 1.0 / 2000.0
        + 1.0 / 25.0
    )
    cfg.k_vap = 1.0 / (
        delt / lambda_ins
        + 1.0 / 200.0
        + 1.0 / 25.0
    )
    cfg.T_amb = params["T_amb"]
    cfg.T_env = params["T_amb"]
    cfg.liq_vol_ratio = params["initial_fill"]

    builder = getattr(
        cfg,
        spec["builder_name"],
    )

    return builder(
        list(sc.SLOSH_TIMES),
        spec["evap_flow"],
        spec["super_flow"],
        *spec["extra_args"],
    )


# =============================================================================
# OUTPUT HELPERS
# =============================================================================

def ensure_output_directories() -> None:
    for directory in (
        RESULTS_DIR,
        TIMEPROFILE_DIR,
        RAW_RESULTS_DIR,
        SUMMARY_DIR,
    ):
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )


def flatten_column_name(column: Any) -> str:
    if isinstance(column, tuple):
        parts = [
            str(part).strip()
            for part in column
            if str(part).strip()
        ]
        return "__".join(parts) if parts else "value"

    return str(column).strip() or "value"


def make_unique_names(names: list[str]) -> list[str]:
    counters: dict[str, int] = {}
    unique: list[str] = []

    for name in names:
        count = counters.get(name, 0)
        counters[name] = count + 1

        if count == 0:
            unique.append(name)
        else:
            unique.append(f"{name}__{count + 1}")

    return unique


def save_raw_model_results(
    system,
    output_path: Path,
) -> None:
    """
    Save every variable available in system.results.

    The model uses a two-level pandas column index. The names are flattened to
    component__variable so the CSV remains easy to inspect and reuse.
    """
    raw = system.results.copy()

    flattened = [
        flatten_column_name(column)
        for column in raw.columns
    ]
    raw.columns = make_unique_names(flattened)

    if "time" in raw.columns:
        raw = raw.rename(columns={"time": "time_s"})

    raw.to_csv(
        output_path,
        index=False,
    )


def add_identification_fields(
    metrics: dict[str, Any],
    profile: str,
    system_key: str,
    source: str,
) -> dict[str, Any]:
    result = {
        "profile": profile,
        "system_key": system_key,
        "configuration": SYSTEM_LABELS[system_key],
        "source": source,
    }
    result.update(metrics)
    return result


def run_one_case(
    system_key: str,
    profile_key: str,
    engine: sc.Engine,
) -> dict[str, Any]:
    print(
        f"\nRunning {SYSTEM_LABELS[system_key]} "
        f"for the {profile_key} profile..."
    )

    system = build_system(
        system_key,
        BASE_PARAMS,
    )
    system.calculate(engine)

    file_label = SYSTEM_FILE_LABELS[system_key]
    profile_label = profile_key.replace("-", "_")

    curated_path = (
        TIMEPROFILE_DIR
        / f"{file_label}_{profile_label}_timeprofile.csv"
    )
    raw_path = (
        RAW_RESULTS_DIR
        / f"{file_label}_{profile_label}_raw_results.csv"
    )
    summary_path = (
        SUMMARY_DIR
        / f"{file_label}_{profile_label}_summary.csv"
    )

    sc.save_time_profile(
        system,
        engine,
        system_key,
        curated_path,
    )
    save_raw_model_results(
        system,
        raw_path,
    )

    metrics = sc.extract_metrics(
        system,
        engine,
        system_key,
    )

    identified = add_identification_fields(
        metrics=metrics,
        profile=profile_key,
        system_key=system_key,
        source="simulation",
    )

    planned_voyage_duration_days = float(engine.times[-1]) / 86400.0
    identified["planned_voyage_duration_days"] = (
        planned_voyage_duration_days
    )
    actual_voyage_duration_days = float(
        identified["voyage_simulation_duration_days"]
    )
    identified["completed_planned_window"] = bool(
        actual_voyage_duration_days
        >= planned_voyage_duration_days
        - 1.5 * sc.DT / 86400.0
    )

    pd.DataFrame([identified]).to_csv(
        summary_path,
        index=False,
    )

    print(f"  Comprehensive time profile: {curated_path}")
    print(f"  Raw model results:          {raw_path}")
    print(f"  Summary:                    {summary_path}")

    return identified


# =============================================================================
# COMBINED MANUSCRIPT OUTPUTS
# =============================================================================

def save_combined_summaries(
    lng_rows: list[dict[str, Any]],
    dual_rows: list[dict[str, Any]],
) -> pd.DataFrame:
    lng_df = pd.DataFrame(lng_rows)
    dual_df = pd.DataFrame(dual_rows)

    lng_path = SUMMARY_DIR / "LNG_base_case_summary.csv"
    dual_path = SUMMARY_DIR / "dual_base_case_summary.csv"
    all_path = SUMMARY_DIR / "all_base_case_summary.csv"

    lng_df.to_csv(lng_path, index=False)
    dual_df.to_csv(dual_path, index=False)

    combined = pd.concat(
        [lng_df, dual_df],
        ignore_index=True,
        sort=False,
    )
    combined.to_csv(all_path, index=False)

    print(f"\nCombined LNG summary:       {lng_path}")
    print(f"Combined dual summary:      {dual_path}")
    print(f"Combined all-profile file:  {all_path}")

    return combined


def save_manuscript_table(
    combined: pd.DataFrame,
) -> Path:
    """
    Create a wide, manuscript-oriented CSV.

    The first two columns contain the indicator and unit. Remaining columns are:
        LNG — FGSS 1
        LNG — FGSS 2
        LNG — FGSS 3
        Dual-fuel — FGSS 1
        Dual-fuel — FGSS 2
        Dual-fuel — FGSS 3
    """
    profile_order = [
        "LNG",
        "dual-fuel",
    ]
    system_order = [
        "PBU",
        "Pump",
        "Compressor",
    ]

    output_rows: list[dict[str, Any]] = []

    for metric, label, unit, factor in MANUSCRIPT_METRICS:
        row: dict[str, Any] = {
            "performance_indicator": label,
            "unit": unit,
        }

        for profile in profile_order:
            for system_key in system_order:
                selection = combined.loc[
                    (combined["profile"] == profile)
                    & (
                        combined["system_key"]
                        == system_key
                    )
                ]

                column_name = (
                    f"{profile} — "
                    f"{SYSTEM_LABELS[system_key]}"
                )

                if (
                    len(selection) == 1
                    and metric in selection.columns
                ):
                    value = selection.iloc[0][metric]
                    if pd.notna(value):
                        row[column_name] = (
                            float(value)
                            * factor
                        )
                    else:
                        row[column_name] = np.nan
                else:
                    row[column_name] = np.nan

        output_rows.append(row)

    output = pd.DataFrame(output_rows)
    output_path = (
        SUMMARY_DIR
        / "manuscript_base_case_table.csv"
    )
    output.to_csv(
        output_path,
        index=False,
    )

    print(f"Manuscript-oriented table:  {output_path}")
    return output_path


def validate_summary_row(
    row: dict[str, Any],
) -> None:
    configuration = row["configuration"]
    profile = row["profile"]

    termination = str(
        row.get("termination_reason", "")
    )
    if termination not in {
        "minimum_heel",
        "MAWP",
        "profile_end",
    }:
        print(
            f"WARNING: {configuration}, {profile}: "
            f"unexpected termination reason {termination!r}."
        )

    shortfall = row.get(
        "lng_supply_shortfall_kg",
        np.nan,
    )
    if pd.notna(shortfall) and float(shortfall) > 1.0:
        print(
            f"WARNING: {configuration}, {profile}: "
            f"LNG-supply shortfall = {float(shortfall):.3f} kg."
        )

    balance_error = row.get(
        "bog_balance_error_kg",
        np.nan,
    )
    if (
        pd.notna(balance_error)
        and abs(float(balance_error))
        > sc.BOG_BALANCE_TOL_KG
    ):
        print(
            f"WARNING: {configuration}, {profile}: "
            f"BOG balance error = {float(balance_error):.3f} kg."
        )


def save_manifest(
    lng_rows: list[dict[str, Any]],
    dual_rows: list[dict[str, Any]],
) -> Path:
    manifest = {
        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "project_root": ".",
        "settings": {
            "run_lng_base_cases": RUN_LNG_BASE_CASES,
            "lng_profile_days": LNG_PROFILE_DAYS,
            "dual_fuel_profile_days": DUAL_FUEL_PROFILE_DAYS,
            "pbu_startup_excluded_from_voyage_clock": True,
            "baseline_parameters": BASE_PARAMS,
            "time_step_s": sc.DT,
            "mawp_bar": sc.MAWP_BAR,
            "minimum_heel_fraction":
                sc.HEEL_FRACTION,
        },
        "lng_sources": [
            row.get("source")
            for row in lng_rows
        ],
        "dual_sources": [
            row.get("source")
            for row in dual_rows
        ],
        "output_directory": "results/base_case",
    }

    output_path = (
        RESULTS_DIR
        / "base_case_run_manifest.json"
    )
    output_path.write_text(
        json.dumps(
            manifest,
            indent=2,
        ),
        encoding="utf-8",
    )
    return output_path


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    ensure_output_directories()

    lng_rows: list[dict[str, Any]] = []
    dual_rows: list[dict[str, Any]] = []

    if RUN_LNG_BASE_CASES:
        for system_key in (
            "PBU",
            "Pump",
            "Compressor",
        ):
            engine = build_lng_engine(
                fuel_flow_scale=1.0,
                max_days=LNG_PROFILE_DAYS,
            )
            result = run_one_case(
                system_key=system_key,
                profile_key="LNG",
                engine=engine,
            )
            validate_summary_row(result)
            lng_rows.append(result)
    else:
        raise RuntimeError(
            "RUN_LNG_BASE_CASES is False. "
            "No LNG base-case summaries would be available."
        )

    for system_key in (
        "PBU",
        "Pump",
        "Compressor",
    ):
        engine = sc.build_engine(
            fuel_flow_scale=1.0,
            max_days=DUAL_FUEL_PROFILE_DAYS,
        )
        result = run_one_case(
            system_key=system_key,
            profile_key="dual-fuel",
            engine=engine,
        )
        validate_summary_row(result)
        dual_rows.append(result)

    combined = save_combined_summaries(
        lng_rows=lng_rows,
        dual_rows=dual_rows,
    )
    save_manuscript_table(combined)
    manifest_path = save_manifest(
        lng_rows=lng_rows,
        dual_rows=dual_rows,
    )

    print(f"Run manifest:               {manifest_path}")
    print(f"\nAll outputs saved under:\n  {RESULTS_DIR.resolve()}")


if __name__ == "__main__":
    main()
