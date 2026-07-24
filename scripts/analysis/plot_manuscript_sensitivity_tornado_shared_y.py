#!/usr/bin/env python3
"""
Create the manuscript tornado figures from the FGSS sensitivity summaries.

The script uses a fixed parameter order, marks MAWP-limited cases, and creates
one 2x2 figure for each FGSS configuration.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# =============================================================================
# USER SETTINGS
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

# Edit this if your CSV files are elsewhere.
DATA_DIR = PROJECT_ROOT / "results" / "sensitivity"
OUTPUT_DIR = PROJECT_ROOT / "results" / "sensitivity" /"figures"

SAVE_PNG = True
SAVE_PDF = True
SAVE_SVG = False
PNG_DPI = 600
SHOW_FIGURES = False

# Slightly wider than the previous version to eliminate label collisions.
FIGSIZE = (9.2, 6.6)

# Manuscript-like styling.
plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "font.size": 11.0,
        "axes.titlesize": 13.0,
        "axes.labelsize": 12.0,
        "xtick.labelsize": 11.0,
        "ytick.labelsize": 11.0,
        "legend.fontsize": 11.0,
        "axes.linewidth": 0.9,
        "savefig.bbox": "tight",
    }
)

# Colors chosen to resemble the manuscript time-profile figure.
COLOR_LOW = "#2b6cb0"       # blue
COLOR_HIGH = "#c73e2d"      # red
COLOR_BASELINE = "#5a5a5a"  # neutral gray
GRID_COLOR = "#cfcfcf"


# =============================================================================
# DEFINITIONS
# =============================================================================

PARAMETER_ORDER = [
    "tank_ins_scale",
    "initial_fill",
    "T_amb",
    "fuel_flow_scale",
]

PARAMETER_LABELS = {
    "tank_ins_scale": r"Insulation, $\lambda/\lambda_0$",
    "initial_fill": "Initial tank filling",
    "T_amb": "Ambient temperature",
    "fuel_flow_scale": "Fuel-demand scaling",
}

PANEL_LABELS = ["(a)", "(b)", "(c)", "(d)"]


@dataclass(frozen=True)
class MetricSpec:
    column: str
    title: str
    axis_label: str
    mawp_fallback_column: str | None = None


@dataclass(frozen=True)
class SystemSpec:
    key: str
    name: str
    csv_candidates: tuple[str, ...]
    output_stem: str
    metrics: tuple[MetricSpec, MetricSpec, MetricSpec, MetricSpec]


SYSTEMS = (
    SystemSpec(
        key="PBU",
        name="FGSS 1 (PBU)",
        csv_candidates=(
            "PBU_dual_sensitivity_results.csv",
            "PBU_dual_sensitivity_results(3).csv",
            "PBU_dual_sensitivity_results(2).csv",
            "PBU_dual_sensitivity_results(1).csv",
        ),
        output_stem="FGSS1_PBU_tornado_final",
        metrics=(
            MetricSpec(
                column="voyage_time_to_min_heel_days",
                title="Voyage time to minimum heel",
                axis_label="Voyage time (d)",
                mawp_fallback_column="time_to_mawp_days",
            ),
            MetricSpec(
                column="pressure_margin_to_mawp_bar",
                title="Pressure margin to MAWP",
                axis_label="Pressure margin (bar)",
            ),
            MetricSpec(
                column="bog_excess_kg",
                title="Excess BOG requiring management",
                axis_label="Excess BOG (t)",
            ),
            MetricSpec(
                column="pbu_thermal_energy_kWh",
                title="PBU thermal energy",
                axis_label="PBU thermal energy (kWh)",
            ),
        ),
    ),
    SystemSpec(
        key="Pump",
        name="FGSS 2 (Pump)",
        csv_candidates=(
            "Pump_dual_sensitivity_results.csv",
            "Pump_dual_sensitivity_results(3).csv",
            "Pump_dual_sensitivity_results(2).csv",
            "Pump_dual_sensitivity_results(1).csv",
        ),
        output_stem="FGSS2_Pump_tornado_final",
        metrics=(
            MetricSpec(
                column="voyage_time_to_min_heel_days",
                title="Voyage time to minimum heel",
                axis_label="Voyage time (d)",
                mawp_fallback_column="time_to_mawp_days",
            ),
            MetricSpec(
                column="pressure_margin_to_mawp_bar",
                title="Pressure margin to MAWP",
                axis_label="Pressure margin (bar)",
            ),
            MetricSpec(
                column="bog_excess_kg",
                title="Excess BOG requiring management",
                axis_label="Excess BOG (t)",
            ),
            MetricSpec(
                column="pump_energy_kWh",
                title="Pump mechanical energy",
                axis_label="Pump energy (kWh)",
            ),
        ),
    ),
    SystemSpec(
        key="Compressor",
        name="FGSS 3 (Compressor)",
        csv_candidates=(
            "Compressor_dual_sensitivity_results.csv",
            "Compressor_dual_sensitivity_results(2).csv",
            "Compressor_dual_sensitivity_results(1).csv",
        ),
        output_stem="FGSS3_Compressor_tornado_final",
        metrics=(
            MetricSpec(
                column="voyage_time_to_min_heel_days",
                title="Voyage time to minimum heel",
                axis_label="Voyage time (d)",
                mawp_fallback_column="time_to_mawp_days",
            ),
            MetricSpec(
                column="pressure_margin_to_mawp_bar",
                title="Pressure margin to MAWP",
                axis_label="Pressure margin (bar)",
            ),
            MetricSpec(
                column="bog_used_as_fuel_kg",
                title="BOG used as fuel",
                axis_label="BOG used as fuel (t)",
            ),
            MetricSpec(
                column="mechanical_energy_kWh",
                title="Total mechanical energy",
                axis_label="Mechanical energy (kWh)",
            ),
        ),
    ),
)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def first_existing_valid_file(folder: Path, candidates: Iterable[str], required_columns: set[str]) -> Path:
    existing = []
    for name in candidates:
        path = folder / name
        if path.is_file():
            existing.append(path)
            try:
                preview = pd.read_csv(path, nrows=1)
            except Exception:
                continue
            if required_columns.issubset(set(preview.columns)):
                return path
    if existing:
        raise ValueError(
            "CSV files were found, but none contains the required columns "
            f"{sorted(required_columns)}.\nChecked:\n  " + "\n  ".join(str(p) for p in existing)
        )
    raise FileNotFoundError(
        "None of the expected CSV files were found in\n"
        f"  {folder}\nCandidates:\n  " + "\n  ".join(candidates)
    )


def read_results(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"parameter", "level", "termination_reason"}
    missing = needed.difference(df.columns)
    if missing:
        raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")
    return df


def get_baseline_row(df: pd.DataFrame) -> pd.Series:
    rows = df.loc[df["parameter"] == "baseline"]
    if len(rows) != 1:
        raise ValueError("Each CSV must contain exactly one baseline row.")
    return rows.iloc[0]


def get_case_row(df: pd.DataFrame, parameter: str, level: str) -> pd.Series:
    rows = df.loc[(df["parameter"] == parameter) & (df["level"] == level)]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one row for parameter={parameter!r}, level={level!r}; got {len(rows)}"
        )
    return rows.iloc[0]


def convert_value(metric: MetricSpec, value: float) -> float:
    if not np.isfinite(value):
        return np.nan
    if metric.column in {"bog_excess_kg", "bog_used_as_fuel_kg"}:
        return value / 1000.0
    return value


def get_metric_value(row: pd.Series, metric: MetricSpec) -> tuple[float, bool, bool]:
    """
    Returns:
        plotted_value, is_mawp_limited, used_fallback_value
    """
    term = str(row.get("termination_reason", "")).strip().upper()
    mawp = term == "MAWP"

    raw = row.get(metric.column, np.nan)
    value = float(raw) if pd.notna(raw) else np.nan
    used_fallback = False

    if (not np.isfinite(value)) and mawp and metric.mawp_fallback_column:
        fb = row.get(metric.mawp_fallback_column, np.nan)
        if pd.notna(fb):
            value = float(fb)
            used_fallback = True

    return convert_value(metric, value), mawp, used_fallback


def prepare_panel_data(df: pd.DataFrame, metric: MetricSpec) -> pd.DataFrame:
    baseline_row = get_baseline_row(df)
    baseline = convert_value(metric, float(baseline_row[metric.column]))

    records = []
    for param in PARAMETER_ORDER:
        low_row = get_case_row(df, param, "low")
        high_row = get_case_row(df, param, "high")

        low_val, low_mawp, low_fb = get_metric_value(low_row, metric)
        high_val, high_mawp, high_fb = get_metric_value(high_row, metric)

        records.append(
            {
                "parameter": param,
                "label": PARAMETER_LABELS[param],
                "baseline": baseline,
                "low_value": low_val,
                "high_value": high_val,
                "low_mawp": low_mawp,
                "high_mawp": high_mawp,
                "low_fallback": low_fb,
                "high_fallback": high_fb,
            }
        )

    return pd.DataFrame.from_records(records)


def set_xlim_with_padding(ax: plt.Axes, panel_data: pd.DataFrame) -> None:
    values = [float(panel_data["baseline"].iloc[0])]
    for col in ("low_value", "high_value"):
        values.extend(panel_data.loc[np.isfinite(panel_data[col]), col].astype(float).tolist())

    vmin = min(values)
    vmax = max(values)
    span = vmax - vmin
    if np.isclose(span, 0.0):
        span = max(abs(vmax), 1.0) * 0.15
    pad = 0.10 * span
    ax.set_xlim(vmin - pad, vmax + pad)


# =============================================================================
# PLOTTING
# =============================================================================

def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontweight="bold",
        fontsize=11.0,
        clip_on=False,
    )


def plot_panel(
    ax: plt.Axes,
    panel_data: pd.DataFrame,
    metric: MetricSpec,
    panel_label: str,
    show_y_labels: bool,
) -> None:
    # Keep the same parameter order everywhere; first parameter at the top.
    y = np.arange(len(panel_data), dtype=float)
    baseline = panel_data["baseline"].to_numpy(dtype=float)
    low_value = panel_data["low_value"].to_numpy(dtype=float)
    high_value = panel_data["high_value"].to_numpy(dtype=float)

    low_left = np.minimum(low_value, baseline)
    high_left = np.minimum(high_value, baseline)
    low_width = np.abs(low_value - baseline)
    high_width = np.abs(high_value - baseline)

    low_bars = ax.barh(
        y - 0.17,
        np.nan_to_num(low_width, nan=0.0),
        left=np.nan_to_num(low_left, nan=baseline),
        height=0.30,
        color=COLOR_LOW,
        edgecolor=COLOR_LOW,
        label="Low parameter level",
        zorder=3,
    )
    high_bars = ax.barh(
        y + 0.17,
        np.nan_to_num(high_width, nan=0.0),
        left=np.nan_to_num(high_left, nan=baseline),
        height=0.30,
        color=COLOR_HIGH,
        edgecolor=COLOR_HIGH,
        label="High parameter level",
        zorder=3,
    )

    for patch, is_mawp in zip(low_bars.patches, panel_data["low_mawp"]):
        if bool(is_mawp):
            patch.set_hatch("///")
            patch.set_linewidth(0.9)
            patch.set_edgecolor("black")
    for patch, is_mawp in zip(high_bars.patches, panel_data["high_mawp"]):
        if bool(is_mawp):
            patch.set_hatch("///")
            patch.set_linewidth(0.9)
            patch.set_edgecolor("black")

    ax.axvline(
        float(panel_data["baseline"].iloc[0]),
        color=COLOR_BASELINE,
        linewidth=1.3,
        zorder=4,
    )

    ax.set_yticks(y)
    if show_y_labels:
        ax.set_yticklabels(panel_data["label"])
        ax.tick_params(axis="y", which="both", left=True, labelleft=True, pad=5)
    else:
        # Panels (b) and (d) share the same ordered parameters as the panels
        # immediately to their left, so repeating the labels is unnecessary.
        ax.tick_params(axis="y", which="both", left=False, labelleft=False)

    ax.invert_yaxis()
    ax.set_title(metric.title, pad=8)
    ax.set_xlabel(metric.axis_label)
    ax.grid(axis="x", color=GRID_COLOR, linewidth=0.7, alpha=0.85)
    ax.set_axisbelow(True)

    set_xlim_with_padding(ax, panel_data)
    add_panel_label(ax, panel_label)

    # Concise annotation for MAWP fallback in the voyage-time panel.
    if metric.mawp_fallback_column is not None:
        x0, x1 = ax.get_xlim()
        dx = 0.012 * (x1 - x0)
        for idx, row in panel_data.iterrows():
            if row["low_fallback"] and np.isfinite(row["low_value"]):
                x = float(row["low_value"])
                ax.text(x - dx, y[idx] - 0.17, "MAWP", ha="right", va="center", fontsize=8.0)
            if row["high_fallback"] and np.isfinite(row["high_value"]):
                x = float(row["high_value"])
                ax.text(x - dx, y[idx] + 0.17, "MAWP", ha="right", va="center", fontsize=8.0)


def save_figure(fig: plt.Figure, out_stem: Path) -> None:
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    if SAVE_PNG:
        fig.savefig(out_stem.with_suffix(".png"), dpi=PNG_DPI)
    if SAVE_PDF:
        fig.savefig(out_stem.with_suffix(".pdf"))
    if SAVE_SVG:
        fig.savefig(out_stem.with_suffix(".svg"))


def make_figure(system: SystemSpec, df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=FIGSIZE, sharey=True)
    axes_flat = axes.ravel()

    for i, (ax, metric) in enumerate(zip(axes_flat, system.metrics)):
        pdata = prepare_panel_data(df, metric)
        plot_panel(
            ax=ax,
            panel_data=pdata,
            metric=metric,
            panel_label=PANEL_LABELS[i],
            show_y_labels=(i % 2 == 0),
        )

    fig.suptitle(
        f"One-at-a-time sensitivity analysis — {system.name}",
        y=0.985,
        fontsize=14.0,
    )

    legend_handles = [
        Patch(facecolor=COLOR_LOW, edgecolor=COLOR_LOW, label="Low parameter level"),
        Patch(facecolor=COLOR_HIGH, edgecolor=COLOR_HIGH, label="High parameter level"),
        Line2D([0], [0], color=COLOR_BASELINE, linewidth=1.3, label="Baseline"),
    ]

    if df["termination_reason"].astype(str).str.strip().str.upper().eq("MAWP").any():
        legend_handles.append(
            Patch(facecolor="white", edgecolor="black", hatch="///", label="MAWP-limited case")
        )

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=len(legend_handles),
        frameon=False,
        bbox_to_anchor=(0.5, 0.01),
        columnspacing=1.5,
        handlelength=2.3,
    )

    # Spacing tuned to prevent any text/image overlap.
    fig.subplots_adjust(
        left=0.22,
        right=0.97,
        top=0.88,
        bottom=0.14,
        wspace=0.20,
        hspace=0.52,
    )

    save_figure(fig, OUTPUT_DIR / system.output_stem)

    if SHOW_FIGURES:
        plt.show()
    else:
        plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for system in SYSTEMS:
        csv_path = first_existing_valid_file(
            DATA_DIR, system.csv_candidates, {"parameter", "level", "termination_reason"}
        )
        df = read_results(csv_path)
        make_figure(system, df)
        print(f"Created {system.name} from {csv_path.name}")

    print(f"\nSaved figures to:\n  {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
