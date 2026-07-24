#!/usr/bin/env python3
"""
Publication-ready tornado plots for the FGSS sensitivity analysis.

The script reads the three one-at-a-time sensitivity-result CSV files and
creates one four-panel figure for each FGSS configuration:

    (a) voyage time to minimum heel,
    (b) pressure margin to MAWP,
    (c) BOG-management indicator,
    (d) system-specific equipment energy.

The figures are intended for the manuscript, while the complete numerical
results can be reported in the Supplementary Materials.

VS Code use
-----------
1. Put this script in the project folder.
2. Set DATA_DIR below to the folder containing the three CSV files.
3. Run the script using the "Run Python File" button in VS Code.

Expected CSV file names
-----------------------
PBU_dual_sensitivity_results.csv
Pump_dual_sensitivity_results.csv
Compressor_dual_sensitivity_results.csv

Outputs
-------
figures/sensitivity/FGSS1_PBU_sensitivity_tornado.pdf
figures/sensitivity/FGSS1_PBU_sensitivity_tornado.png
figures/sensitivity/FGSS2_Pump_sensitivity_tornado.pdf
figures/sensitivity/FGSS2_Pump_sensitivity_tornado.png
figures/sensitivity/FGSS3_Compressor_sensitivity_tornado.pdf
figures/sensitivity/FGSS3_Compressor_sensitivity_tornado.png
figures/sensitivity/tornado_plot_data.csv
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


# =============================================================================
# USER SETTINGS
# =============================================================================

# Folder containing this script.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
# Edit this line if the CSV files are stored elsewhere.
# Example:
# DATA_DIR = Path(r"C:\Users\Danijel\Documents\FGSS\results")
DATA_DIR = PROJECT_ROOT / "results" / "sensitivity"

PBU_CSV = DATA_DIR / "PBU_dual_sensitivity_results.csv"
PUMP_CSV = DATA_DIR / "Pump_dual_sensitivity_results.csv"
COMPRESSOR_CSV = DATA_DIR / "Compressor_dual_sensitivity_results.csv"

OUTPUT_DIR = PROJECT_ROOT / "results" / "sensitivity" /"figures"

SAVE_PDF = True
SAVE_PNG = True
SAVE_SVG = False
PNG_DPI = 600
SHOW_FIGURES = False

# Full-width manuscript figure. Adjust only if required by the journal template.
FIGURE_SIZE_INCHES = (7.2, 5.8)

# Use a generic serif family so the script remains portable.
plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
        "axes.linewidth": 0.8,
        "savefig.bbox": "tight",
    }
)


# =============================================================================
# MANUSCRIPT DEFINITIONS
# =============================================================================

PARAMETER_ORDER = [
    "tank_ins_scale",
    "initial_fill",
    "T_amb",
    "fuel_flow_scale",
]

PARAMETER_LABELS = {
    "tank_ins_scale": r"Insulation conductivity, $\lambda/\lambda_0$",
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
    decimals: int
    mawp_fallback_column: str | None = None


@dataclass(frozen=True)
class SystemSpec:
    key: str
    name: str
    csv_path: Path
    output_stem: str
    metrics: tuple[MetricSpec, MetricSpec, MetricSpec, MetricSpec]


SYSTEMS = (
    SystemSpec(
        key="PBU",
        name="FGSS 1 (PBU)",
        csv_path=PBU_CSV,
        output_stem="FGSS1_PBU_sensitivity_tornado",
        metrics=(
            MetricSpec(
                column="voyage_time_to_min_heel_days",
                title="Voyage time to minimum heel",
                axis_label="Voyage time (d)",
                decimals=2,
                mawp_fallback_column="time_to_mawp_days",
            ),
            MetricSpec(
                column="pressure_margin_to_mawp_bar",
                title="Pressure margin to MAWP",
                axis_label="Pressure margin (bar)",
                decimals=2,
            ),
            MetricSpec(
                column="bog_excess_kg",
                title="Excess BOG requiring management",
                axis_label="Excess BOG (t)",
                decimals=2,
            ),
            MetricSpec(
                column="pbu_thermal_energy_kWh",
                title="PBU thermal energy",
                axis_label="PBU thermal energy (kWh)",
                decimals=0,
            ),
        ),
    ),
    SystemSpec(
        key="Pump",
        name="FGSS 2 (Pump)",
        csv_path=PUMP_CSV,
        output_stem="FGSS2_Pump_sensitivity_tornado",
        metrics=(
            MetricSpec(
                column="voyage_time_to_min_heel_days",
                title="Voyage time to minimum heel",
                axis_label="Voyage time (d)",
                decimals=2,
                mawp_fallback_column="time_to_mawp_days",
            ),
            MetricSpec(
                column="pressure_margin_to_mawp_bar",
                title="Pressure margin to MAWP",
                axis_label="Pressure margin (bar)",
                decimals=2,
            ),
            MetricSpec(
                column="bog_excess_kg",
                title="Excess BOG requiring management",
                axis_label="Excess BOG (t)",
                decimals=2,
            ),
            MetricSpec(
                column="pump_energy_kWh",
                title="Pump mechanical energy",
                axis_label="Pump energy (kWh)",
                decimals=0,
            ),
        ),
    ),
    SystemSpec(
        key="Compressor",
        name="FGSS 3 (Compressor)",
        csv_path=COMPRESSOR_CSV,
        output_stem="FGSS3_Compressor_sensitivity_tornado",
        metrics=(
            MetricSpec(
                column="voyage_time_to_min_heel_days",
                title="Voyage time to minimum heel",
                axis_label="Voyage time (d)",
                decimals=2,
                mawp_fallback_column="time_to_mawp_days",
            ),
            MetricSpec(
                column="pressure_margin_to_mawp_bar",
                title="Pressure margin to MAWP",
                axis_label="Pressure margin (bar)",
                decimals=2,
            ),
            MetricSpec(
                column="bog_used_as_fuel_kg",
                title="BOG used as fuel",
                axis_label="BOG used as fuel (t)",
                decimals=2,
            ),
            MetricSpec(
                column="mechanical_energy_kWh",
                title="Total mechanical energy",
                axis_label="Mechanical energy (kWh)",
                decimals=0,
            ),
        ),
    ),
)


# =============================================================================
# DATA HANDLING
# =============================================================================

def read_results(path: Path) -> pd.DataFrame:
    """Read and validate one sensitivity-result CSV file."""
    if not path.is_file():
        raise FileNotFoundError(
            f"\nCSV file not found:\n  {path}\n\n"
            "Edit DATA_DIR in the USER SETTINGS section of the script."
        )

    df = pd.read_csv(path)

    required_columns = {
        "parameter",
        "level",
        "value",
        "termination_reason",
    }
    missing = required_columns.difference(df.columns)
    if missing:
        raise ValueError(
            f"{path.name} is missing required columns: {sorted(missing)}"
        )

    return df


def get_baseline_row(df: pd.DataFrame) -> pd.Series:
    """Return the unique explicit baseline row."""
    baseline_rows = df.loc[df["parameter"] == "baseline"]

    if len(baseline_rows) != 1:
        raise ValueError(
            "Each CSV must contain exactly one row with parameter='baseline'."
        )

    return baseline_rows.iloc[0]


def get_case_row(
    df: pd.DataFrame,
    parameter: str,
    level: str,
) -> pd.Series:
    """Return one low- or high-level sensitivity case."""
    rows = df.loc[
        (df["parameter"] == parameter)
        & (df["level"] == level)
    ]

    if len(rows) != 1:
        raise ValueError(
            f"Expected one row for parameter='{parameter}', level='{level}', "
            f"but found {len(rows)}."
        )

    return rows.iloc[0]


def convert_for_plot(metric: MetricSpec, value: float) -> float:
    """Convert mass indicators from kg to t; leave other quantities unchanged."""
    if not np.isfinite(value):
        return np.nan

    if metric.column in {"bog_excess_kg", "bog_used_as_fuel_kg"}:
        return value / 1000.0

    return value


def obtain_plot_value(
    row: pd.Series,
    metric: MetricSpec,
) -> tuple[float, bool, str | None]:
    """
    Return the plotted value, MAWP flag, and optional annotation.

    For voyage time, a case reaching MAWP before the minimum heel has no valid
    time-to-heel value. The time to MAWP is plotted instead and hatched.
    """
    termination = str(row["termination_reason"]).strip().upper()
    mawp_limited = termination == "MAWP"

    raw_value = row.get(metric.column, np.nan)
    value = float(raw_value) if pd.notna(raw_value) else np.nan
    annotation = None

    if (
        not np.isfinite(value)
        and mawp_limited
        and metric.mawp_fallback_column is not None
    ):
        fallback = row.get(metric.mawp_fallback_column, np.nan)
        if pd.notna(fallback):
            value = float(fallback)
            annotation = "MAWP"

    return convert_for_plot(metric, value), mawp_limited, annotation


def prepare_metric_data(
    df: pd.DataFrame,
    metric: MetricSpec,
) -> pd.DataFrame:
    """Prepare and rank the low/high effects for one tornado panel."""
    baseline_row = get_baseline_row(df)
    baseline_raw = float(baseline_row[metric.column])
    baseline = convert_for_plot(metric, baseline_raw)

    if not np.isfinite(baseline):
        raise ValueError(
            f"Baseline value for '{metric.column}' is not finite."
        )

    records: list[dict[str, object]] = []

    for parameter in PARAMETER_ORDER:
        low_row = get_case_row(df, parameter, "low")
        high_row = get_case_row(df, parameter, "high")

        low_value, low_mawp, low_note = obtain_plot_value(low_row, metric)
        high_value, high_mawp, high_note = obtain_plot_value(high_row, metric)

        finite_deviations = [
            abs(value - baseline)
            for value in (low_value, high_value)
            if np.isfinite(value)
        ]
        impact = max(finite_deviations, default=0.0)

        records.append(
            {
                "parameter": parameter,
                "label": PARAMETER_LABELS[parameter],
                "baseline": baseline,
                "low_value": low_value,
                "high_value": high_value,
                "low_mawp": low_mawp,
                "high_mawp": high_mawp,
                "low_note": low_note,
                "high_note": high_note,
                "impact": impact,
            }
        )

    # Ascending order is used because the last item appears at the top of barh.
    return (
        pd.DataFrame.from_records(records)
        .sort_values("impact", ascending=True)
        .reset_index(drop=True)
    )


# =============================================================================
# PLOTTING
# =============================================================================

def add_panel_label(ax: plt.Axes, label: str) -> None:
    """Place a manuscript panel label above the upper-left corner."""
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontweight="bold",
        clip_on=False,
    )


def apply_mawp_hatch(
    bar_container,
    mawp_flags: pd.Series,
) -> None:
    """Hatch bars corresponding to MAWP-limited simulations."""
    for patch, is_mawp in zip(bar_container.patches, mawp_flags):
        if bool(is_mawp):
            patch.set_hatch("///")
            patch.set_linewidth(0.8)


def annotate_mawp_cases(
    ax: plt.Axes,
    panel_data: pd.DataFrame,
    y_positions: np.ndarray,
) -> None:
    """Add a concise MAWP label to fallback voyage-time bars."""
    x_min, x_max = ax.get_xlim()
    x_span = x_max - x_min
    offset = 0.012 * x_span

    for index, row in panel_data.iterrows():
        if row["low_note"] == "MAWP" and np.isfinite(row["low_value"]):
            x = float(row["low_value"])
            ha = "left" if x >= row["baseline"] else "right"
            dx = offset if ha == "left" else -offset
            ax.text(
                x + dx,
                y_positions[index] - 0.17,
                "MAWP",
                ha=ha,
                va="center",
                fontsize=7.2,
            )

        if row["high_note"] == "MAWP" and np.isfinite(row["high_value"]):
            x = float(row["high_value"])
            ha = "left" if x >= row["baseline"] else "right"
            dx = offset if ha == "left" else -offset
            ax.text(
                x + dx,
                y_positions[index] + 0.17,
                "MAWP",
                ha=ha,
                va="center",
                fontsize=7.2,
            )


def set_reasonable_xlim(
    ax: plt.Axes,
    panel_data: pd.DataFrame,
) -> None:
    """Set symmetric visual padding around all finite plotted values."""
    values = [float(panel_data["baseline"].iloc[0])]

    for column in ("low_value", "high_value"):
        values.extend(
            panel_data.loc[
                np.isfinite(panel_data[column]),
                column,
            ].astype(float)
        )

    value_min = min(values)
    value_max = max(values)
    span = value_max - value_min

    if np.isclose(span, 0.0):
        span = max(abs(value_max), 1.0) * 0.2

    padding = 0.10 * span
    ax.set_xlim(value_min - padding, value_max + padding)


def plot_tornado_panel(
    ax: plt.Axes,
    panel_data: pd.DataFrame,
    metric: MetricSpec,
    panel_label: str,
):
    """Draw one publication-style tornado panel."""
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
        label="Low parameter level",
    )
    high_bars = ax.barh(
        y + 0.17,
        np.nan_to_num(high_width, nan=0.0),
        left=np.nan_to_num(high_left, nan=baseline),
        height=0.30,
        label="High parameter level",
    )

    apply_mawp_hatch(low_bars, panel_data["low_mawp"])
    apply_mawp_hatch(high_bars, panel_data["high_mawp"])

    ax.axvline(
        float(panel_data["baseline"].iloc[0]),
        linewidth=1.0,
        label="Baseline",
    )

    ax.set_yticks(y)
    ax.set_yticklabels(panel_data["label"])
    ax.set_title(metric.title, pad=6)
    ax.set_xlabel(metric.axis_label)
    ax.grid(axis="x", alpha=0.28, linewidth=0.6)
    ax.set_axisbelow(True)

    set_reasonable_xlim(ax, panel_data)
    annotate_mawp_cases(ax, panel_data, y)
    add_panel_label(ax, panel_label)

    return low_bars, high_bars


def save_figure(
    fig: plt.Figure,
    output_stem: str,
) -> None:
    """Save the figure in manuscript and preview formats."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    base_path = OUTPUT_DIR / output_stem

    if SAVE_PDF:
        fig.savefig(base_path.with_suffix(".pdf"))
    if SAVE_PNG:
        fig.savefig(base_path.with_suffix(".png"), dpi=PNG_DPI)
    if SAVE_SVG:
        fig.savefig(base_path.with_suffix(".svg"))


def plot_system_figure(
    system: SystemSpec,
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Create one four-panel sensitivity figure for an FGSS configuration."""
    fig, axes = plt.subplots(
        2,
        2,
        figsize=FIGURE_SIZE_INCHES,
    )
    axes_flat = axes.ravel()

    export_frames: list[pd.DataFrame] = []
    legend_handles = None

    for index, (ax, metric) in enumerate(zip(axes_flat, system.metrics)):
        panel_data = prepare_metric_data(df, metric)

        low_bars, high_bars = plot_tornado_panel(
            ax=ax,
            panel_data=panel_data,
            metric=metric,
            panel_label=PANEL_LABELS[index],
        )

        if legend_handles is None:
            legend_handles = [
                low_bars.patches[0],
                high_bars.patches[0],
            ]

        export_data = panel_data.copy()
        export_data.insert(0, "system", system.name)
        export_data.insert(1, "metric", metric.column)
        export_frames.append(export_data)

    fig.suptitle(
        f"One-at-a-time sensitivity analysis — {system.name}",
        y=0.995,
        fontsize=10.5,
    )

    baseline_handle = plt.Line2D(
        [0],
        [0],
        linewidth=1.0,
        label="Baseline",
    )

    figure_handles = [
        legend_handles[0],
        legend_handles[1],
        baseline_handle,
    ]
    figure_labels = [
        "Low parameter level",
        "High parameter level",
        "Baseline",
    ]

    has_mawp_case = (
        df["termination_reason"]
        .astype(str)
        .str.strip()
        .str.upper()
        .eq("MAWP")
        .any()
    )
    if has_mawp_case:
        figure_handles.append(
            Patch(
                fill=False,
                hatch="///",
                label="MAWP-limited case",
            )
        )
        figure_labels.append("MAWP-limited case")

    fig.legend(
        handles=figure_handles,
        labels=figure_labels,
        loc="lower center",
        ncol=len(figure_handles),
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
    )

    fig.subplots_adjust(
        left=0.17,
        right=0.98,
        top=0.90,
        bottom=0.12,
        wspace=0.50,
        hspace=0.48,
    )

    save_figure(fig, system.output_stem)

    if SHOW_FIGURES:
        plt.show()
    else:
        plt.close(fig)

    return pd.concat(export_frames, ignore_index=True)


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    """Read all CSV files and generate the final manuscript figures."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_plot_data: list[pd.DataFrame] = []

    for system in SYSTEMS:
        df = read_results(system.csv_path)
        plot_data = plot_system_figure(system, df)
        all_plot_data.append(plot_data)

        print(
            f"Created {system.name}: "
            f"{(OUTPUT_DIR / system.output_stem).resolve()}"
        )

    combined = pd.concat(all_plot_data, ignore_index=True)
    combined.to_csv(
        OUTPUT_DIR / "tornado_plot_data.csv",
        index=False,
    )

    print(f"\nAll outputs saved to:\n{OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
