# LNG tank pressure management

Research code accompanying the manuscript *System-Level Analysis of LNG Tank
Pressure Management in Dual-Fuel Ships*. The model compares three low-pressure
fuel gas supply system (FGSS) configurations: pressure build-up (PBU), LNG
pump, and BOG compressor.

## Setup

Python 3.12 is recommended (tested with Python 3.12.0). From the repository
root, create an isolated environment and install the recorded dependencies:

```bash
python -m venv .venv
```

On Windows:

```powershell
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Linux or macOS:

```bash
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Reproduce the analyses

Run all commands from the repository root.

Base cases for the six-day LNG and ten-day dual-fuel profiles:

```bash
python scripts/analysis/run_base_case_analysis.py
```

One-at-a-time sensitivity analysis for the three FGSS configurations:

```bash
python scripts/analysis/run_all_sensitivity.py
```

Regenerate the manuscript tornado figures and appendix workbook after the
sensitivity analysis:

```bash
python scripts/analysis/plot_manuscript_sensitivity_tornado_shared_y.py
python scripts/analysis/create_appendix_sensitivity_tables_compact.py
```

The simulations use a fixed 10 s time step and are computationally intensive.
Detailed time profiles and raw model tables are regenerated locally and can
occupy several gigabytes; they are intentionally excluded from version
control. Compact CSV summaries, publication figures, and the appendix workbook
under `results/` are retained as reference outputs.

## Repository layout

- `configurations/`: parameters and factory functions for the three FGSS
  architectures.
- `src/`: tank, thermodynamic-property, heat-exchanger, piping, and machinery
  models.
- `scripts/analysis/`: base-case, sensitivity, figure, and table scripts.
- `scripts/validation/`: tank pressurisation validation and source data.
- `results/`: compact reference results used in the manuscript.
- `Systems_paper.docx`: accompanying manuscript.

Model and operating assumptions are defined in the configuration modules and
near the top of the two analysis drivers. The tracked base-case manifest records
the principal run settings. The validation comparison can be run with
`python scripts/validation/tank_validation.py`; it opens a Matplotlib window and
reports the mean absolute relative deviation.

## License

This project is distributed under the terms of the [MIT License](LICENSE).
