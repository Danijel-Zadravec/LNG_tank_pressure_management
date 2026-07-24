Final fixed-window base-case runner

Place run_base_case_analysis_fixed_windows.py beside sensitivity_common.py in:
  scripts/analysis/

The script repeats all six base cases:
  - LNG profile: 6 voyage days
  - dual-fuel profile: 10 voyage days

For FGSS 1, PBU startup is performed before the voyage clock. The PBU model
stops earlier if the tank reaches the 5% minimum heel.

Run:
  python scripts/analysis/run_base_case_analysis_fixed_windows.py

Outputs are written under:
  results/base_case/
