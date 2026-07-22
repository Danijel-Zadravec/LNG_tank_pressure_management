# -*- coding: utf-8 -*-
"""
One-at-a-time (OAT) sensitivity analysis for the PBU dual-fuel FGSS.

Base case: scripts/analysis/PBU_dual.py

Swept parameters (one at a time, others held at base value):
    tank_ins_scale  : insulation thermal conductivity scale, \u03bb/\u03bb0  {0.5, 1.0, 1.5}
    initial_fill    : initial liquid volume ratio         {0.75, 0.83, 0.90}
    T_amb           : ambient temperature, K              {278.15, 293.15, 308.15}
    fuel_flow_scale : fuel flow scaling                   {0.75, 1.00, 1.25}

Reports (per case): termination reason (minimum_heel / MAWP / profile_end),
time to minimum heel, maximum tank pressure, pressure at termination, margin
to MAWP, total BOG removed, excess BOG requiring management, BOG used as
fuel, and mechanical (pump + compressor) energy -- always 0 for PBU, which
has no rotating machinery and instead consumes thermal energy from the
glycol heating medium. Saves a long-form CSV and a publication-style tornado
diagram (one panel per metric).

NOTE: this script is prepared but intentionally not executed here.
"""
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

import sensitivity_common as sc

SYSTEM_KEY = 'PBU'

if __name__ == '__main__':
    df, base_metrics = sc.run_oat_sensitivity(SYSTEM_KEY)
    out_png = sc.RESULTS_DIR / f'{SYSTEM_KEY}_dual_tornado.png'
    sc.plot_tornado(
        df, base_metrics, sc.SYSTEM_SPECS[SYSTEM_KEY]['label'], out_png,
        na_metrics={'mechanical_energy_kWh'},
        na_note='PBU uses thermal energy from the glycol\nheating medium, not mechanical work',
    )
