# -*- coding: utf-8 -*-
"""Run the dual-fuel OAT sensitivity analysis for all three FGSS designs."""
import sys
from pathlib import Path

ANALYSIS_DIR = Path(__file__).resolve().parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

import sensitivity_common as sc

NA_METRICS_BY_SYSTEM = {
    'PBU': {
        'mechanical_energy_kWh',
        'specific_mechanical_energy_kWh_per_t',
    },
}

NA_NOTES_BY_SYSTEM = {
    'PBU': (
        'Not applicable: FGSS 1 has no pump or compressor. '
        'The PBU requires thermal energy from the heating medium.'
    ),
}
if __name__ == '__main__':
    for system_key in ('PBU', 'Pump', 'Compressor'):
        df, base_metrics = sc.run_oat_sensitivity(system_key)
        out_png = sc.RESULTS_DIR / f'{system_key}_dual_tornado.png'
        sc.plot_tornado(
            df,
            base_metrics,
            sc.SYSTEM_SPECS[system_key]['label'],
            out_png,
            na_metrics=NA_METRICS_BY_SYSTEM.get(system_key),
            na_note=NA_NOTES_BY_SYSTEM.get(
                system_key,
                'Not applicable for this FGSS',
            ),
        )
