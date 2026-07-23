from pathlib import Path

import pandas as pd

from sensitivity_common import PARAM_META, run_case


base_params = {
    name: meta['base']
    for name, meta in PARAM_META.items()
}

# for system_key in ('PBU', 'Pump', 'Compressor'):
for system_key in ('Pump', 'Compressor'):

    print(f'Running {system_key} baseline...')

    profile_path = Path(
        f'results/sensitivity/test_{system_key}_baseline.csv'
    )

    metrics = run_case(
        system_key,
        base_params,
        profile_path=profile_path,
    )

    for key, value in metrics.items():
        print(f'{key}: {value}')
    if system_key == 'PBU':
        df = pd.read_csv(profile_path)

        idx = df['pbu_liquid_reynolds'].idxmax()

        print(
            'Maximum PBU Reynolds number:',
            df.loc[idx, 'pbu_liquid_reynolds']
        )
        print(
            'Occurred at day:',
            df.loc[idx, 'day']
        )
        print(
            'Tank pressure [bar]:',
            df.loc[idx, 'pressure_bar']
        )

    print()


