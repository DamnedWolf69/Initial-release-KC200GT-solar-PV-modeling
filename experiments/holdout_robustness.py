"""Four predefined interleaved held-out folds, no refitting to held-out samples.
Run after `python -m pip install -e .`: python experiments/holdout_robustness.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from kc200gt.fit import fit_model
from kc200gt.model import current

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'kc200gt_stc_cotfas2021_table21.csv'
OUT = ROOT / 'results' / 'holdout_robustness.csv'

def main():
    data = pd.read_csv(DATA)
    voltage = data['voltage_V'].to_numpy()
    measured = data['current_A'].to_numpy()
    results = []
    for fold in range(4):
        test = np.arange(len(data)) % 4 == fold
        for model in ('SDM', 'DDM'):
            params, info = fit_model(model, voltage[~test], measured[~test])
            predictions = current(voltage[test], params)
            rmse = float(np.sqrt(np.mean((predictions-measured[test])**2)))
            results.append({'fold': fold+1, 'model': model,
                            'train_points': int((~test).sum()),
                            'test_points': int(test.sum()),
                            'test_RMSE_A': rmse})
    OUT.parent.mkdir(exist_ok=True)
    pd.DataFrame(results).to_csv(OUT, index=False)
    print(pd.DataFrame(results).to_string(index=False, float_format=lambda x: f'{x:.8f}'))
    print('\nMean across folds:')
    print(pd.DataFrame(results).groupby('model')['test_RMSE_A'].agg(['mean','std']).to_string())
    print(f'\nSaved to {OUT.relative_to(ROOT)}')

if __name__ == '__main__':
    main()
