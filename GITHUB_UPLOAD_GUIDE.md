# Publish this project on GitHub

**Repository name:** `kc200gt-pv-validation`

**Description:** Reproducible Python/MATLAB solar PV modeling with published I–V data, single/double-diode fitting, NOCT power comparison, and four-fold robustness tests.

## Recommended, no-terminal workflow

1. Download the ZIP and extract it.
2. Install GitHub Desktop: https://desktop.github.com/ and sign in.
3. File → Add Local Repository → choose the extracted `kc200gt-pv-validation` folder.
4. If prompted, Create a Repository Here; accept Git initialization.
5. Commit with `Initial release: KC200GT PV modeling and validation`.
6. Publish repository → uncheck **Keep this code private**.
7. Check the GitHub README image and **Actions** tab.

## Reproduction commands

From the extracted repository root, in a terminal with Python 3.10+:

```bash
python -m venv .venv
# On macOS/Linux: source .venv/bin/activate
# On Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -e '.[dev]'
python -m pytest -q
python -m kc200gt.analysis
python experiments/holdout_robustness.py
```

Expected: 8 tests pass. `results/` and `figures/` are reproducible outputs. Native MATLAB script has not been independently executed in MATLAB.
