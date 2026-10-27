"""Deterministic multistart nonlinear least-squares parameter estimation."""
import numpy as np
from scipy.optimize import least_squares
from .model import decode, current

# Optimization uses log10(I0) and log10(Rsh) to avoid scaling problems.
CONFIG = {
    'SDM': {
        'lower': [7.8, -13.0, 1.0, 0.0, 1.0],
        'upper': [8.6, -4.0, 2.5, 0.8, 4.0],
        'seeds': [
            [8.19, -9.0, 1.1, 0.2, 2.2],
            [8.19, -8.5, 1.3, 0.3, 2.3],
            [8.19, -7.0, 1.7, 0.2, 2.0]],
    },
    'DDM': {
        'lower': [7.8, -13.0, -13.0, 1.0, 1.5, 0.0, 1.0],
        'upper': [8.6, -4.0, -4.0, 1.5, 3.0, 0.8, 4.0],
        'seeds': [
            [8.19, -9.0, -9.0, 1.05, 2.0, 0.2, 2.2],
            [8.19, -10.0, -8.0, 1.1, 2.2, 0.25, 2.2],
            [8.19, -8.8, -11.0, 1.1, 2.5, 0.3, 2.0],
            [8.19, -9.0, -10.0, 1.3, 1.6, 0.25, 2.4]],
    },
}


def fit_model(kind: str, V: np.ndarray, I: np.ndarray):
    config = CONFIG[kind]
    solutions = []
    for seed in config['seeds']:
        result = least_squares(
            lambda x: current(V, decode(x, kind)) - I,
            x0=seed, bounds=(config['lower'], config['upper']),
            method='trf', x_scale='jac', max_nfev=6000,
            ftol=2e-12, xtol=2e-12, gtol=2e-12)
        if not result.success and not result.cost < 1e-3:
            raise RuntimeError(f'Optimizer did not converge: {kind}: {result.message}')
        solutions.append(result)
    best = min(solutions, key=lambda r: r.cost)
    return decode(best.x, kind), {'optimization_coordinates': best.x.tolist(),
                                 'training_SSE_A2': float(2 * best.cost),
                                 'n_starts': len(solutions),
                                 'optimizer_message': best.message}
