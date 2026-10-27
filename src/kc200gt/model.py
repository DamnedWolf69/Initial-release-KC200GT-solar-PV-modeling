"""Single/double-diode IV models for a 54-cell KC200GT module.

The numerical current solver uses a damped Newton iteration and only uses a
bracketed scalar fallback if Newton fails. All currents are module currents.
Temperature means *cell* temperature, not ambient temperature.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq, minimize_scalar

K_BOLTZMANN = 1.380649e-23  # J/K
Q_ELECTRON = 1.602176634e-19  # C
K_EV = K_BOLTZMANN / Q_ELECTRON  # eV/K


@dataclass(frozen=True)
class Module:
    cells: int = 54
    irradiance_ref: float = 1000.0
    cell_temp_ref_C: float = 25.0
    alpha_isc_A_per_C: float = 3.18e-3
    eg_ref_eV: float = 1.12  # silicon approximation, NOT a KC200GT-fitted value
    voc_datasheet_V: float = 32.9
    isc_datasheet_A: float = 8.21
    vmp_datasheet_V: float = 26.3
    imp_datasheet_A: float = 7.61
    noct_cell_temp_C: float = 47.0
    noct_irradiance: float = 800.0
    noct_pmax_datasheet_W: float = 142.0  # directly quoted, rounded nameplate
    noct_vmp_datasheet_V: float = 23.2
    noct_imp_datasheet_A: float = 6.13


MODULE = Module()


@dataclass(frozen=True)
class Params:
    iph_A: float
    i0_A: tuple[float, ...]
    n: tuple[float, ...]
    rs_ohm: float
    rsh_ohm: float


def decode(x: np.ndarray, kind: str) -> Params:
    x = np.asarray(x, dtype=float)
    if kind == 'SDM':
        assert x.size == 5
        return Params(x[0], (10.0 ** x[1],), (x[2],), x[3], 10.0 ** x[4])
    if kind == 'DDM':
        assert x.size == 7
        return Params(x[0], (10.0 ** x[1], 10.0 ** x[2]),
                      (x[3], x[4]), x[5], 10.0 ** x[6])
    raise ValueError('kind must be SDM or DDM')


def at_condition(ref: Params, G: float, Tcell_C: float,
                 module: Module = MODULE) -> Params:
    if G <= 0 or Tcell_C <= -273.15:
        raise ValueError('Require positive irradiance and Kelvin temperature')
    T = Tcell_C + 273.15
    Tref = module.cell_temp_ref_C + 273.15
    EgT = module.eg_ref_eV * (1.0 - 0.0002677 * (Tcell_C - module.cell_temp_ref_C))
    iph = (G / module.irradiance_ref) * (
        ref.iph_A + module.alpha_isc_A_per_C * (Tcell_C - module.cell_temp_ref_C))
    i0 = tuple(i0ref * (T / Tref) ** 3 *
               np.exp((module.eg_ref_eV / Tref - EgT / T) / (ni * K_EV))
               for i0ref, ni in zip(ref.i0_A, ref.n))
    return Params(iph, i0, ref.n, ref.rs_ohm,
                  ref.rsh_ohm * module.irradiance_ref / G)


def _voltage_f(I: np.ndarray, V: np.ndarray, p: Params, Tcell_C: float,
               module: Module) -> tuple[np.ndarray, np.ndarray]:
    T = Tcell_C + 273.15
    vd = V + I * p.rs_ohm
    F = p.iph_A - I - vd / p.rsh_ohm
    derivative = np.zeros_like(I) - (1.0 + p.rs_ohm / p.rsh_ohm)
    for i0, ni in zip(p.i0_A, p.n):
        a = ni * module.cells * K_BOLTZMANN * T / Q_ELECTRON
        ez = np.exp(np.clip(vd / a, -100.0, 100.0))
        F -= i0 * (ez - 1.0)
        derivative -= i0 * ez * p.rs_ohm / a
    return F, derivative


def current(V: np.ndarray | float, ref: Params, G: float = 1000.0,
            Tcell_C: float = 25.0, module: Module = MODULE) -> np.ndarray | float:
    v = np.asarray(V, dtype=float)
    shape = v.shape
    vf = np.atleast_1d(v).ravel()
    p = at_condition(ref, G, Tcell_C, module)
    I = np.full(vf.shape, p.iph_A, dtype=float) - vf / p.rsh_ohm
    for _ in range(90):
        F, dF = _voltage_f(I, vf, p, Tcell_C, module)
        step = np.clip(F / dF, -4.0, 4.0)
        I -= step
        if np.max(np.abs(step)) < 1e-12:
            break
    F, _ = _voltage_f(I, vf, p, Tcell_C, module)
    bad = ~np.isfinite(I) | ~np.isfinite(F) | (np.abs(F) > 1e-7)
    for idx in np.where(bad)[0]:
        # F(I) is strictly decreasing; these bounds bracket the physical root.
        fn = lambda ii: float(_voltage_f(np.array([ii]), np.array([vf[idx]]),
                                         p, Tcell_C, module)[0][0])
        low, high = -100.0, max(100.0, p.iph_A + 100.0)
        I[idx] = brentq(fn, low, high, xtol=1e-12)
    out = I.reshape(shape)
    return float(out) if v.ndim == 0 else out


def _open_circuit_f(V: float, p: Params, Tcell_C: float,
                    module: Module) -> float:
    T = Tcell_C + 273.15
    f = p.iph_A - V / p.rsh_ohm
    for i0, ni in zip(p.i0_A, p.n):
        a = ni * module.cells * K_BOLTZMANN * T / Q_ELECTRON
        f -= i0 * np.expm1(np.clip(V / a, -100.0, 100.0))
    return float(f)


def mpp(ref: Params, G: float = 1000.0, Tcell_C: float = 25.0,
        module: Module = MODULE) -> dict[str, float]:
    p = at_condition(ref, G, Tcell_C, module)
    voc = brentq(lambda v: _open_circuit_f(v, p, Tcell_C, module),
                 0.0, 80.0, xtol=1e-11)
    solution = minimize_scalar(lambda v: -v * current(v, ref, G, Tcell_C, module),
                               method='bounded', bounds=(0, voc),
                               options={'xatol': 1e-9})
    if not solution.success:
        raise RuntimeError('Pmax optimization failed')
    vmp = float(solution.x)
    imp = float(current(vmp, ref, G, Tcell_C, module))
    isc = float(current(0.0, ref, G, Tcell_C, module))
    return {'Pmax_W': vmp * imp, 'Vmp_V': vmp, 'Imp_A': imp,
            'Voc_V': voc, 'Isc_A': isc, 'FF': vmp * imp / (voc * isc)}
