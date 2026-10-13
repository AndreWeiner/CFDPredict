#!/usr/bin/env python3
"""Reduce the delta_p / delta_T function objects to per-cycle statistics.

Pure postProcessing arithmetic, no process control: the worker
(py_microchannel.py) reports the last complete cycle, and the solver loop
(solver.run_case) stops the run once is_steady().
"""
from __future__ import annotations

from glob import glob
from pathlib import Path

import numpy as np

# postProcessing function-object names (case_template/system/controlDict)
FO_DP = "delta_p"                 # multiFieldValue: p(inlet) - p(outlet)
FO_DT = "delta_T"                 # multiFieldValue: T(inlet) - T(outlet)

# Reversal transients are discarded before averaging each half-pass. dp sits
# on a flat plateau (mean over a wide trim); dT drifts monotonically across
# the pass, so it gets a plain mean over a narrower trim.
# Fractions of the period T: the inlet ramp tanh(sin(2*pi*t/T) * Ts) scales
# with T (99 % reached 0.043*T after a reversal at Ts = 10), so a fixed time
# trim would cut into the plateau at high N_rot and keep the ramp at low N_rot.
# The values reproduce the former 0.08 s / 0.05 s at N_rot = 75 (T = 0.8 s).
DP_TRIM_FRAC = 0.10      # [-] of T, discarded at each end of a pass
DT_TRIM_FRAC = 0.0625    # [-] of T

# Quasi-steady stop: kill the solver once each output has changed by less
# than its ABSOLUTE tolerance against the previous cycle.
STEADY_TOL = {"dP_avg_change": 0.01,         # [Pa]
              "dT_forward_change": 0.01,     # [K]
              "dT_reverse_change": 0.01}     # [K]


def series(case, name):
    """Concatenate every postProcessing/<name>/*/*.dat, sorted by time."""
    T, V = [], []
    for f in sorted(glob(f"{case}/postProcessing/{name}/*/*.dat")):
        for ln in open(f):
            if ln.startswith("#"):
                continue
            p = ln.split()
            if len(p) >= 2:
                try:
                    T.append(float(p[0])); V.append(float(p[1]))
                except ValueError:
                    pass
    if not T:
        return None, None
    t = np.array(T); v = np.array(V)
    o = np.argsort(t)
    return t[o], v[o]


def per_pass(t, v, T, trim, reduce_):
    """Apply `reduce_` to each half-pass of every full cycle.

    Cycle k spans [k*T, (k+1)*T): the forward pass is its first half, the
    reverse its second. `trim` seconds are dropped at both ends of each pass
    to exclude the flow-reversal transient.

    Returns (cycle index, forward, reverse)."""
    cyc, fwd, rev = [], [], []
    half = T / 2.0
    k = 0
    while (k + 1) * T <= t.max() + 1e-9:
        f = v[(t >= k * T + trim) & (t < k * T + half - trim)]
        r = v[(t >= k * T + half + trim) & (t < (k + 1) * T - trim)]
        if len(f) > 5 and len(r) > 5:
            cyc.append(k)
            fwd.append(reduce_(f))
            rev.append(reduce_(r))
        k += 1
    return np.array(cyc), np.array(fwd), np.array(rev)


def _change(values: np.ndarray, i: int) -> float | None:
    """Change of cycle i against cycle i-1 (absolute, in the value's units);
    None without a previous one."""
    if i < 1:
        return None
    return float(values[i] - values[i - 1])


def steady_cycle_stats(case_dir: Path, period: float) -> dict:
    """dP_avg / dT_forward / dT_reverse over the last complete cycle, plus
    that cycle's change against the one before (is it settled?). All None
    when the run did not complete a cycle."""
    out = dict(dP_avg=None, dP_forward=None, dP_reverse=None,
               dT_forward=None, dT_reverse=None, steady_cycle=None,
               n_cycles=0, t_end=None,
               dP_avg_change=None, dT_forward_change=None,
               dT_reverse_change=None)
    t_p, v_p = series(case_dir, FO_DP)
    t_T, v_T = series(case_dir, FO_DT)
    if t_p is None or t_T is None:
        return out
    out["t_end"] = float(min(t_p.max(), t_T.max()))
    # dp is +plateau forward, -plateau reverse: average its magnitude
    cp, p_f, p_r = per_pass(t_p, v_p, period, DP_TRIM_FRAC * period,
                            lambda x: float(np.abs(x).mean()))
    cT, T_f, T_r = per_pass(t_T, v_T, period, DT_TRIM_FRAC * period,
                            lambda x: float(x.mean()))
    common = sorted(set(cp.tolist()) & set(cT.tolist()))
    if not common:
        return out
    k = common[-1]
    i = int(np.flatnonzero(cp == k)[0])
    j = int(np.flatnonzero(cT == k)[0])
    p_avg = 0.5 * (p_f + p_r)
    out.update(
        steady_cycle=int(k), n_cycles=len(common),
        dP_forward=float(p_f[i]), dP_reverse=float(p_r[i]),
        dP_avg=float(p_avg[i]),
        dT_forward=float(T_f[j]), dT_reverse=float(T_r[j]),
        dP_avg_change=_change(p_avg, i),
        dT_forward_change=_change(T_f, j),
        dT_reverse_change=_change(T_r, j),
    )
    return out


def outputs(stats: dict, T_hot: float, T_cold: float) -> dict:
    """The reported outputs of steady_cycle_stats(): dP_* as they are, the
    two dT replaced by the effectiveness of the last complete cycle,

        eps_forward = dT_forward / (T_hot - T_cold)   (hot stream)
        eps_reverse = dT_reverse / (T_hot - T_cold)   (cold stream)
        eps         = (eps_forward + eps_reverse) / 2

    Convergence stays tracked on the dT (dT_forward_change /
    dT_reverse_change, the STEADY_TOL criterion)."""
    out = {k: v for k, v in stats.items()
           if k not in ("dT_forward", "dT_reverse")}
    span = T_hot - T_cold
    fwd, rev = stats["dT_forward"], stats["dT_reverse"]
    if fwd is None or rev is None or span <= 0:
        out["eps"] = None
    else:
        eps_forward, eps_reverse = fwd / span, rev / span
        out["eps"] = 0.5 * (eps_forward + eps_reverse)
    return out


def is_steady(stats: dict) -> bool:
    """Every output changed by less than its STEADY_TOL against the previous
    cycle (needs two complete cycles)."""
    return all(stats[k] is not None and abs(stats[k]) < tol
               for k, tol in STEADY_TOL.items())
