#!/usr/bin/env python3
"""Reference Bayes-loop driver for the A2 mother workflow.

Multi-objective optimisation (ParEGO-style scalarisation) of:
    primary:   uniformity at outlet  (maximize)
    secondary: kinematic pressure drop dP  (minimize)

over the 6-dim mother parameter space:
    (vane_r1, vane_r2, vane_r3, vane_ext1, vane_ext2, vane_ext3).

This script is a **reference implementation**, not a production driver.
Adapt freely:
  - Switch to qExpectedHypervolumeImprovement for proper Pareto-aware BO.
  - Replace the penalty with a constraint-aware classifier (Gardner 2014).
  - Submit batches to SLURM instead of running sequentially.
  - Add more design parameters (R, Re, ...) -- edit interface.json + this file.

See BAYES_INTEGRATION.md for the design rationale.

Usage:
    python examples/bayes_skeleton.py --n-init 20 --n-iter 30 \\
                                       --workdir ./runs/bayes_demo

Requirements (in addition to the workflow's requirements.txt):
    pip install botorch gpytorch torch numpy scipy
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

# Make tools/run_workflow.py importable when this script is run from anywhere
HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
sys.path.insert(0, str(REPO_ROOT))

from tools.run_workflow import run_workflow_headless

# Defer heavy imports so --help works without botorch installed
def _require_botorch():
    try:
        import torch
        import botorch
        import gpytorch
        import numpy as np
        from scipy.stats import qmc
        return torch, botorch, gpytorch, np, qmc
    except ImportError as e:
        print("This skeleton needs botorch + gpytorch + torch + numpy + scipy.",
              file=sys.stderr)
        print(f"  pip install botorch gpytorch torch numpy scipy", file=sys.stderr)
        print(f"  (missing: {e.name})", file=sys.stderr)
        sys.exit(2)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
@dataclass
class BayesConfig:
    # Geometry that stays fixed across the loop (matches mother defaults)
    H: float = 1.0
    W: float = 1.0
    R: float = 1.5
    L_in: float = 3.0
    L_out: float = 5.0
    # Parameter space (per-vane)
    N_max: int = 3
    eps: float = 1e-3                # margin from the wall
    # Loop config
    n_init: int = 20                 # Latin-hypercube initial design size
    n_iter: int = 30                 # post-init Bayes iterations
    seed: int = 0
    timeout_s: float = 1800.0        # per-case
    # Penalties for failed cases (see BAYES_INTEGRATION.md)
    penalty_uniformity: float = 0.5  # well below typical 0.9+
    penalty_dP: float = 10.0         # ~25x typical default dP=0.4

    @property
    def inner_r(self) -> float: return self.R - self.H / 2 + self.eps
    @property
    def outer_r(self) -> float: return self.R + self.H / 2 - self.eps

    def to_bounds(self):
        """Return torch.Tensor of shape (2, d) -- botorch bounds convention."""
        import torch
        d = 2 * self.N_max
        lo = [self.inner_r] * self.N_max + [0.0] * self.N_max
        hi = [self.outer_r] * self.N_max + [self.L_out] * self.N_max
        return torch.tensor([lo, hi], dtype=torch.double)


# ---------------------------------------------------------------------------
# Settings construction
# ---------------------------------------------------------------------------
def base_settings(cfg: BayesConfig) -> dict:
    """Mother-workflow settings dict with everything except the design vars."""
    return {
        "project name": "bayes_loop",
        "channel height H [m]": cfg.H,
        "channel depth W [m]": cfg.W,
        "bend center radius R [m]": cfg.R,
        "inlet length L_in [m]": cfg.L_in,
        "outlet length L_out [m]": cfg.L_out,
        "number of vanes": cfg.N_max,
        "vane_r1 [m]": cfg.R + cfg.H / 4,         # placeholder; will be overwritten
        "vane_r2 [m]": cfg.R,
        "vane_r3 [m]": cfg.R - cfg.H / 4,
        "vane_ext1 [m]": 0.0,
        "vane_ext2 [m]": 0.0,
        "vane_ext3 [m]": 0.0,
        "Reynolds number": 1.0e5,
        "turbulent intensity": 0.05,
        "cells inlet axial nx_in": 60,
        "cells bend angular nx_bend": 80,
        "cells outlet axial nx_out": 100,
        "cells across channel ny": 40,
        "cells in z (nz)": 1,                     # quasi-2D for cheap demo
        "auto layers": "true",
        "y+ target": 40.0,
        "wall-layer splits": 0,
        "layer-split thickness fraction": 0.5,
        "number of processors": 4,
        "contact email": "",
        "comments": "Bayes-loop iteration (see iterations.csv).",
    }


def x_to_settings(x, cfg: BayesConfig, base: dict) -> dict:
    """Inject design vector x = (r1, r2, r3, ext1, ext2, ext3) into base settings.
    The worker validates ordering / range; we don't re-check here."""
    s = dict(base)
    radii = [float(v) for v in x[:cfg.N_max]]
    exts = [float(v) for v in x[cfg.N_max:]]
    # Sort radii descending so vane_r1 > vane_r2 > vane_r3 (mother convention).
    # Ext indices follow the (unsorted) acquisition order -- this means the GP
    # sees the ext-i variables as exchangeable across vanes. If you want
    # per-position ext semantics, drop the sort here and let the GP learn it.
    radii_sorted = sorted(radii, reverse=True)
    for i, r in enumerate(radii_sorted, 1):
        s[f"vane_r{i} [m]"] = r
    for i, e in enumerate(exts, 1):
        s[f"vane_ext{i} [m]"] = e
    return s


# ---------------------------------------------------------------------------
# Evaluation: workflow -> (uniformity, dP)
# ---------------------------------------------------------------------------
def evaluate(x, cfg: BayesConfig, base: dict, iter_id: int,
             workdir_root: Path) -> dict:
    """Run one workflow case. Returns dict with keys for the iterations log."""
    settings = x_to_settings(x, cfg, base)
    settings["project name"] = f"bayes_iter_{iter_id:04d}"
    workdir = workdir_root / f"iter_{iter_id:04d}"
    t0 = time.time()
    result = run_workflow_headless(
        "A2_leitbleche", settings, workdir,
        tail=False, timeout_s=cfg.timeout_s,
    )
    wall = time.time() - t0

    uniformity = dP = vorticity = iters = None
    if result["exit_code"] == 0 and result.get("summary"):
        cases = result["summary"].get("cases", [])
        if cases:
            c = cases[0]
            uniformity = c.get("uniformity")
            dP = c.get("vorticity")  # placeholder; corrected below
            dP = c.get("dP")
            vorticity = c.get("vorticity")
            iters = c.get("iters")

    failed = (uniformity is None or dP is None)
    if failed:
        uniformity = cfg.penalty_uniformity
        dP = cfg.penalty_dP

    return {
        "iter_id": iter_id,
        "vane_r1": settings["vane_r1 [m]"], "vane_r2": settings["vane_r2 [m]"],
        "vane_r3": settings["vane_r3 [m]"], "vane_ext1": settings["vane_ext1 [m]"],
        "vane_ext2": settings["vane_ext2 [m]"], "vane_ext3": settings["vane_ext3 [m]"],
        "uniformity": uniformity, "dP": dP,
        "vorticity": vorticity, "iters": iters,
        "wallclock_s": wall, "exit_code": result["exit_code"],
        "reason": result.get("reason", ""),
        "penalized": failed,
        "workdir": str(workdir),
    }


# ---------------------------------------------------------------------------
# Initial design (Latin-hypercube + ordering projection)
# ---------------------------------------------------------------------------
def initial_design(cfg: BayesConfig):
    _, _, _, np, qmc = _require_botorch()
    import torch
    lhs = qmc.LatinHypercube(d=2 * cfg.N_max, seed=cfg.seed)
    raw = lhs.random(n=cfg.n_init)
    # Map to physical bounds (radii get sorted in x_to_settings, so any
    # ordering of the first N_max columns is acceptable here).
    lo, hi = cfg.to_bounds().numpy()
    samples = lo[None, :] + raw * (hi - lo)[None, :]
    return torch.tensor(samples, dtype=torch.double)


# ---------------------------------------------------------------------------
# Botorch ParEGO-style acquisition
# ---------------------------------------------------------------------------
def fit_gp_and_propose(X, Y, cfg: BayesConfig):
    """Fit a multi-output GP, scalarise with random weights (Tchebycheff /
    ParEGO), maximise EI to pick the next point. Returns x_next as torch
    tensor of shape (d,)."""
    torch, botorch, gpytorch, np, _ = _require_botorch()
    from botorch.models import SingleTaskGP
    from botorch.fit import fit_gpytorch_mll
    from botorch.acquisition import ExpectedImprovement
    from botorch.optim import optimize_acqf
    from gpytorch.mlls import ExactMarginalLogLikelihood

    # Y is (n, 2): col 0 = uniformity (maximize), col 1 = dP (minimize).
    # ParEGO trick: random Tchebycheff weights -> scalar surrogate.
    # We minimise -uniformity + dP for the EI fit; sign-flip uniformity.
    y_flipped = Y.clone()
    y_flipped[:, 0] = -y_flipped[:, 0]            # now both "to be minimised"

    # Normalise so weights are scale-invariant
    y_norm = (y_flipped - y_flipped.mean(dim=0)) / (y_flipped.std(dim=0) + 1e-9)
    rng = np.random.default_rng(cfg.seed + len(X))
    w = rng.dirichlet([1.0, 1.0])                  # random simplex weight
    rho = 0.05
    # Augmented Tchebycheff
    aug = (w * y_norm).max(dim=1).values + rho * (w * y_norm).sum(dim=1)

    gp = SingleTaskGP(X, aug.unsqueeze(-1).double())
    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
    fit_gpytorch_mll(mll)

    # EI on -aug (we want to minimise aug → maximise -aug)
    best_neg_aug = -aug.min().item()
    acq = ExpectedImprovement(model=gp, best_f=best_neg_aug, maximize=False)

    bounds = cfg.to_bounds()
    x_next, _ = optimize_acqf(
        acq, bounds=bounds, q=1, num_restarts=10, raw_samples=256,
    )
    return x_next.squeeze(0)


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
LOG_COLUMNS = [
    "iter_id",
    "vane_r1", "vane_r2", "vane_r3",
    "vane_ext1", "vane_ext2", "vane_ext3",
    "uniformity", "dP", "vorticity", "iters",
    "wallclock_s", "exit_code", "reason", "penalized", "workdir",
]


def init_log(path: Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerow(LOG_COLUMNS)


def append_log(path: Path, row: dict) -> None:
    with open(path, "a", newline="", encoding="utf-8") as f:
        csv.writer(f).writerow([row.get(k, "") for k in LOG_COLUMNS])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-init", type=int, default=20)
    ap.add_argument("--n-iter", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--workdir", type=Path, default=Path("./runs/bayes_demo"))
    ap.add_argument("--timeout", type=float, default=1800.0,
                    help="Per-case timeout in seconds.")
    ap.add_argument("--Re", type=float, default=1.0e5)
    ap.add_argument("--nz", type=int, default=1,
                    help="1 = quasi-2D (cheap demo). 20 = 3D (production).")
    args = ap.parse_args()

    cfg = BayesConfig(
        n_init=args.n_init, n_iter=args.n_iter,
        seed=args.seed, timeout_s=args.timeout,
    )
    base = base_settings(cfg)
    base["Reynolds number"] = args.Re
    base["cells in z (nz)"] = args.nz

    args.workdir.mkdir(parents=True, exist_ok=True)
    log_path = args.workdir / "iterations.csv"
    init_log(log_path)
    (args.workdir / "config.json").write_text(
        json.dumps({"cfg": cfg.__dict__, "base_settings": base}, indent=2),
        encoding="utf-8")

    torch, _, _, np, _ = _require_botorch()

    print(f"[bayes] Initial design: {cfg.n_init} Latin-hypercube samples")
    X = initial_design(cfg)
    Y_rows = []
    for i, x in enumerate(X):
        print(f"[bayes] iter {i:04d} (init) ", end="", flush=True)
        row = evaluate(x.tolist(), cfg, base, iter_id=i, workdir_root=args.workdir)
        Y_rows.append((row["uniformity"], row["dP"]))
        append_log(log_path, row)
        print(f"gamma={row['uniformity']:.4f}  dP={row['dP']:.4f}  "
              f"[{row['wallclock_s']:.1f}s, reason={row['reason']}]")

    Y = torch.tensor(Y_rows, dtype=torch.double)

    for k in range(cfg.n_iter):
        i = cfg.n_init + k
        x_next = fit_gp_and_propose(X, Y, cfg)
        print(f"[bayes] iter {i:04d} (EI)   ", end="", flush=True)
        row = evaluate(x_next.tolist(), cfg, base, iter_id=i,
                       workdir_root=args.workdir)
        X = torch.cat([X, x_next.unsqueeze(0)], dim=0)
        Y = torch.cat([Y, torch.tensor([[row["uniformity"], row["dP"]]],
                                       dtype=torch.double)], dim=0)
        append_log(log_path, row)
        print(f"gamma={row['uniformity']:.4f}  dP={row['dP']:.4f}  "
              f"[{row['wallclock_s']:.1f}s, reason={row['reason']}]")

    # Pareto front summary
    Y_np = Y.numpy()
    # uniformity is to-maximize, dP is to-minimize
    pareto_mask = []
    for i in range(len(Y_np)):
        dominated = False
        for j in range(len(Y_np)):
            if i == j: continue
            if Y_np[j, 0] >= Y_np[i, 0] and Y_np[j, 1] <= Y_np[i, 1] \
               and (Y_np[j, 0] > Y_np[i, 0] or Y_np[j, 1] < Y_np[i, 1]):
                dominated = True
                break
        pareto_mask.append(not dominated)

    print()
    print(f"[bayes] done. {sum(pareto_mask)} / {len(Y_np)} Pareto-optimal points.")
    print(f"[bayes] best uniformity: {Y_np[:, 0].max():.4f}")
    print(f"[bayes] best dP:         {Y_np[:, 1].min():.4f}")
    print(f"[bayes] iterations.csv:  {log_path}")
    print(f"[bayes] full run dirs:   {args.workdir}/iter_XXXX/")


if __name__ == "__main__":
    main()
