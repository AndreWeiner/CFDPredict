"""A3 nozzle -- Zhang Fig.3 Multi-Case U_theta(r) chart.

Sammelt die `swirlU` FunctionObject-Outputs (raw surface format) aus mehreren
A3_nozzle case-Verzeichnissen einer Serie und plottet U_theta(r) am Orifice-
Exit fuer alle Cases gemeinsam (vergleichend) -- Zhang Fig. 3 style.

Drall-Konvention: swirlU schreibt x,y,z,Ux,Uy,Uz auf die orificePlane. Die
Tangentialgeschwindigkeit ist U_theta = (-Ux*y + Uy*x)/r. Da unser
geflipptes Frame (sketch3d.py FLIP x->-x) die Drall-Handedness umkehrt,
plotten wir gegen -U_theta*sign(x) damit der Vergleich mit Zhang Fig. 3
(positive Werte bei +x) direkt geht.

Usage:
    python plot_uTheta_series.py <series_dir> [--zhang-ref <json>] [--out <png>]

<series_dir> enthaelt 1+ Unterverzeichnisse `case_<name>/postProcessing/swirlU/<t>/orificePlane*.raw`.
Default-Output: <series_dir>/progress/zhang_fig3_uTheta_series.png
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _read_swirl_raw(case_dir: Path):
    """Find the latest orificePlane*.raw under postProcessing/swirlU/ and return
    (r [m], U_theta [m/s]) one-dim radial profile.

    Raw format: header line, then x y z Ux Uy Uz columns. We bin by |x|
    (anti-symm dedupe) since the slice is one cell thick around y=0.
    """
    base = case_dir / "postProcessing" / "swirlU"
    if not base.is_dir():
        return None
    times = sorted([d for d in base.iterdir() if d.is_dir()],
                   key=lambda p: float(p.name) if p.name.replace(".", "").isdigit() else 0)
    if not times:
        return None
    raws = list(times[-1].glob("orificePlane*.raw"))
    if not raws:
        return None
    raw = raws[0]
    data = np.loadtxt(raw, comments="#")
    if data.size == 0:
        return None
    x = data[:, 0]
    y = data[:, 1]
    Ux = data[:, 3]
    Uy = data[:, 4]
    r = np.hypot(x, y)
    # U_theta about z-axis: (-Ux*y + Uy*x)/r ; handle r=0 cell-center carefully.
    r_safe = np.where(r > 1e-9, r, 1e-9)
    u_theta = (-Ux * y + Uy * x) / r_safe
    # FLIP-Convention: plot -U_theta*sign(x) so Zhang Fig.3 (+x -> +U_theta) matches.
    s = np.where(x >= 0, 1.0, -1.0)
    u_theta_plot = -u_theta * s
    # Sort by r for line plot.
    idx = np.argsort(r)
    return r[idx], u_theta_plot[idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("series_dir", help="Series dir with case_<name>/ subdirs")
    ap.add_argument("--zhang-ref",
                    default="applications/A3_drallduese/zhang_figure_3_110w_Resolution.json",
                    help="Path to Zhang Fig.3 reference data (JSON)")
    ap.add_argument("--out", default=None, help="Output PNG path")
    ap.add_argument("--title", default="U_theta(r) at orifice exit -- Zhang Fig.3 style")
    args = ap.parse_args()

    series_dir = Path(args.series_dir).resolve()
    if not series_dir.is_dir():
        print("ERROR: %s not a dir" % series_dir, file=sys.stderr)
        sys.exit(1)

    cases = sorted([p for p in series_dir.iterdir()
                    if p.is_dir() and p.name.startswith("case_")])
    if not cases:
        print("ERROR: no case_<name>/ subdirs in %s" % series_dir, file=sys.stderr)
        sys.exit(1)

    fig, ax = plt.subplots(figsize=(9, 6))

    # Zhang reference
    zhang_path = Path(args.zhang_ref)
    if zhang_path.is_file():
        try:
            zdata = json.loads(zhang_path.read_text(encoding="utf-8"))
            # JSON contains a list of {x, y} marker dicts; skip [0] = double-click outlier.
            markers = zdata.get("data", zdata) if isinstance(zdata, dict) else zdata
            if isinstance(markers, list) and len(markers) > 1:
                xs = np.array([m.get("x", m.get("r")) for m in markers[1:]])
                ys = np.array([m.get("y", m.get("U_theta")) for m in markers[1:]])
                ax.plot(xs, ys, "ko", label="Zhang Fig.3 (110W)", markersize=5,
                        markerfacecolor="none")
        except Exception as e:
            print("WARNING: failed to read Zhang ref %s: %s" % (zhang_path, e))

    # Each case as its own line.
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(cases)))
    for c, case in zip(colors, cases):
        prof = _read_swirl_raw(case)
        if prof is None:
            print("WARNING: no swirlU data in %s" % case.name)
            continue
        r, ut = prof
        # Convert r [m] -> mm so axes match Zhang.
        ax.plot(r * 1e3, ut, "-", color=c, lw=1.6,
                label=case.name.replace("case_", ""))

    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xlabel("r [mm]")
    ax.set_ylabel("U_theta [m/s]")
    ax.set_title(args.title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()

    out = Path(args.out) if args.out else (series_dir / "progress" /
                                           "zhang_fig3_uTheta_series.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120)
    print("[plot_uTheta_series] saved %s" % out)


if __name__ == "__main__":
    main()
