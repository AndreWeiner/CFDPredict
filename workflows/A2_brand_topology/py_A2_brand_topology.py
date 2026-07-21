#!/usr/bin/env python3
"""A2 brand_topology worker -- parallel-translated 90deg bend with full-
height guide vanes (CFDPredict A2 Visualisierungs-Variante).

Per case:  blockMesh -> (topoSet + createBaffles when vanes) ->
refineWallLayer -> potentialFoam -> simpleFoam (k-omega-SST) via the
case_template Allrun (mesh-only) + worker-driven solver.

Single case OR mini-series: both "number of vanes" and "cells across
inlet N_inlet" accept lists (e.g. "[0, 3, 5]" or "[12, 16, 20]") -> one
case per entry (cartesian product when both are lists).

Worker contract (WORKFLOW_AUTHORING.md section 4):
    python py_A2_brand_topology.py <series_dir>
  - series_dir contains interface.json (6-tuple, values injected because
    inject_schema: true), progress/, owner, workflow.
  - progress/<step>_info.txt for the live log; result PNG into progress/.
  - touch <series_dir>/command_finished as the very last step.
  - honour <series_dir>/command_kill: exit cleanly + touch command_finished.

NOTE: OpenFOAM is NOT installed locally; the build/parse logic is
offline-testable (see test_build_offline.py), but the actual solver run
is exercised only on the VM with OF v2512.
"""
from __future__ import annotations

import ast
import os
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import foam_dictionary
import py_OF_utils
from py_OF_utils import (touch, run_OF_utility, run_solver_copy_progress,
                         clean_processor_directories, get_latestTime,
                         run_pvbatch_utility)

from gen_blockmesh import BrandTopology
from build_case import build

# --- fixed physics, matching the A2 default --------------------------------
NU = 1.5e-5            # m^2/s, air 20 C
# OpenFOAM bashrc is resolved lazily via py_OF_utils.resolve_openfoam_bashrc()
# at solver-spawn time, not at import time -- so the workflow can also be
# imported on machines without OpenFOAM (offline tests, Bayes-loop drivers
# that only need the build half).


# ---------------------------------------------------------------------------
# input handling
# ---------------------------------------------------------------------------
def check_args() -> str:
    if len(sys.argv) != 2:
        print("Usage:")
        print("python %s /path/to/series_dir" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    return sys.argv[1]


def json_to_unrolled(series_path: str, json_fname: str) -> str:
    import json
    with open(os.path.join(series_path, json_fname), "r", encoding="utf-8") as f:
        json_input = json.load(f)
    json_unrolled = foam_dictionary.unroll_dict(json_input)
    out_name = json_fname.replace(".json", "Unrolled2.json")
    with open(os.path.join(series_path, out_name), "w", encoding="utf-8") as f:
        f.write(json.dumps(json_unrolled))
    return out_name


def read_unrolled(series_path: str, fname: str) -> dict:
    import json
    with open(os.path.join(series_path, fname), "r", encoding="utf-8") as f:
        return json.load(f)


def parse_int_series(raw, min_val: int = 1, fallback: int | None = None):
    """Parse a field that may be a single number, a JSON-ish list or a
    comma/space-separated list. Returns ints >= min_val. Mirrors the
    zetaBrand parser."""
    s = str(raw).strip()
    out: list[int] = []
    if s:
        try:
            val = ast.literal_eval(s)
        except (ValueError, SyntaxError, NameError):
            val = None
        if isinstance(val, (list, tuple)):
            for v in val:
                try:
                    out.append(int(round(float(v))))
                except (ValueError, TypeError):
                    pass
        elif isinstance(val, (int, float)):
            out = [int(round(float(val)))]
        else:
            toks = (s.replace("[", " ").replace("]", " ")
                     .replace(",", " ").split())
            for t in toks:
                try:
                    out.append(int(round(float(t))))
                except ValueError:
                    pass
    out = [v for v in out if v >= min_val]
    if out:
        return out
    return [fallback if fallback is not None else min_val]


# Each vane spec is either:
#   int                            -> N equidistant vanes (Mode 1)
#   list[float]                    -> explicit fractions, one case (Mode 2)
#   ("P", int) | ("P", list[int])  -> Brand-Pos single-vane case(s) (Mode 3)
VaneSpec = int  # type: ignore  # documentation alias

# Brand 2020 wind-tunnel convention (S. 114): channel side length 4 m,
# vane-position grid 100 mm. Pos N -> distance from inner wall = (N+1)*100 mm
# -> fraction (in normalised H=1) = (N+1)/40.
BRAND_H_PHYS_MM = 4000.0
BRAND_POS_STEP_MM = 100.0


def brand_pos_to_fraction(pos: int) -> float:
    """Brand-Pos index -> Pathline-Bruchteil in (0, 1).
    Pos N = (N+1)*100 mm Abstand zur Innenwand bei H_phys=4 m.
    """
    f = (pos + 1) * BRAND_POS_STEP_MM / BRAND_H_PHYS_MM
    return round(float(f), 8)


def parse_vane_specs(raw) -> list:
    """Parse the 'number of vanes' field into a list of per-case vane specs.

    Mode 1 (equidistant counts) -- backwards-compatible:
        "3"          -> [3]                     # one case, 3 equidist vanes
        "[0, 3, 5]"  -> [0, 3, 5]               # 3-case sweep over N
        "0"          -> [0]                     # empty bend

    Mode 2 (explicit fractions, single case):
        "[0.27, 0.5, 0.7]"  -> [[0.27, 0.5, 0.7]]   # one case, 3 explicit

    Mode 3 (Brand-Pos, prefix 'P'):
        "P11"               -> [("P", 11)]                    # 1 case, 1 vane at Pos 11
        "P[1, 5, 11]"       -> [("P", 1), ("P", 5), ("P", 11)]  # 3-case Pos sweep, each w/ 1 vane
        Brand Fig. 7.14/7.15 single-vane positions study.

    Detection rules:
      - String starts with 'P' (or 'p') followed by a digit / '[' -> Mode 3.
      - List with any non-integer float in (0, 1) -> Mode 2.
      - Otherwise -> Mode 1 (counts).
    """
    s = str(raw).strip()
    if not s:
        return [0]

    # Mode 3: Brand-Pos prefix 'P'
    if s[0] in ("P", "p") and (len(s) > 1 and (s[1].isdigit() or s[1] == "[")):
        body = s[1:].strip()
        try:
            val = ast.literal_eval(body)
        except (ValueError, SyntaxError, NameError):
            val = None
        if isinstance(val, (list, tuple)):
            return [("P", int(round(float(v)))) for v in val
                    if int(round(float(v))) >= 0]
        if isinstance(val, (int, float)):
            return [("P", int(round(float(val))))]
        # fallback: split tokens like "P 1 5 11"
        toks = (body.replace("[", " ").replace("]", " ")
                    .replace(",", " ").split())
        specs: list = []
        for t in toks:
            try:
                specs.append(("P", int(round(float(t)))))
            except ValueError:
                pass
        return specs if specs else [0]

    try:
        val = ast.literal_eval(s)
    except (ValueError, SyntaxError, NameError):
        val = None

    def _is_fractional(x) -> bool:
        try:
            f = float(x)
        except (TypeError, ValueError):
            return False
        return f != int(f)

    if isinstance(val, (list, tuple)):
        if any(_is_fractional(v) for v in val):
            fractions = sorted({round(float(v), 8) for v in val
                                if 0.0 < float(v) < 1.0})
            return [list(fractions)] if fractions else [0]
        # all integer-valued -> counts series
        counts = [int(round(float(v))) for v in val
                  if int(round(float(v))) >= 0]
        return counts if counts else [0]

    if isinstance(val, (int, float)):
        if _is_fractional(val) and 0.0 < float(val) < 1.0:
            return [[float(val)]]
        n = int(round(float(val)))
        return [max(0, n)]

    # manual fallback for malformed input -- try comma/space split,
    # interpret as integer counts
    toks = (s.replace("[", " ").replace("]", " ")
             .replace(",", " ").split())
    counts: list[int] = []
    for t in toks:
        try:
            n = int(round(float(t)))
            if n >= 0:
                counts.append(n)
        except ValueError:
            pass
    return counts if counts else [0]


def update_progress(progress_path: str, info_list: list[str],
                    step: str = "1") -> None:
    try:
        with open(os.path.join(progress_path, "%s_info.txt" % step),
                  "w", encoding="utf-8") as f:
            for line in info_list:
                f.write("%s\n" % line)
    except Exception as e:
        print("update_progress failed (non-fatal): %s" % e)


def killed(series_path: str) -> bool:
    return os.path.isfile(os.path.join(series_path, "command_kill"))


# ---------------------------------------------------------------------------
# per-case build + run + parse
# ---------------------------------------------------------------------------
def equidistant_fractions(N: int) -> list[float]:
    """Pathline-Bruchteile in (0,1) fuer N gleichmaessig verteilte Bleche.
    Konvention aus dem brand_topology-Handoff:  i/(N+1) fuer i=1..N."""
    if N <= 0:
        return []
    return [i / (N + 1) for i in range(1, N + 1)]


def resolve_vane_spec(spec) -> list[float]:
    """Spec -> sorted list of pathline fractions (inner=0 .. outer=1).
        int            -> equidistant counts
        list[float]    -> explicit fractions, filtered to (0,1), sorted
        ("P", int)     -> Brand-Pos single vane: [brand_pos_to_fraction(N)]
        ("P", list)    -> Brand-Pos: multiple vanes from list of pos values
    """
    if isinstance(spec, tuple) and len(spec) == 2 and spec[0] == "P":
        pos_val = spec[1]
        if isinstance(pos_val, (list, tuple)):
            fr = sorted({brand_pos_to_fraction(int(p)) for p in pos_val
                         if 0.0 < brand_pos_to_fraction(int(p)) < 1.0})
            return fr
        f = brand_pos_to_fraction(int(pos_val))
        return [f] if 0.0 < f < 1.0 else []
    if isinstance(spec, (list, tuple)):
        return sorted({round(float(f), 8) for f in spec
                       if 0.0 < float(f) < 1.0})
    return equidistant_fractions(int(spec))


def build_case_dir(case_dir: Path, *, H: float, W: float, R: float,
                   L_in: float, L_out: float, L_trail: float,
                   fractions: list[float], N_inlet: int, nz: int,
                   Re: float, ti: float, n_proc: int,
                   n_layer_splits: int, layer_thickness: float) -> dict:
    """Build one OF case via the brand_topology build() pipeline.
    `fractions` is the explicit list of Pathline-Bruchteile for the case."""
    geo = BrandTopology(
        H=H, W=W, R=R, L_in=L_in, L_trail=L_trail, L_out=L_out,
        r_bend=None,
        inner_arc_fractions=list(fractions),
        N_inlet=N_inlet, nz=nz,
    )
    info = build(case_dir, geo,
                 Re=Re, nu=NU, turbulent_intensity=ti,
                 n_layer_splits=n_layer_splits,
                 layer_thickness=layer_thickness)

    # set the decomposition rank count (worker runs mpirun -np n_proc)
    import re as _re
    dpd = case_dir / "system" / "decomposeParDict"
    try:
        txt = dpd.read_text(encoding="utf-8")
        txt = _re.sub(r"numberOfSubdomains\s+\d+\s*;",
                      "numberOfSubdomains  %d;" % n_proc, txt)
        dpd.write_text(txt, encoding="utf-8", newline="\n")
    except Exception as e:
        print("WARNING: cannot patch decomposeParDict in %s: %s" % (dpd, e))

    return info


def parse_dat_last(case_dir: Path, fo_name: str, col: int = 1) -> float | None:
    """Read the last numeric value in the given column of a surfaceFieldValue
    log (column 0 = time). Returns None if missing."""
    dat = (case_dir / "postProcessing" / fo_name / "0"
           / "surfaceFieldValue.dat")
    if not dat.is_file():
        return None
    try:
        rows = [ln for ln in dat.read_text(encoding="utf-8").splitlines()
                if ln.strip() and not ln.lstrip().startswith("#")]
        if not rows:
            return None
        parts = rows[-1].split()
        if len(parts) <= col:
            return None
        return float(parts[col])
    except Exception as e:
        print("parse_dat_last(%s,%d) failed: %s" % (fo_name, col, e))
        return None


def parse_iters(case_dir: Path) -> int | None:
    try:
        lt = get_latestTime(str(case_dir))
        if lt is not None and float(lt) > 0:
            return int(round(float(lt)))
    except Exception:
        pass
    log = case_dir / "logs" / "solver.log"
    if not log.is_file():
        return None
    import re
    try:
        text = log.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None
    m = re.search(r"converged in (\d+) iterations", text, re.IGNORECASE)
    if m:
        return int(m.group(1))
    times = re.findall(r"^Time = (\d+)", text, re.MULTILINE)
    return int(times[-1]) if times else None


def _log_tail(path: Path, n_chars: int = 1500) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")[-n_chars:]
    except OSError:
        return f"(log unreadable: {path})"


def _newest_log_tail(case_dir: Path) -> str:
    """Tail of the most recently written logs/* file. Allrun runs with
    set -e, so the newest log belongs to the step that failed."""
    logs = [p for p in (case_dir / "logs").glob("*") if p.is_file()]
    if not logs:
        return "(no logs/ written)"
    newest = max(logs, key=lambda p: p.stat().st_mtime)
    return f"--- tail of logs/{newest.name} ---\n{_log_tail(newest)}"


def run_case(series_path: str, case_dir: Path, n_proc: int) -> None:
    """Mesh via Allrun (mesh-only) -> pv_meshView preview -> decompose ->
    solver+observer -> reconstruct -> pv_streamlines + pv_3slice_iso post.
    Per-case PNGs are case-prefixed and copied into progress/ so they end up
    in results.zip.

    Fails loud (RuntimeError) on a broken mesh or a solver that never wrote
    a time step -- callers mark the case as failed instead of reporting a
    green job with null results (2026-06-09 incident on A2_leitbleche)."""
    allrun = case_dir / "Allrun"
    try:
        os.chmod(allrun, 0o755)
    except OSError:
        pass
    rc = subprocess.call(["bash", str(allrun)], cwd=str(case_dir))
    if rc != 0 or not (case_dir / "constant" / "polyMesh" / "boundary").is_file():
        raise RuntimeError(
            "case %s FAILED at meshing (Allrun rc=%d, mesh %s)\n%s"
            % (case_dir.name, rc,
               "present" if (case_dir / "constant" / "polyMesh" / "boundary").is_file()
               else "missing",
               _newest_log_tail(case_dir))
        )

    # Pre-solver geometry preview (cheap, lets the user spot a broken
    # configuration before paying for the solver run).
    # The .foam touch is normally done by the solver wrapper later; for the
    # pre-solver pvbatch we touch it ourselves.
    touch(str(case_dir / f"{case_dir.name}.foam"))
    _run_pv_script(case_dir, "pv_meshView.py")
    _copy_pv_pngs(series_path, case_dir, ["pv_meshView.png"])

    run_OF_utility(str(case_dir), "decomposePar", params=["-force"])
    png_list = ["residuals.png", "solverIterations.png"]
    run_solver_copy_progress(series_path, png_list, str(case_dir), "simpleFoam",
                             log_file="logs/solver.log", n_proc=str(n_proc),
                             start_observer=True)
    if parse_iters(case_dir) is None:
        raise RuntimeError(
            "case %s FAILED: solver wrote no 'Time =' step (logs/solver.log)\n%s"
            % (case_dir.name, _newest_log_tail(case_dir))
        )
    run_OF_utility(str(case_dir), "reconstructPar", params=["-latestTime"])
    clean_processor_directories(str(case_dir))

    # Post-solver visualisation
    _run_pv_script(case_dir, "pv_streamlines.py")
    _run_pv_script(case_dir, "pv_3slice_iso.py")
    _copy_pv_pngs(series_path, case_dir,
                  ["pv_streamlines.png", "pv_3slice_iso.png"])

    _sweep_to_plot_pngs(series_path, case_dir)


def _run_pv_script(case_dir: Path, script_name: str) -> None:
    """Invoke a pvbatch script that lives in the case dir. Failures are
    non-fatal: a missing PNG just won't appear in the report."""
    script = case_dir / script_name
    if not script.is_file():
        print(f"WARNING: {script} missing, skipping pvbatch")
        return
    try:
        run_pvbatch_utility(str(case_dir),
                            params=[str(script), "--force-offscreen-rendering"])
    except Exception as e:
        print(f"WARNING: pvbatch {script_name} failed (non-fatal): {e}")


def _copy_pv_pngs(series_path: str, case_dir: Path,
                  filenames: list[str]) -> None:
    """Copy named PNGs from the case dir to progress/ with a case-prefix so
    they end up in results.zip."""
    progress = os.path.join(series_path, "progress")
    for fn in filenames:
        src = case_dir / fn
        if not src.is_file():
            continue
        try:
            shutil.copy(str(src), os.path.join(
                progress, f"{case_dir.name}_{fn}"))
        except Exception as e:
            print(f"WARNING: copy {src} -> progress failed: {e}")


def _sweep_to_plot_pngs(series_path: str, case_dir: Path) -> None:
    import glob
    to_plot = case_dir / "logs" / "to_plot"
    if not to_plot.is_dir():
        return
    progress = os.path.join(series_path, "progress")
    for p in glob.glob(str(to_plot / "*.png")):
        try:
            shutil.copy(p, os.path.join(
                progress, "%s_%s" % (case_dir.name, os.path.basename(p))))
        except Exception:
            pass


# ---------------------------------------------------------------------------
# result chart
# ---------------------------------------------------------------------------
def render_chart(png_path: str, case_stats: list[dict], *,
                 project_name: str, R_over_H: float, Re: float) -> str | None:
    """Two-panel chart: dP (top) + outlet uniformity (bottom) versus the
    sweep variable. The x-axis switches between number-of-vanes (Mode 1/2)
    and Brand-Pos (Mode 3, single-vane positions study) based on the cases.
    One line per N_inlet (mesh resolution)."""
    done = [c for c in case_stats if c.get("dP") is not None
            or c.get("uniformity") is not None]
    if not done:
        return None

    # If every case used Brand-Pos, label the x-axis as Brand-Pos; if mixed,
    # we still use the stored x_axis but label generically.
    all_brand = all(c.get("brand_pos") for c in case_stats) and case_stats
    if all_brand:
        x_label = "Brand-Pos N  (vane radius = (N+1) x 100 mm, H_phys = 4 m)"
        sub = "single-vane Pos-Sweep (Brand Fig. 7.14 / 7.15)"
    else:
        x_label = "number of guide vanes N"
        sub = "vane-count sweep"

    pbs = sorted({c["N_inlet"] for c in case_stats})
    fig, (ax_dp, ax_uni) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    cmap = plt.get_cmap("tab10")

    def _xval(c):
        return c.get("x_axis", c["N"])

    for idx, pb in enumerate(pbs):
        pts_dp = sorted([(_xval(c), c["dP"]) for c in done
                         if c["N_inlet"] == pb and c["dP"] is not None])
        pts_un = sorted([(_xval(c), c["uniformity"]) for c in done
                         if c["N_inlet"] == pb
                         and c["uniformity"] is not None])
        color = cmap(idx % 10)
        if pts_dp:
            ax_dp.plot([p[0] for p in pts_dp], [p[1] for p in pts_dp],
                       "o-", color=color, lw=2, ms=7,
                       label="N_inlet=%d" % pb)
        if pts_un:
            ax_uni.plot([p[0] for p in pts_un], [p[1] for p in pts_un],
                        "s-", color=color, lw=2, ms=7,
                        label="N_inlet=%d" % pb)

    ax_dp.set_ylabel("kinematic pressure loss dP [m^2/s^2]")
    ax_dp.set_title("%s -- A2 brand topology  (R/H=%.2g, Re=%.3g)\n%s"
                    % (project_name, R_over_H, Re, sub))
    ax_dp.grid(alpha=0.3)
    if len(pbs) > 1:
        ax_dp.legend(fontsize=9)

    ax_uni.set_ylabel("outlet uniformity gamma_U [-]")
    ax_uni.set_xlabel(x_label)
    ax_uni.set_ylim(0.0, 1.05)
    ax_uni.grid(alpha=0.3)
    if len(pbs) > 1:
        ax_uni.legend(fontsize=9)

    fig.tight_layout()
    try:
        fig.savefig(png_path, dpi=140)
    finally:
        plt.close(fig)
    return png_path


# ---------------------------------------------------------------------------
# results-zip (thin: charts + logs + postProcessing scalars; mesh/fields/VTK
# stay on the server in /srv/streamlit-portal/app/jobs/<job>/<case>/ for SFTP
# pickup). Size target: < 50 MB so app.py uses st.download_button and does
# not have to hit the broken nginx-signed-link path on a production HPC server.
# ---------------------------------------------------------------------------
def _write_sweep_summary(progress_path: str, *, project_name: str,
                         R_over_H: float, H: float, W: float,
                         L_in: float, L_out: float, L_trail: float,
                         Re: float, ti: float, n_layer_splits: int,
                         case_stats: list[dict]) -> None:
    """Write machine-readable summary of the sweep next to the chart PNG.

    Two files:
      summary.json  -- structured, sweep-level metadata + per-case results
      summary.csv   -- flat per-case rows with sweep metadata duplicated
                       (so two sweeps can be concatenated and grouped by
                       the metadata columns in a notebook/pandas).
    """
    import json
    import csv
    sweep_meta = {
        "workflow": "A2_brand_topology",
        "project_name": project_name,
        "R_over_H": R_over_H,
        "H": H, "W": W,
        "L_in": L_in, "L_out": L_out, "L_trail": L_trail,
        "Re": Re, "turbulent_intensity": ti,
        "n_layer_splits": n_layer_splits,
        "modified_vanes": L_trail > 0.0,  # Brand 7.3.4 convention
    }
    try:
        with open(os.path.join(progress_path, "summary.json"), "w",
                  encoding="utf-8") as f:
            json.dump({"sweep": sweep_meta, "cases": case_stats},
                      f, indent=2, allow_nan=True)
    except Exception as e:
        print("summary.json write failed (non-fatal): %s" % e)

    try:
        # Union of keys to keep CSV header stable across cases
        case_keys = []
        seen = set()
        for c in case_stats:
            for k in c.keys():
                if k not in seen:
                    seen.add(k); case_keys.append(k)
        meta_keys = list(sweep_meta.keys())
        with open(os.path.join(progress_path, "summary.csv"), "w",
                  encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(meta_keys + case_keys)
            for c in case_stats:
                row = [sweep_meta.get(k, "") for k in meta_keys] + \
                      [c.get(k, "") for k in case_keys]
                w.writerow(row)
    except Exception as e:
        print("summary.csv write failed (non-fatal): %s" % e)


def zip_results(series_path: str, results_path: str,
                case_dirs: list[Path]) -> None:
    """Bundle results.zip: progress/ + per-case mesh + reconstructed fields.

    The user needs the case to be openable in their own ParaView after
    download, so we include constant/ (polyMesh + dicts), 0/ (initial
    fields), and any reconstructed time directories (after reconstructPar
    -latestTime there should be just one). processor*/ has been cleaned
    up earlier by clean_processor_directories.
    """
    import zipfile
    series = Path(series_path)
    out_zip = Path(results_path) / "results.zip"
    INCLUDE_SUBDIRS = ("system", "constant", "logs", "postProcessing")
    EXTRA_FILES = ("Allrun", "Allclean")

    def _is_time_dir(d: Path) -> bool:
        if not d.is_dir():
            return False
        try:
            float(d.name)
            return True
        except ValueError:
            return False

    try:
        with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
            prog = series / "progress"
            if prog.is_dir():
                for f in prog.rglob("*"):
                    if f.is_file():
                        zf.write(f, f.relative_to(series))
            wlog = series / "worker.log"
            if wlog.is_file():
                zf.write(wlog, wlog.relative_to(series))
            for case_dir in case_dirs:
                if not case_dir.is_dir():
                    continue
                for sub in INCLUDE_SUBDIRS:
                    d = case_dir / sub
                    if not d.is_dir():
                        continue
                    for f in d.rglob("*"):
                        if f.is_file():
                            zf.write(f, f.relative_to(series))
                # Reconstructed time directories (0/ + latestTime/)
                for child in case_dir.iterdir():
                    if _is_time_dir(child):
                        for f in child.rglob("*"):
                            if f.is_file():
                                zf.write(f, f.relative_to(series))
                # Touch files that let "paraFoam -case ." find the .foam
                for fn in EXTRA_FILES:
                    f = case_dir / fn
                    if f.is_file():
                        zf.write(f, f.relative_to(series))
                for f in case_dir.glob("*.foam"):
                    if f.is_file():
                        zf.write(f, f.relative_to(series))
    except Exception as e:
        print("zip_results failed (non-fatal): %s" % e)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run_worker(series_path: str, *, do_source_of: bool = True,
               do_run: bool = True) -> None:
    series_path = os.path.abspath(series_path)
    while series_path.endswith(("/", "\\")):
        series_path = series_path[:-1]

    progress_path = os.path.join(series_path, "progress")
    results_path = os.path.join(series_path, "results")
    for p in (progress_path, results_path):
        os.makedirs(p, exist_ok=True)

    info: list[str] = ["py_A2_brand_topology.py is running now."]
    update_progress(progress_path, info)

    # ---- read + unroll the injected 6-tuple schema ----
    unrolled_name = json_to_unrolled(series_path, "interface.json")
    j = read_unrolled(series_path, unrolled_name)

    project_name = str(j.get("project name", os.path.basename(series_path)))
    R_over_H = float(j["bend radius ratio R/H"])
    H = float(j["channel height H [m]"])
    W = float(j["channel depth W [m]"])
    L_in = float(j["inlet length [H]"]) * H
    L_out = float(j["outlet length [H]"]) * H
    L_trail = float(j["trail length [H]"]) * H
    Re = float(j["Reynolds number"])
    ti = float(j["turbulent intensity"])

    vane_specs = parse_vane_specs(j["number of vanes"])
    Ninlet_serie = parse_int_series(j["cells across inlet N_inlet"],
                                    min_val=4, fallback=16)
    nz = max(2, int(float(j.get("cells in z (nz)", 8))))
    n_layer_splits = max(0, int(float(j.get("wall-layer splits", 0))))
    layer_thickness = float(j.get("layer-split thickness fraction", 0.5))
    n_proc = max(1, min(128, int(float(j.get("number of processors", 2)))))
    contact = str(j.get("contact email", "")).strip()

    R = R_over_H * H
    U = Re * NU / H

    from itertools import product
    combos = list(product(vane_specs, Ninlet_serie))   # [(spec, N_inlet), ...]

    def _fmt_spec(spec) -> str:
        if isinstance(spec, tuple) and len(spec) == 2 and spec[0] == "P":
            pos_val = spec[1]
            if isinstance(pos_val, (list, tuple)):
                return ("Brand-Pos [%s] (%d vanes)"
                        % (", ".join(str(int(p)) for p in pos_val),
                           len(pos_val)))
            mm = int(round((int(pos_val) + 1) * BRAND_POS_STEP_MM))
            return ("Brand-Pos %d  (1 vane at %d mm, fraction=%.5g)"
                    % (int(pos_val), mm,
                       brand_pos_to_fraction(int(pos_val))))
        if isinstance(spec, (list, tuple)):
            return ("explicit fractions [%s] (%d vanes)"
                    % (", ".join("%.5g" % f for f in spec), len(spec)))
        return "N=%d uniform (i/(N+1))" % int(spec)

    info += [
        "=" * 60,
        "A2 brand_topology -- parallel-translated 90deg bend",
        "=" * 60,
        "Project        : %s" % project_name,
        "Geometry       : H=%.4g m, W=%.4g m, R/H=%.4g (R=%.4g m)"
        % (H, W, R_over_H, R),
        "                 L_in=%.4g m (%.4g H), L_out=%.4g m (%.4g H), "
        "L_trail=%.4g m (%.4g H)"
        % (L_in, L_in / H, L_out, L_out / H, L_trail, L_trail / H),
        "Flow           : Re=%.4g, nu=%.4g m^2/s, U_inlet=%.4g m/s (=Re*nu/H), TI=%.3g"
        % (Re, NU, U, ti),
        "Vane spec(s)   : %d" % len(vane_specs),
    ]
    for s in vane_specs:
        info.append("                  - %s" % _fmt_spec(s))
    info += [
        "Mesh           : N_inlet = %s, nz=%d, refineWallLayer x %d (frac=%.2g)"
        % (Ninlet_serie, nz, n_layer_splits, layer_thickness),
        "Solver         : mpirun -np %d (decomposePar scotch)" % n_proc,
        "Total cases    : %d  (vane-specs x N_inlet)" % len(combos),
        "=" * 60,
    ]
    update_progress(progress_path, info)

    if do_source_of:
        # CFD-Lauf erst auf der VM verifiziert. Lazy resolution: only the
        # actual solver run needs OpenFOAM -- offline builds (Bayes-loop
        # pre-warming, smoke tests) don't.
        of_bashrc, of_source = py_OF_utils.resolve_openfoam_bashrc()
        info.append("OpenFOAM bashrc: %s  (%s)" % (of_bashrc, of_source))
        update_progress(progress_path, info)
        py_OF_utils.source_OF(source_version=of_bashrc)

    png_path = os.path.join(progress_path, "brandTopology_dP_uniformity.png")
    case_stats: list[dict] = []
    case_dirs: list[Path] = []

    for i, (spec, N_inlet) in enumerate(combos):
        if killed(series_path):
            info.append("\n>>> command_kill detected -- stopping cleanly.")
            update_progress(progress_path, info)
            touch(os.path.join(series_path, "command_finished"))
            return

        fractions = resolve_vane_spec(spec)
        N_eff = len(fractions)
        is_brand_pos = (isinstance(spec, tuple)
                        and len(spec) == 2 and spec[0] == "P")
        if is_brand_pos:
            mode = "P"
            pos_val = spec[1]
            if isinstance(pos_val, (list, tuple)):
                pos_tag = "_".join(str(int(p)) for p in pos_val)
            else:
                pos_tag = str(int(pos_val))
            case_dir = (Path(series_path)
                        / ("%d_P%s_N%d_ni%d" % (i, pos_tag, N_eff, N_inlet)))
        else:
            mode = "x" if isinstance(spec, (list, tuple)) else ""
            case_dir = (Path(series_path)
                        / ("%d_N%d%s_ni%d" % (i, N_eff, mode, N_inlet)))
        info.append("\n>>> Case %d / %d   %s   N_inlet=%d"
                    % (i + 1, len(combos), _fmt_spec(spec), N_inlet))
        update_progress(progress_path, info)

        binfo = build_case_dir(
            case_dir, H=H, W=W, R=R, L_in=L_in, L_out=L_out, L_trail=L_trail,
            fractions=fractions, N_inlet=N_inlet, nz=nz, Re=Re, ti=ti,
            n_proc=n_proc,
            n_layer_splits=n_layer_splits, layer_thickness=layer_thickness,
        )
        case_dirs.append(case_dir)
        info.append("    built: U_inlet=%.4g, k=%.4g, omega=%.4g, n_vanes=%d"
                    % (binfo["U_inlet"], binfo["k_init"], binfo["omega_init"],
                       binfo["n_vanes"]))
        update_progress(progress_path, info)

        dP = None
        vort = None
        uni = None
        iters = None
        failed_msg = None
        if do_run:
            # CFD-Lauf erst auf der VM verifiziert
            try:
                run_case(series_path, case_dir, n_proc=n_proc)
            except Exception as exc:
                # Mark-and-continue: one broken case must not kill the
                # remaining sweep cases. The job as a whole still fails
                # loud at the end of run_worker.
                failed_msg = str(exc)
                info.append("    >>> CASE FAILED: %s" % exc)
                update_progress(progress_path, info)
            dP = parse_dat_last(case_dir, "inletPressure", col=1)
            vort = parse_dat_last(case_dir, "outletVorticity", col=1)
            uni = parse_dat_last(case_dir, "outletUniformity", col=1)
            iters = parse_iters(case_dir)
            info.append("    result: dP=%s, vorticity_outlet=%s, "
                        "uniformity_outlet=%s, iters=%s"
                        % ("%.5g" % dP if dP is not None else "n/a",
                           "%.5g" % vort if vort is not None else "n/a",
                           "%.5g" % uni if uni is not None else "n/a",
                           iters if iters is not None else "n/a"))
        else:
            info.append("    (solver run skipped -- offline build only)")
        update_progress(progress_path, info)

        # x-axis value for the chart: Brand-Pos N for Mode 3, vane count N
        # for Modes 1/2 (Mode 2 puts the case at the multi-vane count).
        if is_brand_pos:
            pos_val = spec[1]
            x_axis_val = (float(sum(pos_val)) / len(pos_val)
                          if isinstance(pos_val, (list, tuple)) and pos_val
                          else float(pos_val))
        else:
            x_axis_val = float(N_eff)

        case_stats.append({"case_name": case_dir.name,
                           "N_inlet": N_inlet, "dP": dP,
                           "vorticity": vort, "uniformity": uni,
                           "iters": iters, "N": N_eff,
                           "x_axis": x_axis_val,
                           "brand_pos": is_brand_pos,
                           "custom": isinstance(spec, (list, tuple)),
                           "failed": failed_msg is not None})

        try:
            render_chart(png_path, case_stats,
                         project_name=project_name,
                         R_over_H=R_over_H, Re=Re)
        except Exception as e:
            print("chart render failed (non-fatal): %s" % e)

    info.append("\nAll %d case(s) processed." % len(combos))
    if any(c["dP"] is not None or c["uniformity"] is not None
           for c in case_stats):
        info.append("Summary:")
        for c in case_stats:
            tag = "*" if c.get("custom") else " "
            info.append("  %sN=%-2d  N_inlet=%-3d  dP=%s  vort=%s  "
                        "uni=%s  iters=%s"
                        % (tag, c["N"], c["N_inlet"],
                           "%.5g" % c["dP"] if c["dP"] is not None else "n/a",
                           "%.5g" % c["vorticity"] if c["vorticity"] is not None else "n/a",
                           "%.5g" % c["uniformity"] if c["uniformity"] is not None else "n/a",
                           c["iters"] if c["iters"] is not None else "n/a"))
        info.append("  ('*' marks cases with explicit non-equidistant fractions)")
    update_progress(progress_path, info)

    # Machine-readable summary so two sweeps (e.g. unmodified vs modified
    # vanes) can be merged into a common chart after download.
    _write_sweep_summary(progress_path, project_name=project_name,
                         R_over_H=R_over_H, H=H, W=W, L_in=L_in, L_out=L_out,
                         L_trail=L_trail, Re=Re, ti=ti,
                         n_layer_splits=n_layer_splits,
                         case_stats=case_stats)

    zip_results(series_path, results_path, case_dirs)

    # Fail loud AFTER summary + zip: partial results of the healthy cases
    # stay downloadable, but the job must not present itself as green
    # (2026-06-09 silent-failure incident on A2_leitbleche).
    failed_cases = [c["case_name"] for c in case_stats if c.get("failed")]
    if failed_cases:
        raise RuntimeError(
            "%d of %d case(s) FAILED: %s -- healthy cases are summarised in "
            "results.zip; per-case details in progress/1_info.txt and the "
            "case logs/"
            % (len(failed_cases), len(combos), ", ".join(failed_cases))
        )

    touch(os.path.join(series_path, "command_finished"))
    print("Finished! (contact=%s)" % (contact or "-"))


if __name__ == "__main__":
    series = check_args()
    # Env-var overrides for headless / dry-build use (tools/run_workflow.py
    # --dry-build, offline tests, Bayes-loop pre-warm).
    _env_source = os.environ.get("A2_DO_SOURCE_OF", "1") not in ("0", "false", "False")
    _env_run = os.environ.get("A2_DO_RUN", "1") not in ("0", "false", "False")
    try:
        run_worker(series, do_source_of=_env_source, do_run=_env_run)
    except Exception as exc:
        import traceback
        traceback.print_exc()
        try:
            prog = os.path.join(os.path.abspath(series), "progress")
            os.makedirs(prog, exist_ok=True)
            with open(os.path.join(prog, "9_error.txt"), "w",
                      encoding="utf-8") as f:
                f.write("py_A2_brand_topology.py failed:\n%s\n"
                        % traceback.format_exc())
        except Exception:
            pass
        touch(os.path.join(os.path.abspath(series), "command_finished"))
        sys.exit(1)
