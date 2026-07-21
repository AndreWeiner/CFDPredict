#!/usr/bin/env python3
"""A2 Leitbleche -- Mutter-Topologie (konzentrische Vanes).

Single-case worker. Reads the 6-tuple schema (interface.json) injected by
the portal / tools/run_workflow.py, builds an ElbowGeometry + FlowParams,
calls build_case.build(), runs simpleFoam (k-omega-SST) via the
case_template Allrun, then parses the surfaceFieldValue function objects
for dP_kin, outlet vorticity, outlet uniformity.

Designed to be driven both from:
  - the Streamlit portal (web form -> series_dir/interface.json -> spawn),
  - tools/run_workflow.py (CLI / RunGui-next subprocess),
  - a Bayes-loop driver calling run_workflow_headless() one case at a time
    (returns the summary dict -> objective values for the next GP step).

Worker contract:
    python py_A2_leitbleche.py <series_dir>
  - <series_dir>/interface.json: injected 6-tuple schema (inject_schema:true).
  - writes progress/<step>_info.txt for the live log.
  - touches <series_dir>/command_finished at exit.
  - honours <series_dir>/command_kill (Bayes timeouts / Ctrl-C).

Env-var overrides (consumed at __main__ time):
  A2_DO_SOURCE_OF  -- "0" to skip the bashrc-source (offline tests / Bayes pre-warm)
  A2_DO_RUN        -- "0" to build the case only, no solver
"""
from __future__ import annotations

import ast
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:      # package form (see __init__.py); flat script dir falls back below
    from . import foam_dictionary
    from . import py_OF_utils
    from .py_OF_utils import (touch, run_OF_utility, run_solver_copy_progress,
                             clean_processor_directories, get_latestTime,
                             run_pvbatch_utility)
except ImportError:
    import foam_dictionary
    import py_OF_utils
    from py_OF_utils import (touch, run_OF_utility, run_solver_copy_progress,
                             clean_processor_directories, get_latestTime,
                             run_pvbatch_utility)

try:      # package form (see __init__.py); flat script dir falls back below
    from .gen_blockmesh import ElbowGeometry
    from .build_case import build as build_case
    from .params import FlowParams
    from .wall_resolution import layer_splits_for_Re, predicted_yplus
except ImportError:
    from gen_blockmesh import ElbowGeometry
    from build_case import build as build_case
    from params import FlowParams
    from wall_resolution import layer_splits_for_Re, predicted_yplus


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
    with open(os.path.join(series_path, json_fname), "r", encoding="utf-8") as f:
        json_input = json.load(f)
    json_unrolled = foam_dictionary.unroll_dict(json_input)
    out_name = json_fname.replace(".json", "Unrolled2.json")
    with open(os.path.join(series_path, out_name), "w", encoding="utf-8") as f:
        f.write(json.dumps(json_unrolled))
    return out_name


def read_unrolled(series_path: str, fname: str) -> dict:
    with open(os.path.join(series_path, fname), "r", encoding="utf-8") as f:
        return json.load(f)


def parse_bool(raw, default: bool = False) -> bool:
    s = str(raw).strip().lower()
    if s in ("1", "true", "yes", "y", "on"):
        return True
    if s in ("0", "false", "no", "n", "off"):
        return False
    return default


# ---------------------------------------------------------------------------
# progress helpers
# ---------------------------------------------------------------------------
def update_progress(progress_path: str, lines: list[str]) -> None:
    os.makedirs(progress_path, exist_ok=True)
    with open(os.path.join(progress_path, "1_info.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def killed(series_path: str) -> bool:
    return os.path.isfile(os.path.join(series_path, "command_kill"))


# ---------------------------------------------------------------------------
# result parsing (mirrors A2_brand_topology)
# ---------------------------------------------------------------------------
def parse_dat_last(case_dir: Path, fo_name: str, col: int) -> float | None:
    """Return the last numeric value in column `col` (0-indexed) of the
    surfaceFieldValue .dat file for the given functionObject. None if the
    file is missing or unparseable -- common when the solver crashed before
    its first write or the FO isn't in controlDict."""
    base = case_dir / "postProcessing" / fo_name
    if not base.is_dir():
        return None
    # FO writes one subdir per startTime; we want the most recent
    times = sorted([p for p in base.iterdir() if p.is_dir()],
                   key=lambda p: float(p.name) if p.name else 0)
    if not times:
        return None
    dat = times[-1] / "surfaceFieldValue.dat"
    if not dat.is_file():
        return None
    last = None
    try:
        with open(dat, encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if len(parts) > col:
                    try:
                        last = float(parts[col])
                    except ValueError:
                        pass
    except OSError:
        return None
    return last


def parse_iters(case_dir: Path) -> int | None:
    """Pick the largest 'Time = ' integer in logs/solver.log -- the iteration
    count at the time the solver stopped writing (converged or stopped)."""
    log = case_dir / "logs" / "solver.log"
    if not log.is_file():
        return None
    n = None
    try:
        with open(log, encoding="utf-8", errors="replace") as f:
            for line in f:
                s = line.strip()
                if s.startswith("Time = "):
                    try:
                        n = int(s.split("=", 1)[1].strip())
                    except ValueError:
                        pass
    except OSError:
        return None
    return n


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


def _check_case_health(case_dir: Path, allrun_rc: int) -> None:
    """Fail loud instead of green-with-nulls: a failed Allrun, a missing
    mesh or a solver that never wrote a time step must fail the job.
    (2026-06-09 incident: blockMesh FATAL, but the job still reported
    'Finished!' with an all-null summary and an empty results.zip.)"""
    problems = []
    if allrun_rc != 0:
        problems.append(f"Allrun exited with code {allrun_rc}")
    if not (case_dir / "constant" / "polyMesh" / "boundary").is_file():
        problems.append("no mesh (constant/polyMesh/boundary missing)")
    if parse_iters(case_dir) is None:
        problems.append("solver wrote no 'Time =' step (logs/solver.log)")
    if problems:
        raise RuntimeError(
            "case %s FAILED: %s\n%s"
            % (case_dir.name, "; ".join(problems), _newest_log_tail(case_dir))
        )


# ---------------------------------------------------------------------------
# solver run via case_template Allrun
# ---------------------------------------------------------------------------
def _run_pv_script(case_dir: Path, script_name: str) -> None:
    """Invoke a pvbatch script that lives in the case dir. Failures are
    non-fatal -- missing pvbatch just means no PNG."""
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
    """Copy PNGs from the case dir to progress/ with a case-prefix so
    they end up in results.zip."""
    progress = os.path.join(series_path, "progress")
    for fn in filenames:
        src = case_dir / fn
        if not src.is_file():
            continue
        dst = os.path.join(progress, f"{case_dir.name}_{fn}")
        try:
            shutil.copy2(src, dst)
        except OSError as e:
            print(f"WARNING: copy {src} -> {dst} failed: {e}")


def run_case(series_path: str, case_dir: Path, n_proc: int) -> None:
    """Allrun does blockMesh + (topoSet + createBaffles when vanes) +
    refineWallLayer + potentialFoam + simpleFoam. nProcs is read from
    decomposeParDict, so we patch that file before invoking Allrun.

    Then runs pv_meshView (mesh structure), pv_streamlines (mid-slice flow),
    pv_3slice_iso (3D field overview) for the report PNGs."""
    decomp = case_dir / "system" / "decomposeParDict"
    if decomp.is_file():
        text = decomp.read_text(encoding="utf-8")
        # OpenFOAM dict: `numberOfSubdomains  N;` -- replace the int.
        import re
        new = re.sub(r"(numberOfSubdomains\s+)\d+\s*;",
                     r"\g<1>%d;" % n_proc, text)
        if new != text:
            decomp.write_text(new, encoding="utf-8", newline="\n")

    allrun = case_dir / "Allrun"
    if not allrun.is_file():
        raise FileNotFoundError(f"Allrun missing in {case_dir}")
    # Make sure it's executable; bash on POSIX requires +x.
    try:
        import stat as _stat
        st = allrun.stat()
        allrun.chmod(st.st_mode | _stat.S_IXUSR | _stat.S_IXGRP | _stat.S_IXOTH)
    except Exception:
        pass
    # Pipe Allrun stdout/stderr into logs/ (Allrun already does this per
    # utility, but the top-level echo + any unexpected stderr should land
    # somewhere predictable -- pick allrun.log).
    log_path = case_dir / "logs" / "allrun.log"
    log_path.parent.mkdir(exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as fh:
        rc = subprocess.run(["bash", str(allrun)],
                            cwd=str(case_dir),
                            stdout=fh, stderr=subprocess.STDOUT,
                            check=False).returncode
    _check_case_health(case_dir, rc)

    # *.foam marker for ParaView
    touch(str(case_dir / f"{case_dir.name}.foam"))

    # Post-solver visualisations (non-fatal if pvbatch unavailable)
    _run_pv_script(case_dir, "pv_meshView.py")
    _run_pv_script(case_dir, "pv_streamlines.py")
    _run_pv_script(case_dir, "pv_3slice_iso.py")
    _copy_pv_pngs(series_path, case_dir,
                  ["pv_meshView.png", "pv_streamlines.png",
                   "pv_3slice_iso.png"])


# ---------------------------------------------------------------------------
# zip results (mirrors A2_brand_topology layout so summary.csv files can be
# merged across both workflows in pandas)
# ---------------------------------------------------------------------------
INCLUDE_SUBDIRS = ("system", "constant", "logs", "postProcessing")


def _is_time_dir(name: str) -> bool:
    """OpenFOAM time directories are decimal numbers (0, 100, 500, 0.001, ...).
    Used to include 0/ + reconstructed latestTime/ in the zip."""
    try:
        float(name)
        return True
    except ValueError:
        return False


def zip_results(series_path: str, results_path: str,
                case_dirs: list[Path]) -> None:
    import zipfile
    out = os.path.join(results_path, "results.zip")
    try:
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
            # 1. progress/
            prog = os.path.join(series_path, "progress")
            if os.path.isdir(prog):
                for root, _, files in os.walk(prog):
                    for fn in files:
                        full = os.path.join(root, fn)
                        zf.write(full,
                                 os.path.relpath(full, series_path))
            # 2. per-case dirs
            for cd in case_dirs:
                if not cd.is_dir():
                    continue
                # *.foam marker
                for foam in cd.glob("*.foam"):
                    zf.write(foam, foam.relative_to(Path(series_path)))
                # standard subdirs
                for sub in cd.iterdir():
                    if sub.is_dir() and (sub.name in INCLUDE_SUBDIRS
                                         or _is_time_dir(sub.name)):
                        # skip processor*/ -- redundant with reconstructed time
                        if sub.name.startswith("processor"):
                            continue
                        for root, _, files in os.walk(sub):
                            for fn in files:
                                full = os.path.join(root, fn)
                                zf.write(full,
                                         os.path.relpath(full, series_path))
    except Exception as e:
        print("zip_results failed (non-fatal): %s" % e)


def _write_summary(progress_path: str, *, project_name: str,
                   geom_meta: dict, flow_meta: dict, mesh_meta: dict,
                   case_stats: list[dict]) -> None:
    """summary.json / summary.csv for Bayes-loop consumption + cross-sweep
    pandas merging. Schema matches A2_brand_topology where possible."""
    import csv
    sweep_meta = {
        "workflow": "A2_leitbleche",
        "project_name": project_name,
        **geom_meta,
        **flow_meta,
        **mesh_meta,
    }
    try:
        with open(os.path.join(progress_path, "summary.json"), "w",
                  encoding="utf-8") as f:
            json.dump({"sweep": sweep_meta, "cases": case_stats},
                      f, indent=2, allow_nan=True)
    except Exception as e:
        print("summary.json write failed (non-fatal): %s" % e)

    try:
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

    info: list[str] = ["py_A2_leitbleche.py is running now."]
    update_progress(progress_path, info)

    # ---- read + unroll the injected 6-tuple schema ----
    unrolled_name = json_to_unrolled(series_path, "interface.json")
    j = read_unrolled(series_path, unrolled_name)

    project_name = str(j.get("project name", os.path.basename(series_path)))

    # Geometry
    H = float(j["channel height H [m]"])
    W = float(j["channel depth W [m]"])
    R = float(j["bend center radius R [m]"])
    L_in = float(j["inlet length L_in [m]"])
    L_out = float(j["outlet length L_out [m]"])

    # Vanes
    n_vanes = max(0, min(3, int(float(j.get("number of vanes", 3)))))
    all_radii = [float(j[f"vane_r{i} [m]"]) for i in (1, 2, 3)]
    all_exts = [float(j[f"vane_ext{i} [m]"]) for i in (1, 2, 3)]
    vane_radii = all_radii[:n_vanes]
    vane_ext = all_exts[:n_vanes]

    # Validate geometry: R <= H/2 collapses the inner wall (R-H/2 <= 0,
    # blockMesh FATAL "inward-pointing faces") AND silently shifts the
    # valid vane band into the negative, so vane radii would pass the
    # range check below (2026-06-09 incident: H=4 with default R=1.5).
    if R <= H / 2.0:
        raise ValueError(
            f"bend center radius R={R:g} must exceed H/2={H / 2.0:g} "
            f"(inner wall radius R-H/2 would be <= 0). Validated regime "
            f"is R/H >= 0.6, e.g. R={0.6 * H:g} for H={H:g} (Brand 2020)."
        )

    # Validate radii
    inner, outer = R - H / 2.0, R + H / 2.0
    for r in vane_radii:
        if not (inner < r < outer):
            raise ValueError(
                f"vane radius {r} must lie strictly inside ({inner}, {outer}) "
                f"(R={R}, H={H}). Push it back into range or adjust R/H."
            )
    # Ordering r1 > r2 > r3 (the gen_blockmesh sorter handles arbitrary
    # order but we surface a clear error since the form labels imply outer..inner).
    for a, b in zip(vane_radii, vane_radii[1:]):
        if a <= b:
            raise ValueError(
                f"vane radii must be strictly decreasing (outer to inner): {vane_radii}"
            )
    # Ext range
    for e in vane_ext:
        if not (0.0 <= e <= L_out):
            raise ValueError(
                f"vane_ext {e} out of range [0, L_out={L_out}]"
            )

    # Flow
    Re = float(j["Reynolds number"])
    ti = float(j["turbulent intensity"])
    nu = 1.5e-5

    # Mesh
    nx_in = int(float(j.get("cells inlet axial nx_in", 60)))
    nx_bend = int(float(j.get("cells bend angular nx_bend", 80)))
    nx_out = int(float(j.get("cells outlet axial nx_out", 100)))
    ny = int(float(j.get("cells across channel ny", 40)))
    nz = max(1, int(float(j.get("cells in z (nz)", 20))))
    auto_layers = parse_bool(j.get("auto layers", "true"), default=True)
    y_plus_target = float(j.get("y+ target", 40.0))
    n_layer_splits_fixed = max(0, int(float(j.get("wall-layer splits", 0))))
    layer_thickness = float(j.get("layer-split thickness fraction", 0.5))

    n_proc = max(1, min(128, int(float(j.get("number of processors", 4)))))
    contact = str(j.get("contact email", "")).strip()

    # Compose geometry + flow params
    g = ElbowGeometry(
        H=H, W=W, R=R, L_in=L_in, L_out=L_out,
        nx_in=nx_in, nx_bend=nx_bend, nx_out=nx_out, ny=ny, nz=nz,
        vane_radii=vane_radii, vane_ext=vane_ext,
    )
    f = FlowParams(Re=Re, nu=nu, turbulent_intensity=ti)
    U_inlet, k_init, omega_init = f.inlet_state(H)

    # Layers
    if auto_layers:
        n_layer_splits = layer_splits_for_Re(Re, nu, H, y_plus_target, ny)
        yp_pred = predicted_yplus(Re, nu, H, ny, n_layer_splits) \
            if n_layer_splits > 0 else None
        layer_meta = (f"auto -> {n_layer_splits} splits "
                      f"(y+_target={y_plus_target}, y+_pred~"
                      f"{yp_pred:.1f})" if yp_pred is not None
                      else f"auto -> 0 splits (Re too low for refinement)")
    else:
        n_layer_splits = n_layer_splits_fixed
        layer_meta = f"fixed -> {n_layer_splits} splits"

    info += [
        "=" * 60,
        "A2 Leitbleche -- Mother-Topologie (konzentrische Vanes)",
        "=" * 60,
        "Project        : %s" % project_name,
        "Geometry       : H=%.4g, W=%.4g, R=%.4g  (R/H=%.4g)" % (H, W, R, R / H),
        "                 L_in=%.4g, L_out=%.4g" % (L_in, L_out),
        "Vanes          : %d active, r=%s, ext=%s"
        % (n_vanes, vane_radii, vane_ext),
        "Flow           : Re=%.4g, nu=%.4g, U_inlet=%.4g, TI=%.3g"
        % (Re, nu, U_inlet, ti),
        "                 k_init=%.4g, omega_init=%.4g" % (k_init, omega_init),
        "Mesh           : nx=(%d,%d,%d), ny=%d, nz=%d"
        % (nx_in, nx_bend, nx_out, ny, nz),
        "Wall layers    : %s, thickness frac %.2g"
        % (layer_meta, layer_thickness),
        "Solver         : mpirun -np %d (decomposePar scotch)" % n_proc,
        "=" * 60,
    ]
    update_progress(progress_path, info)

    if do_source_of:
        of_bashrc, of_source = py_OF_utils.resolve_openfoam_bashrc()
        info.append("OpenFOAM bashrc: %s  (%s)" % (of_bashrc, of_source))
        update_progress(progress_path, info)
        py_OF_utils.source_OF(source_version=of_bashrc)

    # ---- build the case ----
    if killed(series_path):
        info.append("\n>>> command_kill detected before build -- stopping.")
        update_progress(progress_path, info)
        touch(os.path.join(series_path, "command_finished"))
        return

    case_name = "case_N%d_R%.2f_Re%.0e" % (n_vanes, R, Re)
    case_dir = Path(series_path) / case_name
    info.append("\n>>> Building case: %s" % case_dir.name)
    update_progress(progress_path, info)

    binfo = build_case(case_dir, g, f,
                       n_layer_splits=n_layer_splits,
                       layer_thickness=layer_thickness)
    info.append("    built: U_inlet=%.4g, k=%.4g, omega=%.4g"
                % (binfo["U_inlet"], binfo["k_init"], binfo["omega_init"]))
    update_progress(progress_path, info)

    # marker file so ParaView -> File > Open finds the case
    (case_dir / "case.foam").touch()

    dP = vort = uni = None
    iters = None
    if do_run:
        if killed(series_path):
            info.append("\n>>> command_kill detected before solver -- stopping.")
            update_progress(progress_path, info)
            touch(os.path.join(series_path, "command_finished"))
            return
        info.append("    running Allrun (blockMesh, topoSet/createBaffles, "
                    "refineWallLayer, potentialFoam, simpleFoam) ...")
        update_progress(progress_path, info)
        try:
            run_case(series_path, case_dir, n_proc=n_proc)
        except Exception as exc:
            info.append("\n>>> CASE FAILED: %s" % exc)
            info.append(">>> Job is marked as failed (progress/9_error.txt).")
            update_progress(progress_path, info)
            raise

        # parse results
        dP_in = parse_dat_last(case_dir, "inletPressure", col=1)
        dP_out = parse_dat_last(case_dir, "outletPressure", col=1)
        dP = (dP_in - dP_out) if (dP_in is not None and dP_out is not None) \
            else dP_in
        vort = parse_dat_last(case_dir, "outletVorticity", col=1)
        uni = parse_dat_last(case_dir, "outletUniformity", col=1)
        iters = parse_iters(case_dir)
        info.append("    result: dP_kin=%s, vorticity_outlet=%s, "
                    "uniformity_outlet=%s, iters=%s"
                    % ("%.5g" % dP if dP is not None else "n/a",
                       "%.5g" % vort if vort is not None else "n/a",
                       "%.5g" % uni if uni is not None else "n/a",
                       iters if iters is not None else "n/a"))
    else:
        info.append("    (solver run skipped -- offline build only)")
    update_progress(progress_path, info)

    case_stats = [{
        "case_name": case_dir.name,
        "N": n_vanes,
        "R": R, "H": H,
        "vane_radii": vane_radii,
        "vane_ext": vane_ext,
        "Re": Re,
        "ny": ny, "nz": nz,
        "n_layer_splits": n_layer_splits,
        "auto_layers": auto_layers,
        "y_plus_target": y_plus_target if auto_layers else None,
        "dP": dP,
        "vorticity": vort,
        "uniformity": uni,
        "iters": iters,
    }]

    info.append("\nCase processed.")
    update_progress(progress_path, info)

    _write_summary(
        progress_path,
        project_name=project_name,
        geom_meta={
            "H": H, "W": W, "R": R,
            "L_in": L_in, "L_out": L_out,
            "n_vanes": n_vanes,
        },
        flow_meta={
            "Re": Re,
            "turbulent_intensity": ti,
            "nu": nu,
            "U_inlet": U_inlet,
        },
        mesh_meta={
            "ny": ny, "nz": nz,
            "nx_in": nx_in, "nx_bend": nx_bend, "nx_out": nx_out,
            "n_layer_splits": n_layer_splits,
            "auto_layers": auto_layers,
            "layer_thickness": layer_thickness,
        },
        case_stats=case_stats,
    )

    zip_results(series_path, results_path, [case_dir])

    touch(os.path.join(series_path, "command_finished"))
    print("Finished! (contact=%s)" % (contact or "-"))


if __name__ == "__main__":
    series = check_args()
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
                f.write("py_A2_leitbleche.py failed:\n%s\n"
                        % traceback.format_exc())
        except Exception:
            pass
        touch(os.path.join(os.path.abspath(series), "command_finished"))
        sys.exit(1)
