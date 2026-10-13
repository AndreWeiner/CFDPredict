#!/usr/bin/env python3
"""thermalShell micro-channel -- one channel of a thermal wheel.

Single-case worker. Reads the 6-tuple schema (interface.json) injected by
the portal / tools/run_workflow.py, builds the case with build_case.build()
(vertical-cut 2D mesh from mesher/, extruded in Z by the template's
Allrun), runs buoyantPimpleFoam + the finite-area thermal shell, then
reduces the delta_p / delta_T function objects over the last complete
oscillation cycle (the steady cycle) to the two outputs:

  dP_avg  mean of the forward and reverse |dp| plateaus            [Pa]
  eps     effectiveness, mean of the hot- and cold-side values      [-]
          = (dT_forward + dT_reverse) / (2 (T_hot - T_cold)), with
          dT_forward / dT_reverse the mean inlet-outlet dT over the
          forward (HOT inflow) / reverse (COLD inflow) half-pass

  
Designed to be driven both from:
  - the Streamlit portal (web form -> series_dir/interface.json -> spawn),
  - tools/run_workflow.py (CLI / RunGui-next subprocess),
  - a Bayes-loop driver calling run_workflow_headless() one case at a time
    (returns the summary dict -> objective values for the next GP step).

Worker contract:
    python py_microchannel.py <series_dir>
  - <series_dir>/interface.json: injected 6-tuple schema (inject_schema:true).
  - writes progress/<step>_info.txt for the live log.
  - touches <series_dir>/command_finished at exit.
  - honours <series_dir>/command_kill, also while the solver runs.

Outcomes, matching the A2 workers (docs/BAYES_INTEGRATION.md):
  exit 0 + summary.json   cases[0]["stop_reason"]: "steady" (change per
                          cycle < STEADY_TOL), "endTime" (not steady, see
                          the *_change fields) or None (dry build)
  exit 1, 9_error.txt     shape rejected by a gate, or mesh/solver failure
  exit 1, no summary      command_kill during the solve (timeout / Ctrl-C)
  exit 0, no summary      command_kill before the solver starts

Env-var overrides (consumed at __main__ time; the same names
tools/run_workflow.py --dry-build sets):
  A2_DO_SOURCE_OF  -- "0" to skip the bashrc-source (offline tests)
  A2_DO_RUN        -- "0" to build the case only, no solver
"""
from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

import foam_dictionary
import py_OF_utils
from py_OF_utils import touch
from build_case import build as build_case
from cycle_stats import STEADY_TOL, outputs, steady_cycle_stats
from images import render_images
from params import CaseParams
from solver import killed, run_case, solver_end_time
from true_geometry import check_parameter_values, check_true_geometry


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


# ---------------------------------------------------------------------------
# progress helpers
# ---------------------------------------------------------------------------
def update_progress(progress_path: str, lines: list[str]) -> None:
    os.makedirs(progress_path, exist_ok=True)
    with open(os.path.join(progress_path, "1_info.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def n_cells(case_dir: Path) -> int | None:
    """3D cell count from the polyMesh owner header."""
    owner = case_dir / "constant" / "polyMesh" / "owner"
    if not owner.is_file():
        return None
    with open(owner, "rb") as f:
        head = f.read(4096).decode("utf-8", "replace")
    m = re.search(r"nCells:\s*(\d+)", head)
    return int(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# zip results
# ---------------------------------------------------------------------------
INCLUDE_SUBDIRS = ("system", "constant", "postProcessing", "0")


def zip_results(series_path: str, results_path: str,
                case_dirs: list[Path]) -> None:
    import zipfile
    out = os.path.join(results_path, "results.zip")
    try:
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
            prog = os.path.join(series_path, "progress")
            if os.path.isdir(prog):
                for root, _, files in os.walk(prog):
                    for fn in files:
                        full = os.path.join(root, fn)
                        zf.write(full, os.path.relpath(full, series_path))
            for cd in case_dirs:
                if not cd.is_dir():
                    continue
                for f in (list(cd.glob("*.foam")) + list(cd.glob("log.*"))
                          + list(cd.glob("*.png"))):
                    zf.write(f, f.relative_to(Path(series_path)))
                for sub in cd.iterdir():
                    if sub.is_dir() and sub.name in INCLUDE_SUBDIRS:
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
    """summary.json / summary.csv, same layout as the A2 workflows."""
    import csv
    sweep_meta = {
        "workflow": "microchannel",
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


def _fmt(x, spec=".5g"):
    return format(x, spec) if x is not None else "n/a"


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

    info: list[str] = ["py_microchannel.py is running now."]
    update_progress(progress_path, info)

    # ---- read + unroll the injected 6-tuple schema ----
    unrolled_name = json_to_unrolled(series_path, "interface.json")
    j = read_unrolled(series_path, unrolled_name)

    p = CaseParams.from_interface(j, default_name=os.path.basename(series_path))

    # Gate 1 before anything divides by the inputs: N_ch divides by
    # b*(a+s) and the period by N_rot.
    bad = check_parameter_values(p.shape)
    if bad:
        info.append("\n>>> CASE REJECTED (invalid values): " + "; ".join(bad))
        update_progress(progress_path, info)
        raise ValueError("invalid parameter values:\n  " + "\n  ".join(bad))

    # Gate 2: values that are fine one by one may still form no channel
    # (fillets overlap, flank overhangs the trough, ...). Checked here, before
    # anything is built, and logged the same way as gate 1.
    bad = check_true_geometry(p.a, p.b, p.s, p.R)
    if bad:
        info.append("\n>>> CASE REJECTED (impossible geometry): " + "; ".join(bad))
        update_progress(progress_path, info)
        raise ValueError("no feasible fluid geometry:\n  " + "\n  ".join(bad))

    info += [
        "=" * 60,
        "thermalShell micro-channel",
        "=" * 60,
        "Project        : %s" % p.project_name,
        "Geometry [mm]  : a=%.4g, b=%.4g, R=%.4g, s=%.4g, l=%.4g"
        % (p.a, p.b, p.R, p.s, p.l),
        "Wheel          : N_rot=%.4g rpm -> period T=%.6g s, "
        "D_o=%.4g mm, D_i=%.4g mm"
        % (p.N_rot, p.period, p.D_outer, p.D_inner),
        "Mass flow      : stream %.6g kg/s over N_ch/2 of N_ch=%d channels -> "
        "mdotIn=%.6g kg/s (half channel)" % (p.mdot_stream, p.n_ch, p.mdot),
        "Temperatures   : hot in %.4g K, cold in %.4g K, shell init %.4g K"
        % (p.T_hot, p.T_cold, p.T_init),
        "Mesh           : refinement %.4g" % p.refinement,
        "Solver         : maxCo %.4g, nCorrectors %d, endTime %.4g s, "
        "%d rank(s)" % (p.max_co, p.n_correctors, p.end_time, p.n_proc),
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

    case_dir = Path(series_path) / p.case_name
    info.append("\n>>> Building case: %s" % case_dir.name)
    update_progress(progress_path, info)

    try:
        binfo = build_case(case_dir, p)
    except ValueError as exc:       # gates passed above: a mesher failure
        info.append("\n>>> BUILD FAILED: %s" % exc)
        update_progress(progress_path, info)
        raise
    # extrudeMesh stacks every 2D cell n_z times: total = cells_2d * n_z
    info.append("    built: %d cells" % (binfo["cells_2d"] * binfo["n_z"]))
    update_progress(progress_path, info)

    # marker file so ParaView -> File > Open finds the case
    (case_dir / "case.foam").touch()

    stats = steady_cycle_stats(case_dir, p.period)   # all None until solved
    stop_reason = t_stop = None
    if do_run:
        if killed(series_path):
            info.append("\n>>> command_kill detected before solver -- stopping.")
            update_progress(progress_path, info)
            touch(os.path.join(series_path, "command_finished"))
            return
        info.append("    running %s (extrudeMesh x3, topoSet, createPatch, "
                    "makeFaMesh, buoyantPimpleFoam) ..."
                    % ("Allrun-parallel" if p.n_proc > 1 else "Allrun"))
        update_progress(progress_path, info)
        def on_cycle(st):
            o = outputs(st, p.T_hot, p.T_cold)
            info.append("    cycle %d: dP_avg=%s Pa  eps=%s  "
                        "(change %s Pa / %s K / %s K)"
                        % (o["steady_cycle"], _fmt(o["dP_avg"]),
                           _fmt(o["eps"], ".5f"),
                           _fmt(o["dP_avg_change"], "+.4f"),
                           _fmt(o["dT_forward_change"], "+.4f"),
                           _fmt(o["dT_reverse_change"], "+.4f")))
            update_progress(progress_path, info)

        try:
            stop_reason = run_case(series_path, case_dir, n_proc=p.n_proc,
                                   period=p.period, on_cycle=on_cycle)
        except Exception as exc:
            info.append("\n>>> CASE FAILED: %s" % exc)
            info.append(">>> Job is marked as failed (progress/9_error.txt).")
            update_progress(progress_path, info)
            raise
        if stop_reason == "killed":
            # Same outcome as the A2 workers on a timeout: no summary.json,
            # non-zero exit (run_workflow -> exit_code 1, reason "timeout").
            # Unlike A2, the solver is already stopped (run_case kills the
            # whole process group).
            info.append("\n>>> command_kill detected during the run (t = %s s) "
                        "-- solver stopped, no summary written."
                        % _fmt(solver_end_time(case_dir), ".6g"))
            update_progress(progress_path, info)
            touch(os.path.join(series_path, "command_finished"))
            sys.exit(1)

        t_stop = solver_end_time(case_dir)
        info.append("    solver stopped at t = %s s (endTime %.4g s): %s"
                    % (_fmt(t_stop, ".6g"), p.end_time,
                       "quasi-steady, change < %g Pa / %g K / %g K -- killed"
                       % tuple(STEADY_TOL.values()) if stop_reason == "steady"
                       else "ran to endTime, not steady by the criterion"))
        stats = steady_cycle_stats(case_dir, p.period)
        info.append("    steady cycle k=%s of %d (t_end=%s s):"
                    % (stats["steady_cycle"], stats["n_cycles"],
                       _fmt(stats["t_end"], ".4g")))
        o = outputs(stats, p.T_hot, p.T_cold)
        info.append("      dP_avg=%s Pa  eps=%s"
                    % (_fmt(o["dP_avg"]), _fmt(o["eps"], ".5f")))
        info.append("      change vs previous cycle: dP_avg %s Pa, "
                    "dT_forward %s K, dT_reverse %s K"
                    % (_fmt(o["dP_avg_change"], "+.4f"),
                       _fmt(o["dT_forward_change"], "+.4f"),
                       _fmt(o["dT_reverse_change"], "+.4f")))
        info += render_images(series_path, case_dir, p.period,
                              stats["steady_cycle"])
    else:
        info.append("    (solver run skipped -- offline build only)")
    update_progress(progress_path, info)

    case_stats = [{
        "case_name": case_dir.name,
        **p.shape,
        "N_ch": p.n_ch,
        "mdot": p.mdot,
        "cells_3d": n_cells(case_dir),
        "period": p.period,
        "stop_reason": stop_reason,
        "t_stop": t_stop,
        **outputs(stats, p.T_hot, p.T_cold),
    }]

    info.append("\nCase processed.")
    update_progress(progress_path, info)

    _write_summary(
        progress_path,
        project_name=p.project_name,
        geom_meta={k: getattr(p, k) for k in ("a", "b", "R", "s", "l")},
        flow_meta={
            "N_rot": p.N_rot,
            "D_outer": p.D_outer,
            "D_inner": p.D_inner,
            "mdot_stream": p.mdot_stream,
            "T_hot": p.T_hot, "T_cold": p.T_cold, "T_init_shell": p.T_init,
        },
        mesh_meta={
            "refinement": p.refinement,
            "maxCo": p.max_co,
            "nCorrectors": p.n_correctors,
            "end_time": p.end_time,
            "steady_tol_dP_Pa": STEADY_TOL["dP_avg_change"],
            "steady_tol_dT_K": STEADY_TOL["dT_forward_change"],
            "n_proc": p.n_proc,
        },
        case_stats=case_stats,
    )

    zip_results(series_path, results_path, [case_dir])

    touch(os.path.join(series_path, "command_finished"))
    print("Finished! (contact=%s)" % (p.contact or "-"))


if __name__ == "__main__":
    series_dir = check_args()
    _env_source = os.environ.get("A2_DO_SOURCE_OF", "1") not in ("0", "false", "False")
    _env_run = os.environ.get("A2_DO_RUN", "1") not in ("0", "false", "False")
    try:
        run_worker(series_dir, do_source_of=_env_source, do_run=_env_run)
    except Exception as exc:
        import traceback
        traceback.print_exc()
        try:
            prog = os.path.join(os.path.abspath(series_dir), "progress")
            os.makedirs(prog, exist_ok=True)
            with open(os.path.join(prog, "9_error.txt"), "w",
                      encoding="utf-8") as f:
                f.write("py_microchannel.py failed:\n%s\n"
                        % traceback.format_exc())
        except Exception:
            pass
        touch(os.path.join(os.path.abspath(series_dir), "command_finished"))
        sys.exit(1)
