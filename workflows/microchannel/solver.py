#!/usr/bin/env python3
"""Run the case's Allrun / Allrun-parallel and supervise it: per-cycle
reduction while it runs, the quasi-steady stop, command_kill, and the
fail-loud health check afterwards.
"""
from __future__ import annotations

import os
import re
import signal
import subprocess
import time
from glob import glob
from pathlib import Path

from cycle_stats import FO_DT, is_steady, steady_cycle_stats
from images import prune_image_data

KILL_POLL_S = 5.0


def killed(series_path: str) -> bool:
    return os.path.isfile(os.path.join(series_path, "command_kill"))


def _log_tail(path: Path, n_chars: int = 1500) -> str:
    try:
        size = path.stat().st_size
        with open(path, "rb") as f:
            f.seek(max(0, size - n_chars))
            return f.read().decode("utf-8", errors="replace")
    except OSError:
        return f"(log unreadable: {path})"


def _newest_log_tail(case_dir: Path) -> str:
    """Tail of the most recently written log.* file -- RunFunctions stop at
    the first failing step, so that is the one that failed."""
    logs = [p for p in case_dir.glob("log.*") if p.is_file()]
    if not logs:
        return "(no log.* written)"
    newest = max(logs, key=lambda p: p.stat().st_mtime)
    return f"--- tail of {newest.name} ---\n{_log_tail(newest)}"


def _solver_wrote_steps(case_dir: Path) -> bool:
    log = case_dir / "log.buoyantPimpleFoam"
    return log.is_file() and "\nTime = " in _log_tail(log, 200000)


def _check_case_health(case_dir: Path, allrun_rc: int) -> None:
    """Fail loud instead of green-with-nulls: a failed Allrun, a missing
    mesh or a solver that never wrote a time step must fail the job."""
    problems = []
    if allrun_rc != 0:
        problems.append(f"Allrun exited with code {allrun_rc}")
    if not (case_dir / "constant" / "polyMesh" / "boundary").is_file():
        problems.append("no mesh (constant/polyMesh/boundary missing)")
    if not _solver_wrote_steps(case_dir):
        problems.append("solver wrote no 'Time =' step (log.buoyantPimpleFoam)")
    if problems:
        raise RuntimeError(
            "case %s FAILED: %s\n%s"
            % (case_dir.name, "; ".join(problems), _newest_log_tail(case_dir))
        )


def _latest_time(case_dir: Path, name: str) -> float | None:
    """Last time in the newest postProcessing/<name>/*/*.dat, read from the
    file's tail -- a cheap "has another cycle finished?" probe."""
    files = sorted(glob(f"{case_dir}/postProcessing/{name}/*/*.dat"),
                   key=os.path.getmtime)
    if not files:
        return None
    for ln in reversed(_log_tail(Path(files[-1]), 4096).splitlines()):
        p = ln.split()
        if p and not ln.startswith("#"):
            try:
                return float(p[0])
            except ValueError:
                continue
    return None


def solver_end_time(case_dir: Path) -> float | None:
    """Last 'Time = ' the solver reached, from log.buoyantPimpleFoam."""
    log = case_dir / "log.buoyantPimpleFoam"
    if not log.is_file():
        return None
    hits = re.findall(r"^Time = ([0-9.eE+-]+)", _log_tail(log, 200000), re.M)
    return float(hits[-1]) if hits else None


def _stop(proc: subprocess.Popen) -> None:
    """SIGTERM the whole process group (bash, mpirun, every rank); SIGKILL
    if any of it is still alive after 30 s. Waits for the GROUP, not just
    Allrun's shell -- the ranks can outlive it by a moment."""
    def group_alive() -> bool:
        try:
            os.killpg(proc.pid, 0)
            return True
        except ProcessLookupError:
            return False

    os.killpg(proc.pid, signal.SIGTERM)
    deadline = time.time() + 30
    while time.time() < deadline:
        proc.poll()                     # reap Allrun's shell
        if not group_alive():
            break
        time.sleep(0.2)
    else:
        os.killpg(proc.pid, signal.SIGKILL)
    proc.wait()


def run_case(series_path: str, case_dir: Path, n_proc: int, period: float,
             on_cycle=None) -> str:
    """Allrun (serial) or Allrun-parallel: extrudeMesh x3, topoSet,
    createPatch, (decomposePar,) makeFaMesh, buoyantPimpleFoam.

    While it runs: each time another cycle completes, reduce the cycles so
    far (on_cycle(stats) gets them for the live log) and kill the run once
    is_steady(); kill it on command_kill too.

    Returns "endTime" (ran to the end), "steady" or "killed"."""
    script = case_dir / ("Allrun-parallel" if n_proc > 1 else "Allrun")
    if not script.is_file():
        raise FileNotFoundError(f"{script.name} missing in {case_dir}")
    n_done = 0
    with open(case_dir / "log.Allrun", "w", encoding="utf-8") as fh:
        proc = subprocess.Popen(["bash", str(script)], cwd=str(case_dir),
                                stdout=fh, stderr=subprocess.STDOUT,
                                start_new_session=True)
        while proc.poll() is None:
            if killed(series_path):
                _stop(proc)
                return "killed"
            t = _latest_time(case_dir, FO_DT)
            if t is not None:
                prune_image_data(case_dir, period, t)
            if t is not None and int(t // period) > n_done:
                n_done = int(t // period)
                stats = steady_cycle_stats(case_dir, period)
                if on_cycle is not None and stats["steady_cycle"] is not None:
                    on_cycle(stats)
                if is_steady(stats):
                    _stop(proc)
                    # the kill makes Allrun exit non-zero; only the mesh and
                    # the solver having stepped are checked
                    _check_case_health(case_dir, 0)
                    return "steady"
            time.sleep(KILL_POLL_S)
    _check_case_health(case_dir, proc.returncode)
    return "endTime"
