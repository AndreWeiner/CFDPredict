#!/usr/bin/env python3
"""Report images: keep the slice / shell function-object output of the
steady cycle and render it with the case's pvbatch scripts.

prune_image_data() runs from the solver loop (solver.run_case) to bound
the disk use; render_images() runs once at the end (py_microchannel.py).
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

# postProcessing function-object names (case_template/system/controlDict)
FO_SLICES = "slices"              # surfaces: 4 z-planes every T/4
FO_SHELL = "shell"                # areaWrite: the shell (Ts_outerShell) every T/4

# Report images, rendered from the steady cycle's mid-forward (T/4) and
# mid-reverse (3T/4) output, one PNG per phase:
#   (function object, pvbatch script in the case, PNG name stem)
SLICE_PHASES = {"T4": 0.25, "3T4": 0.75}       # [-] of T
IMAGE_SETS = [(FO_SLICES, "pv_slices.py", "fluid_slices"),
              (FO_SHELL, "pv_shell.py", "shell_temperature")]


def _fo_dirs(case_dir: Path, fo: str) -> list[tuple[float, Path]]:
    """(time, dir) of every postProcessing/<fo>/<time>/, sorted by time."""
    base = case_dir / "postProcessing" / fo
    out = []
    if base.is_dir():
        for d in base.iterdir():
            try:
                out.append((float(d.name), d))
            except ValueError:
                continue
    return sorted(out)


def prune_image_data(case_dir: Path, period: float, t_now: float) -> None:
    """Keep only the IMAGE_SETS output of the cycle in progress and the one
    before it -- the steady cycle is always the last complete one, and the
    kill lands in the cycle after it. Called from the poll loop, so each
    function object never has more than ~2 cycles (8 time dirs) on disk."""
    t_keep = (int(t_now // period) - 1) * period - 1e-6 * period
    for fo, _, _ in IMAGE_SETS:
        for t, d in _fo_dirs(case_dir, fo):
            if t < t_keep:
                shutil.rmtree(d, ignore_errors=True)


def pick_phase_dirs(case_dir: Path, fo: str, period: float,
                    cycle: int | None) -> dict[str, Path]:
    """The <fo> time dir at each SLICE_PHASES phase of `cycle`; for a phase
    that cycle lacks (or cycle None: no complete cycle), the latest dir at
    that phase. Missing phases are left out."""
    dirs = _fo_dirs(case_dir, fo)
    tol = 1e-3 * period
    picked = {}
    for key, phase in SLICE_PHASES.items():
        at_phase = [(t, d) for t, d in dirs
                    if abs((t / period - phase) - round(t / period - phase))
                    * period < tol]
        if cycle is not None:
            exact = [d for t, d in at_phase
                     if abs(t - (cycle + phase) * period) < tol]
            if exact:
                picked[key] = exact[0]
                continue
        if at_phase:
            picked[key] = at_phase[-1][1]
    return picked


def dedupe_vtp(path: Path) -> int:
    """Drop repeated same-name DataArrays from a VTK XML file in place;
    returns how many were dropped. v2406's areaWrite writes every area
    field twice (its DimensionedField pass matches the GeometricField too,
    objectRegistry::names<T> is non-strict), and ParaView 5.12's XML reader
    segfaults on the duplicate. Both copies are identical."""
    text = path.read_text(encoding="utf-8")
    dropped = 0

    def dedupe_block(m):
        nonlocal dropped
        seen = set()

        def keep_first(a):
            nonlocal dropped
            if a.group(1) in seen:
                dropped += 1
                return ""
            seen.add(a.group(1))
            return a.group(0)
        return re.sub(r"\s*<DataArray\b[^>]*\bName='([^']*)'[^>]*>"
                      r".*?</DataArray>", keep_first, m.group(0), flags=re.S)

    # per block: a name may legitimately be both cell and point data
    text = re.sub(r"<(CellData|PointData|FieldData)\b.*?</\1>",
                  dedupe_block, text, flags=re.S)
    if dropped:
        path.write_text(text, encoding="utf-8")
    return dropped


def render_images(series_path: str, case_dir: Path, period: float,
                  cycle: int | None) -> list[str]:
    """For every IMAGE_SETS entry: pick the steady cycle's T/4 + 3T/4 dirs,
    delete every other time dir of that function object, and render one PNG
    per phase with its pvbatch script. Copies the PNGs to progress/ with a
    case prefix. Non-fatal: returns log lines."""
    lines = []
    progress = os.path.join(series_path, "progress")
    for fo, script_name, stem in IMAGE_SETS:
        picked = pick_phase_dirs(case_dir, fo, period, cycle)
        keep = set(picked.values())
        for _, d in _fo_dirs(case_dir, fo):
            if d not in keep:
                shutil.rmtree(d, ignore_errors=True)
        if not picked:
            lines.append("    %s: none written -- no images" % fo)
            continue
        args = []
        for key, d in picked.items():
            for vtp in d.glob("*.vtp"):
                dedupe_vtp(vtp)
            t = float(d.name)
            args += [f"--{key}", d.name]
            lines.append("    %s %s: t = %s s (cycle %d, phase %.3f T)"
                         % (fo, key, d.name, int(t // period),
                            t / period - int(t // period)))
        log = case_dir / f"log.{script_name[:-3]}"
        try:
            with open(log, "w", encoding="utf-8") as fh:
                rc = subprocess.run(["pvbatch", "--force-offscreen-rendering",
                                     str(case_dir / script_name)] + args,
                                    cwd=str(case_dir), stdout=fh,
                                    stderr=subprocess.STDOUT).returncode
        except OSError as e:
            lines.append("    %s failed (non-fatal): %s" % (script_name, e))
            continue
        if rc != 0:
            lines.append("    %s exited %d (non-fatal), see %s"
                         % (script_name, rc, log.name))
        # one PNG per phase (<stem>_<key>.png) or one for all (<stem>.png)
        for png in sorted(case_dir.glob(f"{stem}*.png")):
            shutil.copy2(png, os.path.join(progress,
                                           f"{case_dir.name}_{png.name}"))
            lines.append("    image: %s" % png.name)
    return lines
