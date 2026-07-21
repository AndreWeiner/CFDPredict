#!/usr/bin/env python3
"""Headless CLI + Python API for running a CFDPredict workflow.

This is the "Streamlit-less" entry point: it does exactly what the portal
does on Run-click, just without the web form. Same worker, same series_dir
layout, same kill/finished protocol.

Two surfaces:

  1. CLI (for human users + RunGui-next subprocess spawn):
       python tools/run_workflow.py <workflow_name> <settings.json> \\
                                    [--workdir <dir>] \\
                                    [--name <run_name>] \\
                                    [--no-tail]

  2. Python API (for Bayes loops, RunGui-next in-process):
       from tools.run_workflow import run_workflow_headless, find_workflow
       result = run_workflow_headless(
           "A2_brand_topology",
           settings,                       # dict, flat form values
           workdir=Path("./runs/x42"),
       )
       # result = {"exit_code": int, "series_dir": Path, "summary": dict|None}

What the CLI does, mirroring the Streamlit page:
  - Discover the workflow under workflows/<workflow_name>/.
  - Load its 6-tuple schema (interface.json), validate field names.
  - Inject the flat user settings into slot[1] of the schema.
  - Write the injected schema as <workdir>/interface.json.
  - Spawn the worker (workflow.yaml: entry_script + python_exec).
  - Tail progress/*.txt for live console output (Ctrl-C is honoured:
    touch command_kill, wait for command_finished, then return).
  - On finish, attempt to parse progress/summary.json and surface it.

Exit codes:
  0  worker finished cleanly (command_finished + no 9_error.txt)
  1  worker errored (command_finished present + 9_error.txt or non-zero exit)
  2  user-driven abort (Ctrl-C -> command_kill -> clean stop)
  3  pre-flight failure (workflow not found, schema mismatch, ...)
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
WORKFLOWS_DIR = REPO_ROOT / "workflows"


# ---------------------------------------------------------------------------
# workflow.yaml -- tiny flat parser (avoids hard PyYAML dep for the CLI)
# ---------------------------------------------------------------------------
def _parse_workflow_yaml(path: Path) -> dict[str, Any]:
    """Parse the small subset of YAML our workflow.yaml uses: flat
    `key: value` lines, comments (#), simple booleans / ints / strings.
    Lists and nested maps are not used in the manifest -- if you need
    them, switch to PyYAML.
    """
    out: dict[str, Any] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip() or ":" not in line:
            continue
        k, _, v = line.partition(":")
        key = k.strip()
        val = v.strip()
        if val.startswith(("'", '"')) and val.endswith(val[0]) and len(val) >= 2:
            val = val[1:-1]
        if val.lower() in ("true", "false"):
            out[key] = (val.lower() == "true")
        else:
            try:
                out[key] = int(val)
            except ValueError:
                out[key] = val
    return out


# ---------------------------------------------------------------------------
# Workflow discovery
# ---------------------------------------------------------------------------
def find_workflow(name: str) -> tuple[Path, dict[str, Any]]:
    """Return (workflow_dir, manifest_dict) for the named workflow.
    Raises FileNotFoundError if not found."""
    wdir = WORKFLOWS_DIR / name
    if not wdir.is_dir():
        raise FileNotFoundError(
            f"Workflow {name!r} not found under {WORKFLOWS_DIR}. "
            f"Available: {[p.name for p in WORKFLOWS_DIR.iterdir() if p.is_dir()]}"
        )
    yml = wdir / "workflow.yaml"
    if not yml.is_file():
        raise FileNotFoundError(f"{yml} missing")
    manifest = _parse_workflow_yaml(yml)
    manifest.setdefault("entry_script", f"py_{name}.py")
    manifest.setdefault("input_filename", "input.json")
    manifest.setdefault("schema_file", "interface.json")
    manifest.setdefault("inject_schema", False)
    return wdir, manifest


def list_workflows() -> list[str]:
    """Return the names of all workflows that have a workflow.yaml."""
    if not WORKFLOWS_DIR.is_dir():
        return []
    return sorted(
        p.name for p in WORKFLOWS_DIR.iterdir()
        if p.is_dir() and (p / "workflow.yaml").is_file()
    )


# ---------------------------------------------------------------------------
# Schema injection (mirror of the Streamlit page's inject_schema:true mode)
# ---------------------------------------------------------------------------
def _inject_into_schema(schema: dict, values: dict) -> dict:
    """Return a deep copy of `schema` (6-tuple form) with slot 1 (the
    default value) replaced by `values[field]` for every field the user
    supplied.

    Walks subdict-selectors (slot 4) recursively so flat values from a
    form that exposes nested branches (e.g. A3_nozzle's "Show details" =
    advanced sub-schema) get injected into the right branch. Unknown user
    fields are surfaced as an error so settings typos don't silently get
    ignored."""
    schema = copy.deepcopy(schema)

    def _index(s, out):
        """Recurse through subdict-selector branches; map field_name -> entry."""
        for name, entry in s.items():
            if not isinstance(entry, list) or len(entry) < 5:
                continue
            out[name] = entry
            selector = entry[4]
            if isinstance(selector, dict):
                for branch in selector.values():
                    if (isinstance(branch, list) and len(branch) >= 2
                            and isinstance(branch[1], dict)):
                        _index(branch[1], out)

    field_to_entry: dict = {}
    _index(schema, field_to_entry)

    unknown = [k for k in values.keys() if k not in field_to_entry
               and not k.startswith("_")]
    if unknown:
        raise ValueError(
            f"settings.json has fields not in the workflow's schema: {unknown}. "
            f"Known fields: {sorted(field_to_entry)}"
        )
    for field, value in values.items():
        if field.startswith("_"):
            continue  # metadata, e.g. "_meta": {...}
        entry = field_to_entry[field]
        if entry[2] == "separator":
            raise ValueError(
                f"{field!r} is a separator, not a settable field"
            )
        entry[1] = value
    return schema


# ---------------------------------------------------------------------------
# Worker subprocess management
# ---------------------------------------------------------------------------
def _resolve_python_exec(manifest_python: str) -> str:
    """Pick the Python interpreter for the worker. workflow.yaml may
    hard-code /usr/bin/python3 (which is right on the Linux portal VM
    but wrong on Windows / clusters). Strategy:
      1. If manifest_python is empty -> inherit (this process's exe).
      2. If the literal path exists -> use it.
      3. Otherwise fall back to sys.executable, but warn loudly so the
         user can fix workflow.yaml or set up a venv.
    """
    if not manifest_python:
        return sys.executable
    if Path(manifest_python).is_file():
        return manifest_python
    # PATH lookup as a softer fallback (e.g. manifest has just "python3")
    via_path = shutil.which(manifest_python)
    if via_path:
        return via_path
    print(
        f"[run_workflow] WARN: workflow.yaml python_exec={manifest_python!r} "
        f"is not a file and not on PATH. Falling back to current "
        f"interpreter: {sys.executable}",
        file=sys.stderr,
    )
    return sys.executable


class _Tailer:
    """Best-effort progress tailing: print new content of every
    progress/*.txt file since last poll. Survives files appearing
    mid-run. Not a perfect ordering across files, but adequate -- the
    worker tends to write one file at a time."""

    def __init__(self, progress_dir: Path):
        self.progress_dir = progress_dir
        self.offsets: dict[Path, int] = {}

    def poll(self) -> None:
        if not self.progress_dir.is_dir():
            return
        for f in sorted(self.progress_dir.glob("*.txt")):
            try:
                size = f.stat().st_size
            except OSError:
                continue
            prev = self.offsets.get(f, 0)
            if size <= prev:
                continue
            try:
                with open(f, "rb") as fh:
                    fh.seek(prev)
                    chunk = fh.read(size - prev)
            except OSError:
                continue
            self.offsets[f] = size
            text = chunk.decode("utf-8", errors="replace")
            if text:
                sys.stdout.write(text)
                sys.stdout.flush()


def _wait_for_finish(proc: subprocess.Popen, series_dir: Path,
                     poll_interval: float, tail: bool,
                     timeout_s: float | None) -> tuple[int, str | None]:
    """Block until the worker emits command_finished or the subprocess
    exits. Returns (exit_code, reason). Ctrl-C touches command_kill and
    waits a grace period (30 s) for the worker to clean up."""
    finished_marker = series_dir / "command_finished"
    kill_marker = series_dir / "command_kill"
    tailer = _Tailer(series_dir / "progress") if tail else None
    t0 = time.time()
    abort_requested = False

    def _request_abort(_sig=None, _frame=None):
        nonlocal abort_requested
        if not abort_requested:
            print("\n[run_workflow] Ctrl-C received -- writing command_kill, "
                  "waiting for worker to wind down...", file=sys.stderr)
            kill_marker.touch()
            abort_requested = True

    prev_int = signal.signal(signal.SIGINT, _request_abort)
    try:
        while True:
            if tailer:
                tailer.poll()
            if finished_marker.is_file():
                # let the worker exit its own way
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    pass
                reason = "command_kill" if abort_requested else "normal"
                return (proc.returncode if proc.returncode is not None else 0,
                        reason)
            ret = proc.poll()
            if ret is not None:
                # Worker died without writing command_finished -- treat as error
                if tailer:
                    tailer.poll()
                return ret, "worker_exit_no_marker"
            if timeout_s is not None and (time.time() - t0) > timeout_s:
                print("[run_workflow] timeout reached -- writing command_kill",
                      file=sys.stderr)
                kill_marker.touch()
                # give worker 30s to clean up
                deadline = time.time() + 30
                while time.time() < deadline and proc.poll() is None:
                    if tailer:
                        tailer.poll()
                    time.sleep(poll_interval)
                if proc.poll() is None:
                    proc.terminate()
                return proc.returncode if proc.returncode is not None else 124, "timeout"
            if abort_requested:
                deadline = time.time() + 30
                while time.time() < deadline and proc.poll() is None:
                    if tailer:
                        tailer.poll()
                    if finished_marker.is_file():
                        break
                    time.sleep(poll_interval)
                if proc.poll() is None:
                    proc.terminate()
                return proc.returncode if proc.returncode is not None else 130, "ctrl_c"
            time.sleep(poll_interval)
    finally:
        signal.signal(signal.SIGINT, prev_int)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def run_workflow_headless(
    workflow_name: str,
    settings: dict,
    workdir: Path,
    *,
    tail: bool = True,
    timeout_s: float | None = None,
    poll_interval: float = 0.5,
    extra_env: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Run a workflow without Streamlit. Returns:
        {
          "exit_code":  int,            # 0 ok, 1 worker error, 2 ctrl-c, 3 pre-flight
          "series_dir": Path,           # workdir, with all worker artefacts
          "summary":    dict | None,    # parsed progress/summary.json if present
          "reason":     str,            # "normal" | "ctrl_c" | "timeout" | ...
          "error":      str | None,     # short pre-flight or worker error text
        }
    """
    workdir = Path(workdir).resolve()
    try:
        wdir, manifest = find_workflow(workflow_name)
    except FileNotFoundError as e:
        return {"exit_code": 3, "series_dir": workdir, "summary": None,
                "reason": "workflow_not_found", "error": str(e)}

    schema_file = wdir / manifest["schema_file"]
    if not schema_file.is_file():
        return {"exit_code": 3, "series_dir": workdir, "summary": None,
                "reason": "schema_missing", "error": f"{schema_file} missing"}

    try:
        schema = json.loads(schema_file.read_text(encoding="utf-8"))
    except Exception as e:
        return {"exit_code": 3, "series_dir": workdir, "summary": None,
                "reason": "schema_unreadable", "error": str(e)}

    # Inject (or pass through flat) per inject_schema flag
    if manifest.get("inject_schema"):
        try:
            payload = _inject_into_schema(schema, settings)
        except ValueError as e:
            return {"exit_code": 3, "series_dir": workdir, "summary": None,
                    "reason": "schema_mismatch", "error": str(e)}
    else:
        payload = dict(settings)

    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / manifest["input_filename"]).write_text(
        json.dumps(payload, indent=2), encoding="utf-8")

    # Resolve python + entry script
    py = _resolve_python_exec(manifest.get("python_exec", ""))
    entry = wdir / manifest["entry_script"]
    if not entry.is_file():
        return {"exit_code": 3, "series_dir": workdir, "summary": None,
                "reason": "entry_missing", "error": f"{entry} missing"}

    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)

    print(f"[run_workflow] starting worker: {py} {entry} {workdir}", flush=True)
    proc = subprocess.Popen(
        [py, str(entry), str(workdir)],
        cwd=str(wdir),
        env=env,
    )

    code, reason = _wait_for_finish(proc, workdir, poll_interval, tail, timeout_s)

    # Parse summary.json (best effort) and surface error file if present
    summary = None
    summary_file = workdir / "progress" / "summary.json"
    if summary_file.is_file():
        try:
            summary = json.loads(summary_file.read_text(encoding="utf-8"))
        except Exception:
            pass

    err_file = workdir / "progress" / "9_error.txt"
    err_text = None
    if err_file.is_file():
        try:
            err_text = err_file.read_text(encoding="utf-8")[:2000]
        except Exception:
            pass

    if reason == "ctrl_c":
        exit_code = 2
    elif code != 0 or err_text:
        exit_code = 1
    else:
        exit_code = 0

    return {
        "exit_code": exit_code,
        "series_dir": workdir,
        "summary": summary,
        "reason": reason,
        "error": err_text,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _make_workdir(parent: Path | None, name: str | None,
                  workflow_name: str) -> Path:
    """Derive a workdir path. If --workdir is given, use it. If --name
    is given, use ./runs/<name>. Otherwise default to
    ./runs/<workflow_name>_<timestamp>."""
    if parent is not None:
        return Path(parent).resolve()
    base = Path.cwd() / "runs"
    if name:
        return (base / name).resolve()
    stamp = time.strftime("%Y%m%d_%H%M%S")
    return (base / f"{workflow_name}_{stamp}").resolve()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="run_workflow",
        description="Run a CFDPredict workflow headless from a settings.json.",
    )
    ap.add_argument("workflow_name", nargs="?",
                    help="Workflow under workflows/ (e.g. A2_brand_topology). "
                         "Use --list to see all.")
    ap.add_argument("settings_json", nargs="?",
                    help="Flat JSON of form values "
                         "(field-name -> value, as exported by the portal "
                         "'Save settings' button).")
    ap.add_argument("--workdir", type=Path, default=None,
                    help="Where the series_dir is created. Default: "
                         "./runs/<workflow>_<timestamp>")
    ap.add_argument("--name", default=None,
                    help="Short run name (used to form workdir if --workdir "
                         "is not set).")
    ap.add_argument("--no-tail", action="store_true",
                    help="Don't tail progress/*.txt to stdout while the "
                         "worker runs.")
    ap.add_argument("--timeout", type=float, default=None,
                    help="Hard timeout in seconds (worker is killed cleanly).")
    ap.add_argument("--dry-build", action="store_true",
                    help="Build case dirs but don't source OpenFOAM and don't "
                         "run the solver. Useful for offline smoke tests on "
                         "machines without OpenFOAM. Currently honoured by "
                         "A2_brand_topology (env vars A2_DO_SOURCE_OF / A2_DO_RUN).")
    ap.add_argument("--list", action="store_true",
                    help="List discoverable workflows and exit.")
    args = ap.parse_args(argv)

    if args.list:
        for name in list_workflows():
            print(name)
        return 0

    if not args.workflow_name:
        ap.error("workflow_name is required (unless --list)")
    if not args.settings_json:
        ap.error("settings_json is required (unless --list)")

    sp = Path(args.settings_json)
    if not sp.is_file():
        print(f"settings file not found: {sp}", file=sys.stderr)
        return 3
    try:
        settings = json.loads(sp.read_text(encoding="utf-8"))
    except Exception as e:
        print(f"settings file not parseable: {e}", file=sys.stderr)
        return 3

    workdir = _make_workdir(args.workdir, args.name, args.workflow_name)
    extra_env = None
    if args.dry_build:
        extra_env = {"A2_DO_SOURCE_OF": "0", "A2_DO_RUN": "0"}
    result = run_workflow_headless(
        args.workflow_name, settings, workdir,
        tail=(not args.no_tail), timeout_s=args.timeout,
        extra_env=extra_env,
    )

    print(f"\n[run_workflow] series_dir : {result['series_dir']}")
    print(f"[run_workflow] reason     : {result['reason']}")
    if result.get("summary"):
        s = result["summary"]
        # Workflows expose sweep metadata under "sweep" (A2_brand_topology
        # convention). Accept "sweep_meta" too for future workflows.
        meta = s.get("sweep") or s.get("sweep_meta") or {}
        cases = s.get("cases") or []
        print(f"[run_workflow] workflow   : {meta.get('workflow', '?')}")
        print(f"[run_workflow] cases done : {len(cases)}")
    if result.get("error"):
        print(f"[run_workflow] error tail :\n{result['error'][:1000]}",
              file=sys.stderr)

    return result["exit_code"]


if __name__ == "__main__":
    sys.exit(main())
