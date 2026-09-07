# CFDPredict A2 + A3 -- Standalone Workflow Package

DHCAE Tools GmbH · ZIM-Förderprojekt CFDPredict · 2026-05-27 (A3 added 2026-09)

This package contains the **CFDPredict A2 Leitbleche** and **A3 Pressure-Swirl
Nozzle** workflows in their headless form -- intended for the TU-Dresden
Bayes-optimization loop on the Uni-cluster, as the data-generation pipeline for
the ZIM-project's ROM-training-data work-package.

It's a curated subset of the `streamlit-demo-workflows` repository, with
the Streamlit web layer stripped and replaced by a Python CLI + API.
The Streamlit-Portal at DHCAE serves the same workflows in a different
shape; both produce identical `summary.json` / `summary.csv` outputs, so
results from the cluster and the portal are directly comparable and
mergeable.

## What's in this package

| Path | Content |
|---|---|
| [`workflows/A2_brand_topology/`](workflows/A2_brand_topology/) | Brand 2020-Kap.7 reproduction (parallel-translated vanes, variable count). |
| [`workflows/A2_leitbleche/`](workflows/A2_leitbleche/) | **A2 Mother topology** (concentric vanes + per-vane downstream extension). The Bayes-recommended variant -- see [`BAYES_TOPOLOGY_CHOICE.md`](BAYES_TOPOLOGY_CHOICE.md). |
| [`workflows/A3_nozzle/`](workflows/A3_nozzle/) | **A3 pressure-swirl nozzle** (Zhang et al. 2023). Parametric geometry (Di/Ds/Ls/alpha/Do/Lo/Lk/n_inlet) -> hybrid O-grid + snappyHexMesh -> interFoam VOF. Heavier than A2 -- see the dependency and compute-cost notes below. |
| [`tools/run_workflow.py`](tools/run_workflow.py) | CLI + Python API to drive any workflow without Streamlit. |
| [`examples/`](examples/) | Ready-to-run settings JSON files + a Bayes-loop skeleton (A2) + the full Zhang L16-DOE config set (A3, in `examples/A3_zhang_doe/`). |
| [`SCHEMA.md`](SCHEMA.md) | Field-by-field reference of every form input. |
| [`OUTPUTS.md`](OUTPUTS.md) | What `summary.json` / `results.zip` contain, units, columns. |
| [`BAYES_INTEGRATION.md`](BAYES_INTEGRATION.md) | How to drive the workflow from a GP / EI loop, including crash handling. |
| [`BAYES_TOPOLOGY_CHOICE.md`](BAYES_TOPOLOGY_CHOICE.md) | Why we recommend the mother over the Brand topology for the Bayes loop. |
| [`SYSTEM_REQUIREMENTS.md`](SYSTEM_REQUIREMENTS.md) | OS / OpenFOAM / Python / optional ParaView. |

## Quick start

### 1. System prerequisites

- Python ≥ 3.10
- OpenFOAM ESI v2406 / v2506 / v2512 (the worker is tested on v2512;
  Foundation 12 also works -- see [`SYSTEM_REQUIREMENTS.md`](SYSTEM_REQUIREMENTS.md)
  for the full resolution chain).
- MPI bundled with OpenFOAM (OpenMPI 4.x/5.x).
- Optional: ParaView with `pvbatch` on `PATH` for the `pv_*.py` PNG
  renders. Compute-only cluster nodes don't need it -- the worker
  swallows missing-pvbatch errors and continues.

### 2. Python dependencies

For **A2**, the runtime stack is intentionally tiny -- `matplotlib` is the
only non-stdlib requirement.

```bash
pip install -r workflows/A2_leitbleche/requirements.txt
# (same content as workflows/A2_brand_topology/requirements.txt)
```

**A3 is heavier**: it builds the nozzle geometry in Python (trimesh boolean
unions via the manifold3d engine) before snappy ever runs, so it needs a real
meshing-library stack, not just a plotting library:

```bash
pip install -r workflows/A3_nozzle/requirements.txt
# matplotlib, numpy, scipy, trimesh, manifold3d, mapbox-earcut, shapely, networkx
```

If this goes into a venv (recommended -- see the venv note in
[`SYSTEM_REQUIREMENTS.md`](SYSTEM_REQUIREMENTS.md)), `workflow.yaml`'s
`python_exec: ""` already inherits whichever interpreter launched
`run_workflow.py`, so activating the venv before running is enough -- no
config edit needed.

For the Bayes loop (additional, only needed by `examples/bayes_skeleton.py`):

```bash
pip install botorch gpytorch torch numpy
```

### 3. Smoke-test (no OpenFOAM needed)

```bash
python tools/run_workflow.py --list
# A2_brand_topology
# A2_leitbleche
# (plus older workflows shipped along for completeness)

python tools/run_workflow.py A2_leitbleche examples/mother_default_3vane.json \
                             --name smoke --dry-build --no-tail
# expected output:
#   [run_workflow] starting worker: ...
#   Finished! (contact=)
#   [run_workflow] cases done : 1

cd workflows/A2_leitbleche
python test_build_offline.py
# expected: "all offline smoke checks passed."
```

### 4. End-to-end (requires OpenFOAM)

```bash
# A2 mother with the default 3-vane configuration (~5 min on 4 cores)
python tools/run_workflow.py A2_leitbleche examples/mother_default_3vane.json \
                             --name mother_default --timeout 1800

# Result:
#   ./runs/mother_default/case_N3_R1.50_Re1e+05/   <- ParaView-openable case
#   ./runs/mother_default/progress/summary.json    <- dP, vorticity, uniformity
#   ./runs/mother_default/progress/summary.csv     <- same in CSV
#   ./runs/mother_default/results/results.zip      <- bundled case + progress
```

If `STREAMLIT_OPENFOAM_BASHRC` isn't set, the worker auto-resolves an
OpenFOAM installation -- see the resolution chain in
[`SYSTEM_REQUIREMENTS.md`](SYSTEM_REQUIREMENTS.md).

### 4b. A3 end-to-end -- and its compute cost

```bash
# Zhang's own predicted optimum, A4B1C2D4 (verification case)
python tools/run_workflow.py A3_nozzle examples/A3_zhang_doe/A4B1C2D4.json \
                             --name a4b1c2d4 --timeout 43200

# The full 16-case L16 DOE, sequentially, skipping finished cases on rerun:
bash examples/A3_zhang_doe/run_all.sh
```

⚠ **A3 is not a quick smoke case.** `workflow.yaml` targets 128 cores
(`max_concurrent_per_user: 1` -- these VOF runs are expensive), and the
DOE configs run `end time = 0.150 s` for spray-cone + inlet-pressure
convergence. Budget accordingly before queuing the full 16-case array;
for OpenFOAM-only smoke-testing, `--dry-build` (no solver) is available
the same way as for A2, and `--timeout` can be lowered to bail out of a
run early once a few time-steps have written (useful to confirm the
pipeline works before committing full cluster time).

If `STREAMLIT_OPENFOAM_BASHRC` isn't set for A3 either, only a single
cfdtools-specific fallback path is tried (unlike A2's longer resolution
chain) -- setting the env var explicitly is the safe default on any other
cluster.

### 5. Bayes loop (the actual project use case)

```bash
python examples/bayes_skeleton.py --n-init 6 --n-iter 20 \
                                  --workdir ./runs/bayes_demo
```

The skeleton:
- Generates a Latin-hypercube initial design over the 6-dim mother
  parameter space (`vane_r{1,2,3}` + `vane_ext{1,2,3}`).
- For each candidate point, calls `run_workflow_headless` from
  `tools/run_workflow.py`.
- Reads `progress/summary.json` to extract `(uniformity, dP)`.
- Fits a multi-output GP (botorch), runs EI / ParEGO acquisition, picks
  the next point.
- Logs to `iterations.csv` and a `pareto_front.png`.

See [`BAYES_INTEGRATION.md`](BAYES_INTEGRATION.md) for the design
rationale (penalty values, parallelism, crash handling, etc.).

## What's intentionally NOT in this package

- The Streamlit web app, user accounts, billing, mail dispatch. Those
  are the portal half; this package is the headless half.
- A Job-Queue or HPC scheduler integration. The expectation is that
  Dresden's existing cluster setup (SLURM / PBS / Snakemake / whatever)
  wraps `tools/run_workflow.py` invocations into job arrays. The
  worker is a regular Python subprocess -- nothing exotic.
- An auto-installer of OpenFOAM. The user installs OpenFOAM through
  whatever mechanism their cluster uses (module-system, container,
  manual build).
- A graphical results dashboard. ParaView opens the case directly;
  `summary.csv` opens in any spreadsheet.

## Reporting issues + asking questions

- Workflow bugs / unexpected solver crashes: file an issue against the
  upstream Git mirror (Dresden has access -- ask Martin if not).
- Topology / objective / parameter choice questions:
  [`BAYES_TOPOLOGY_CHOICE.md`](BAYES_TOPOLOGY_CHOICE.md) +
  [`BAYES_INTEGRATION.md`](BAYES_INTEGRATION.md) are the first stop;
  beyond that, the monthly TUD-DHCAE Teams meetings cover open
  decisions (Pareto-front shape, parameter budget, etc.).
- Output schema additions (new function objects, new metadata in
  `summary.json`): coordinate with Martin so the portal + cluster stay
  in sync.

## Version + provenance

This bundle is generated from the `streamlit-demo-workflows` repo at:
- branch / commit: see `BUNDLE_PROVENANCE.txt` (written by the build
  script at packaging time)
- GitHub mirror: `https://github.com/MartinBOF/streamlit-demo-workflows`
