# BAYES_INTEGRATION -- driving the workflow from a GP / EI / ParEGO loop

This document is the design guide for the TU-Dresden Bayes-loop on top
of the A2_leitbleche (mother) workflow. The shipped `examples/bayes_skeleton.py`
is a working botorch-based reference implementation; this file explains
the **decisions baked into it** so you (or your students) can adapt it
without re-deriving them.

If you haven't yet, read [`BAYES_TOPOLOGY_CHOICE.md`](BAYES_TOPOLOGY_CHOICE.md)
first -- it argues for the mother over the Brand topology and is shorter
than this file.

---

## The control flow

Every iteration of the Bayes loop is exactly this:

```
1. Acquisition picks  x = (vane_r1, ..., vane_r{N_max}, vane_ext1, ..., vane_ext{N_max})
2. settings  =  base_settings.copy()
   for i, r in enumerate(x[:N_max], 1):  settings[f"vane_r{i} [m]"] = r
   for i, e in enumerate(x[N_max:], 1):  settings[f"vane_ext{i} [m]"] = e

3. result = run_workflow_headless(
       "A2_leitbleche",  settings,
       workdir=Path(f"./runs/iter_{iter_id:04d}"),
       timeout_s=1800,
   )

4. if result["exit_code"] == 0  and  result["summary"] is not None:
      case = result["summary"]["cases"][0]
      y = (case["uniformity"], case["dP"])    # (maximize, minimize)
   else:
      y = penalty()                            # see below

5. observe (x, y); update GP; iterate.
```

That's the whole pattern. The workflow is a black box: input = flat
dict of field values, output = `(uniformity, dP)` or a penalty.

---

## Parameter space (mother, N_max=3)

| Component | Range | Notes |
|---|---|---|
| `vane_r1` | `(R-H/2+ε, R+H/2-ε)` = `(1.001, 1.999)` at R=1.5/H=1.0 | Outer vane. |
| `vane_r2` | `(1.001, vane_r1)` | Strictly inner of vane_r1 (validated by worker). |
| `vane_r3` | `(1.001, vane_r2)` | Strictly inner of vane_r2. |
| `vane_ext1` | `[0, L_out]` = `[0, 5]` | Downstream extension, parallel to outlet wall. |
| `vane_ext2` | `[0, L_out]` | Same. |
| `vane_ext3` | `[0, L_out]` | Same. |

**Continuous deactivation:** pushing `vane_r_i → R±H/2-ε` collapses the
corresponding radial band geometrically; setting `vane_ext_i=0` removes
the downstream extension. Both are continuous moves the GP can take --
no discrete topology jump. This is the central reason the mother
topology is GP-compatible where the Brand topology isn't.

**Ordering constraint** (`vane_r1 > vane_r2 > vane_r3`): two paths:
1. **Soft -- penalise violations** in the acquisition. The worker still
   sees a ValueError and writes `9_error.txt`, the loop records the
   penalty. Easy but wastes one Latin-hypercube draw per violation.
2. **Hard -- sample from an ordered region.** Re-parameterise via
   `(u₁, u₂, u₃)` with each in `(0, 1)` and map to
   `r_i = inner + (outer-inner) · sorted_descending(u)[i]`. The
   skeleton uses approach (1) for code simplicity; (2) is a cleaner
   GP-input topology if the acquisition's runtime budget can afford the
   reformulation.

**N_max > 3:** for the project's goal range (10-16 design params,
Dresden's parameter-budget recommendation), edit
`workflows/A2_leitbleche/interface.json` to add `vane_r4 [m]`, `vane_ext4 [m]`,
etc. The generator has no N_max ceiling.

---

## Objectives

Three function objects produce metrics in every solver case:

| Key in `summary["cases"][i]` | Direction | Range | What it means |
|---|---|---|---|
| `uniformity` | **maximize** | `[0, 1]`, 1 = uniform | `uniformity(U)` at outlet patch. Customer-billable metric for the propeller-anströmung use case (ZIM application target). |
| `vorticity` | minimize | `[0, ∞)` | `areaAverage(magVorticity)` at outlet. Secondary-flow diagnostic. |
| `dP` | minimize | `[0, ∞)` (typically) | Kinematic pressure drop (m²/s²). Bend loss; convertible to ζ-coefficient via `2·dP/U²`. |

**Recommended objective setup for ParEGO:**

- Primary: `uniformity ↑` (ZIM-customer-billable).
- Secondary: `dP ↓`.

`vorticity` is correlated with `uniformity` (Pearson r ≈ -0.83, R² ≈ 0.69
in our 2D-Iteration-7 sample) -- not independent enough to be a useful
third Pareto axis, but useful as a diagnostic / sanity-check column in
the iterations log.

If you want a third Pareto axis, consider a downstream-engineering
metric we haven't shipped yet (e.g. peak local velocity at outlet -- a
proxy for cavitation risk).

---

## Penalty handling (when a case fails)

A case can fail for several reasons, mapped to `result["reason"]` /
`result["exit_code"]` from `run_workflow_headless`:

| Failure mode | `exit_code` | `reason` | Suggested penalty |
|---|---|---|---|
| Solver crash mid-run (divergence, blockMesh error, ...) | `1` | `"worker_exit_no_marker"` or `"normal"` + `error` set | Worst observed objective + small δ (so the GP knows "bad here, don't revisit") |
| Schema-mismatch (e.g. unknown field) | `3` | `"schema_mismatch"` | Bug in the driver -- raise an exception, don't penalise. |
| Out-of-range radius / wrong ordering | `1` | `"normal"` + `error` references "must lie strictly inside" / "strictly decreasing" | Worst observed objective + δ (same as crash; the GP shouldn't propose this region again) |
| Timeout | `2` (Ctrl-C) or `1` + `reason=timeout` | -- | Same penalty + flag "needs more wallclock"; consider raising `timeout_s`. |

**Concrete penalty choice in the skeleton:**

```python
WORST_OBSERVED_UNIFORMITY = 0.7    # any sample worse than 0.7 in early sweeps
PENALTY_DP = 10.0                  # 10x the typical default dP ~ 0.4

def evaluate(x):
    result = run_workflow_headless(...)
    if result["exit_code"] == 0 and result["summary"]:
        case = result["summary"]["cases"][0]
        if case["uniformity"] is None or case["dP"] is None:
            return PENALTY_UNIFORMITY, PENALTY_DP    # solver wrote no FO output
        return case["uniformity"], case["dP"]
    return PENALTY_UNIFORMITY, PENALTY_DP            # any other failure
```

This biases the GP **away** from failing regions without poisoning the
posterior. For pathological cases (e.g. a whole region of the parameter
space crashes), inspect `progress/9_error.txt` of one failing case to
understand whether it's a workflow bug or a "don't propose here"
finding.

**Optional refinement:** treat the crash event as an additional
**feasibility classifier** input (constraint-aware BO, see Gardner et
al. 2014). Botorch supports this via `OutcomeTransform` /
`InfeasibilityTransform`. The skeleton does NOT do this -- penalty is
the simpler baseline.

---

## Parallelism

The mother worker runs **one case per call** -- internally
`mpirun -np N` for the OpenFOAM solver, but only one case at a time per
invocation. To exploit cluster parallelism:

### Pattern 1: q-batch acquisition (recommended)

Botorch's `qExpectedHypervolumeImprovement` (or `qExpectedImprovement`
for single-objective) proposes `q` candidate points per acquisition
step. Submit all `q` as parallel cluster jobs, gather results, update
GP, repeat. Wallclock per iteration ≈ slowest of the `q` jobs.

```python
# Pseudo
batch_x = q_acquisition.propose(q=4)
futures = [submit_to_slurm(x_i, workdir=f"runs/iter{i:04d}") for i, x_i in enumerate(batch_x)]
results = [f.result() for f in futures]      # blocks
gp.update(batch_x, results)
```

`tools/run_workflow.py` is process-safe: each call writes to its own
`workdir`, no global state. Submitting it as a SLURM job array is
straightforward.

### Pattern 2: Asynchronous BO

If wallclock varies a lot across the parameter space (it does for the
mother -- wall-resolved layers are 5-10x more expensive than
wall-function), async BO (Snoek et al. 2012, Kandasamy et al. 2018)
overlaps acquisition + evaluation better. Botorch supports this via
`pending_x` argument to the acquisition function.

This is a research direction; we don't ship a reference implementation.

### Pattern 3: One case per cluster node, many nodes

For initial design (e.g. 32-point Latin-hypercube before any GP
fitting), just SLURM-array-submit and forget. No coordination needed.

---

## Reproducibility / determinism

- **blockMesh + topoSet + createBaffles** are bit-deterministic.
- **refineWallLayer** is deterministic (single-pass edge-split).
- **simpleFoam** with default `linearUpwind` is deterministic given the
  same mesh + initial fields + nProcs **and** the same MPI
  decomposition seed.
- **decomposePar scotch** is deterministic given the same
  `numberOfSubdomains` and METIS/SCOTCH compiled with deterministic
  flags. If your cluster's MPI/scotch differs from DHCAE's, results
  will agree to ~5 significant figures but not bit-exact.

For ROM-training-data: pin `numberOfSubdomains` to a single value
across your entire training set (e.g. always `n_proc=4`). The
A2_leitbleche worker passes `n_proc` through from the form / settings
JSON -- don't randomize it.

For Bayes inner-loop: bit-exactness across cluster nodes is not
required. The objective noise is dominated by truncation error
(O(10⁻⁴) on `dP`, O(10⁻³) on `uniformity`), not MPI scatter.

---

## Wallclock-cost estimation

On `a production HPC server` (128-core production server, OpenFOAM v2512):

| Case | Mesh | y+ target | Wallclock (8 cores) |
|---|---|---|---|
| Empty bend (N=0), 2D, default ny=40 | ~22k cells | 40 (auto, 2 splits) | ~30 s |
| 3-vane default, 2D | ~30k cells | 40 (auto, 2 splits) | ~60 s |
| 3-vane default, 3D nz=20 | ~600k cells | 40 (auto, 2 splits) | ~3 min |
| 3-vane default, 3D nz=20 | ~600k cells | 1 (auto, 8 splits) | ~15 min |

A typical Bayes loop with **20 initial + 30 EI iterations = 50 cases**,
all 3D, y+~40, runs in ~3 hours of cluster wallclock on a single
8-core node (or ~10 minutes on a 100-node cluster with q=10 batches).

Budget for the full ZIM-AP2 training set is on the order of
**100-300 cases** -- comfortable on the Dresden cluster overnight.

---

## Recommended initial design

Latin-hypercube sampling over the 6-dim mother parameter space, projected
onto the ordering constraint (vane_r1 > vane_r2 > vane_r3):

```python
import scipy.stats.qmc as qmc
lhs = qmc.LatinHypercube(d=6, seed=0)
samples = lhs.random(n=20)                        # uniform in [0,1]^6
inner, outer = R - H/2 + eps, R + H/2 - eps
for s in samples:
    r_sorted = sorted([inner + s[i]*(outer-inner) for i in range(3)], reverse=True)
    ext = [s[3]*L_out, s[4]*L_out, s[5]*L_out]
    yield (*r_sorted, *ext)
```

20 points is on the small side for a 6-dim space; 30-40 is more
conservative. The skeleton defaults to `n_init = 20`.

---

## What to log per iteration

Append one row to `iterations.csv` per evaluation:

```
iter_id, vane_r1, vane_r2, vane_r3, vane_ext1, vane_ext2, vane_ext3, uniformity, dP, vorticity, iters, wallclock_s, exit_code, reason
```

The skeleton does this. After the loop, you have the full evaluation
history independent of the GP -- useful for sanity-plotting and for
re-running the GP fit with different hyperpriors.

Plus: every iteration's full `<series_dir>` is on disk under
`./runs/iter_XXXX/`, so the **raw CFD result** is recoverable for any
point in the Pareto front (or any point you want to re-postprocess).

---

## ROM-Trainingsdaten

Once the Bayes loop converges (or hits the budget), the same workflow
serves as a **data-generation engine** for the ROM:

1. Sample the parameter space (Latin-hypercube + Pareto-augmentation)
   to get a snapshot matrix.
2. For each sample, run the workflow → get the reconstructed solver
   state at `<latestTime>/`.
3. Read U, p, k, omega fields via PyFoam / paraview-python /
   foamFile-direct into a `(n_snapshots, n_cells, n_fields)` numpy
   tensor.
4. POD / kernel-POD / autoencoder onto that tensor (Weiner lab's
   choice).

This isn't shipped here -- the snapshot extraction script is a
Dresden-side artefact -- but the workflow's `<latestTime>/` output is
the input format you'd consume.

---

## Open questions for the next TUD-DHCAE meeting

These are listed in `BAYES_TOPOLOGY_CHOICE.md` already; reproduced here
in case you only have this file:

1. **N_max for the mother:** 3 (current MVP), 6 (Brand-asymptote), or 8
   (safety margin)? Cost trade-off: each extra vane = 2 extra Bayes
   variables.
2. **Vane-deactivation strategy:** rely on `vane_ext_i = 0` + radius-
   push (continuous), or add an explicit `active_i ∈ {0,1}` (mixed
   discrete/continuous GP)? Current recommendation: continuous-only.
3. **Pareto-axis count:** uniformity↑ + dP↓ (the skeleton's default),
   or add a third axis (peak-velocity, secondary-flow, ...)? The
   skeleton extends gracefully -- just add another column to the
   training data.
