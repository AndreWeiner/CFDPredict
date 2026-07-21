# OUTPUTS -- what a workflow produces

A single workflow invocation produces a **series directory** with a
fixed layout, regardless of whether it was launched via the Streamlit
portal or `tools/run_workflow.py`:

```
<series_dir>/
├── interface.json                    # injected 6-tuple schema (input record)
├── interfaceUnrolled2.json           # unrolled flat key→value (worker-internal)
├── command_finished                  # zero-byte marker: worker is done
├── progress/
│   ├── 1_info.txt                    # live log (overwritten as worker progresses)
│   ├── summary.json                  # machine-readable results (Bayes consumes this)
│   ├── summary.csv                   # flat CSV (pandas-mergeable)
│   ├── 9_error.txt                   # only present if the worker crashed
│   ├── brandTopology_dP_uniformity.png   # (A2_brand_topology only) sweep chart
│   ├── <case>_pv_meshView.png        # (when pvbatch is available)
│   ├── <case>_pv_streamlines.png
│   └── <case>_pv_3slice_iso.png
├── results/
│   └── results.zip                   # bundled progress/ + each case dir (ParaView-openable)
└── <case_name>/                      # one or more, depending on workflow (see below)
    ├── system/                       # blockMeshDict, topoSetDict, createBafflesDict, controlDict, ...
    ├── constant/                     # turbulenceProperties, transportProperties
    ├── 0/                            # initial U, p, k, omega, nut
    ├── <latestTime>/                 # reconstructed solver state at convergence
    ├── postProcessing/               # surfaceFieldValue .dat output (see below)
    ├── logs/                         # blockMesh, solver, etc. (one file per OF utility)
    ├── *.foam                        # empty marker file: ParaView File > Open finds it
    └── Allrun, Allclean              # the build / cleanup scripts
```

## Number of `<case>/` directories per workflow

| Workflow | Cases per call | Why |
|---|---|---|
| `A2_leitbleche` | **1** | Single-case worker; series come from the Bayes-loop driver. |
| `A2_brand_topology` | **1..M** | Cartesian product of the vane-spec list (Mode-1/2/3) and the `cells across inlet` list. M ≤ 64 in practice. |

Case dir names encode the configuration:

- `A2_brand_topology`: `<i>_N<n_vanes>[_P<pos>][_x]_ni<N_inlet>` -- e.g. `0_N0_ni16`, `5_N5_ni16`, `0_P11_N1_ni8`, `0_N3x_ni20`.
  - `P<pos>` = Brand-Pos (Mode 3), `x` = custom-fractions case (Mode 2).
- `A2_leitbleche`: `case_N<n_vanes>_R<R>_Re<Re>` -- e.g. `case_N3_R1.50_Re1e+05`.

## `progress/summary.json` schema

```jsonc
{
  "sweep": {
    "workflow": "A2_leitbleche" | "A2_brand_topology",
    "project_name": "20260527_mother_default_3vane",

    // Geometry
    "H": 1.0, "W": 1.0, "R": 1.5,    // mother
    "L_in": 3.0, "L_out": 5.0,
    "L_trail": 0.0,                  // brand only (0 = unmodified, 0.125 = mod. vanes)
    "R_over_H": 0.75,                // brand only
    "n_vanes": 3,                    // mother only

    // Flow
    "Re": 1.0e5,
    "turbulent_intensity": 0.05,
    "nu": 1.5e-5,
    "U_inlet": 1.5,                  // mother (derived)

    // Mesh
    "ny": 40, "nz": 1,               // mother
    "nx_in": 60, "nx_bend": 80, "nx_out": 100,  // mother
    "n_layer_splits": 2,
    "auto_layers": true,             // mother only
    "layer_thickness": 0.5,

    // Brand-only modified-vanes flag
    "modified_vanes": false
  },

  "cases": [
    {
      "case_name": "case_N3_R1.50_Re1e+05",
      "N": 3, "R": 1.5, "H": 1.0,                          // mother
      "vane_radii": [1.7, 1.5, 1.3],                       // mother
      "vane_ext":   [0.0, 0.0, 0.0],                       // mother
      "N_inlet": 16,                                       // brand
      "brand_pos": false, "custom": false, "x_axis": 3.0,  // brand
      "Re": 1.0e5,
      "ny": 40, "nz": 1,
      "n_layer_splits": 2,
      "auto_layers": true,
      "y_plus_target": 40.0,                               // mother

      // Solver results (null if do_run=False or solver crashed)
      "dP":         0.43,        // m²/s² (kinematic pressure drop, areaAverage(p) inlet [- outlet])
      "vorticity":  6.22,        // 1/s (areaAverage(magVorticity) at outlet)
      "uniformity": 0.95,        // dimensionless (uniformity(U) at outlet, 1 = uniform)
      "iters":      1247         // SIMPLE iterations to convergence
    }
    // ...one entry per case
  ]
}
```

**Bayes-loop consumption** reads `summary["cases"][0]` (mother) or
iterates `summary["cases"]` (brand) and pulls the `dP` + `uniformity`
keys. See [`BAYES_INTEGRATION.md`](BAYES_INTEGRATION.md) for the
penalty convention when these are `null`.

## `progress/summary.csv` schema

Same data as `summary.json`, flattened to one row per case with the
`sweep`-level metadata duplicated across rows. Header columns are the
union of `sweep` keys + `cases[i]` keys, preserving first-seen order.

Two cross-workflow merges that just work:

```python
import pandas as pd
df1 = pd.read_csv("runs/brand_kap_7_3_3/progress/summary.csv")
df2 = pd.read_csv("runs/brand_kap_7_3_4/progress/summary.csv")
combined = pd.concat([df1, df2])
for L_trail, sub in combined.groupby("L_trail"):
    plt.plot(sub["N"], sub["dP"], "-o", label=f"L_trail={L_trail}")
```

```python
# Mother sweep (Bayes iterations or manual)
iters = pd.concat([
    pd.read_csv(f"runs/iter_{i:04d}/progress/summary.csv")
    for i in range(N_iter)
])
iters[["vane_radii", "vane_ext", "uniformity", "dP"]].plot.scatter("uniformity", "dP")
```

## Function-object `.dat` files (per case)

`postProcessing/<fo_name>/<startTime>/surfaceFieldValue.dat` -- ASCII,
space-separated, header starts with `#`. Columns:

| FO | Columns | Used by worker as |
|---|---|---|
| `inletPressure`  | `time   areaAverage(p)` | dP_in |
| `outletPressure` | `time   areaAverage(p)` | dP_out (mother only -- brand uses inletPressure direct) |
| `outletVorticity`| `time   areaAverage(magVorticity)` | vorticity |
| `outletUniformity`| `time  uniformity(U)` | uniformity (γ) |
| `yPlus`          | logged into solver.log, not .dat | diagnostic only |

Each FO writes every 50 timesteps (controlDict). The worker parses the
**last** numeric row, i.e. the converged value.

If you add a function object to `case_template/system/controlDict`,
remember to also teach the worker to parse it (`parse_dat_last(case_dir,
"<fo_name>", col=1)`) and add the result to the case-stats dict that
flows into `summary.json`.

## `9_error.txt`

Present iff the worker raised. Contains the Python traceback. The
caller's `result["error"]` from `run_workflow_headless` is the first
2000 chars of this file.

## Live-log `1_info.txt`

Plain text, **overwritten** (not appended) on every progress update.
Polling via file-mtime or whole-file re-read both work; the Streamlit
portal does the latter.

Format is human-readable, not machine-parseable -- use `summary.json`
for any automated consumer.
