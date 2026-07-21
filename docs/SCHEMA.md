# SCHEMA -- A2 workflow form fields

Reference for the **flat settings JSON** consumed by
`tools/run_workflow.py`. The Streamlit portal renders the same fields
from the underlying `workflows/<name>/interface.json` schema (a 6-tuple
form). When the worker runs, the flat user values are injected into
slot 1 of each 6-tuple entry before the worker is spawned -- this
happens transparently in `tools/run_workflow.py`.

> If the key name and type matches the table below, your settings.json
> is valid. Unknown keys cause a `schema_mismatch` error from the CLI
> (with the list of valid field names attached). Type-coercion is
> permissive: "1e5" or 100000 both work for a `float` field.

The 6-tuple definition is at `workflows/<name>/interface.json`.

---

## Conventions across both workflows

- `_meta` (object) at the top of a settings.json is **ignored** by the
  worker -- safe place for human-readable scenario notes, dates,
  reference figure IDs.
- Field names use spaces and physical units (e.g. `"channel height H [m]"`)
  for portal-form readability. Quote them in Python:
  ```python
  settings["channel height H [m]"] = 1.0
  ```
- All length parameters are in metres, kinematic pressure in m²/s², no
  CGS or imperial anywhere. Re is the only non-dimensional input;
  `nu = 1.5e-5 m²/s` (air at 20 °C) is hardcoded in the workers.
- Both workflows accept `"true"` / `"false"` strings as well as Python
  bools for boolean-typed fields (`auto layers` etc.).

---

## A2_brand_topology

Parallel-translated 90° bend with full-height guide vanes. Reproduces
Brand 2020 Kap. 7.3 exactly (R/H=0.75, equidistant vanes at fractional
pathline positions, optional trail-block for modified vanes).

### Project metadata

| Field | Type | Default | Meaning |
|---|---|---|---|
| `project name` | string | `20260524_brandTopology_demo` | Case folder name. |

### Geometry

| Field | Type | Default | Range / notes |
|---|---|---|---|
| `bend radius ratio R/H` | float | `0.75` | 0.75 reproduces Brand 2020 (S. 114: bend-curvature radius = 0.25 H at H_phys=4 m). 1.5 / 2.0 / 3.0 are softer bends. |
| `channel height H [m]` | float | `1.0` | Reference length; Re = U·H/ν. |
| `channel depth W [m]` | float | `1.0` | 3D extrusion. The generator requires nz ≥ 2 (no 2D-empty patch yet). |
| `inlet length [H]` | float | `3.0` | Multiples of H, not metres. |
| `outlet length [H]` | float | `5.0` | Total outlet (includes trail block if any). |
| `trail length [H]` | float | `0.0` | **Brand-7.3.4 control.** 0 = unmodified, 0.125 = 500 mm at H_phys=4 m (Brand "modifizierte Leitbleche"). 0 ≤ trail ≤ outlet. |

### Flow

| Field | Type | Default | Range / notes |
|---|---|---|---|
| `Reynolds number` | float | `1e5` | U_inlet derived as Re·ν/H. |
| `turbulent intensity` | float | `0.05` | Inlet k-ω-SST BC: k = 1.5·(TI·U)². |

### Guide vanes (Brand-specific 3-mode input)

| Field | Type | Modes | Notes |
|---|---|---|---|
| `number of vanes` | string | **Mode 1**: int / int-list (equidistant). <br> **Mode 2**: float-list in (0,1) (explicit pathline fractions, single case). <br> **Mode 3**: prefix `P` for Brand-Pos (1 vane per case at Pos N = (N+1)·100 mm). | See examples below. |

```
"number of vanes": "3"               # Mode 1: one case with 3 equidist vanes
"number of vanes": "[0, 3, 5, 7]"    # Mode 1 sweep: 4 cases (Brand Fig. 7.17)
"number of vanes": "[0.05, 0.15, 0.30]"  # Mode 2: one case with 3 explicit fractions
"number of vanes": "P11"             # Mode 3: 1 case, one vane at Pos 11 (1200 mm)
"number of vanes": "P[1, 5, 11, 19, 29]"  # Mode 3 sweep: 5 single-vane cases (Brand Fig. 7.14/15)
```

### Mesh & wall treatment

| Field | Type | Default | Notes |
|---|---|---|---|
| `cells across inlet N_inlet` | string | `16` | int = 1 case; int-list = grid-convergence sweep. |
| `cells in z (nz)` | integer | `8` | Must be ≥ 2 (no 2D-empty patch). Keep small for demos. |
| `wall-layer splits` | integer | `0` | `refineWallLayer` invocations at walls + vanes (each halves the wall-adjacent cell). For y+ ~40 at Re=1e5: 2 splits. |
| `layer-split thickness fraction` | float | `0.5` | Fraction of the wall-adjacent cell at each split. |

### Run-time

| Field | Type | Default | Notes |
|---|---|---|---|
| `number of processors` | integer | `2` | MPI ranks (decomposePar scotch + `mpirun -np N`). Capped at 128 internally. |
| `contact email` | string | (auto-prefilled in portal) | Job-done mail target. Ignored by CLI. |
| `comments` | text_area | empty | Free-text notes. |

---

## A2_leitbleche (Mother)

Concentric guide vanes around a common bend center; per-vane downstream
extension via `boxToFace` topoSet. Bayes-/ROM-recommended topology.

### Project metadata

| Field | Type | Default | Meaning |
|---|---|---|---|
| `project name` | string | `20260527_A2_mother_demo` | Case folder name. |

### Geometry

| Field | Type | Default | Notes |
|---|---|---|---|
| `channel height H [m]` | float | `1.0` | Reference; Re = U·H/ν. |
| `channel depth W [m]` | float | `1.0` | nz=1 + W small → quasi-2D with empty front/back; nz>1 + W=H → full 3D. |
| `bend center radius R [m]` | float | `1.5` | Strict > H/2 (no inner-wall singularity). R/H ∈ [0.75, 3] sensible. |
| `inlet length L_in [m]` | float | `3.0` | 3·H sufficient for fully-developed k-ω-SST inflow. |
| `outlet length L_out [m]` | float | `5.0` | Must be ≥ max(vane_ext_i). |

### Guide vanes (MVP: N_max = 3)

| Field | Type | Default | Range / notes |
|---|---|---|---|
| `number of vanes` | integer | `3` | 0..3 active. Uses the first N entries of vane_r{i}/vane_ext{i}. 0 = empty bend (baseline). |
| `vane_r1 [m]` | float | `1.7` | Outer vane radius. **Must** be in (R-H/2, R+H/2). |
| `vane_r2 [m]` | float | `1.5` | Middle vane. **Must** be strictly less than `vane_r1`. |
| `vane_r3 [m]` | float | `1.3` | Inner vane. **Must** be strictly less than `vane_r2`. |
| `vane_ext1 [m]` | float | `0.0` | Downstream extension of vane 1 into the outlet (parallel to outlet wall). 0 = bend-only. ≤ L_out. |
| `vane_ext2 [m]` | float | `0.0` | Downstream extension of vane 2. |
| `vane_ext3 [m]` | float | `0.0` | Downstream extension of vane 3. |

Worker validates: out-of-range radius and wrong ordering raise
`schema_mismatch`. **Bayes-loop guidance:** to deactivate a vane
continuously without a topology jump, push its radius towards a wall
(r_i → R ± H/2 − ε) and/or set `vane_ext_i = 0`. See
[`BAYES_TOPOLOGY_CHOICE.md`](BAYES_TOPOLOGY_CHOICE.md) for the rationale.

### Flow

| Field | Type | Default | Notes |
|---|---|---|---|
| `Reynolds number` | float | `1e5` | U_inlet derived as Re·ν/H. ν = 1.5·10⁻⁵ m²/s hardcoded. |
| `turbulent intensity` | float | `0.05` | k = 1.5·(TI·U)². |

### Mesh & wall treatment

| Field | Type | Default | Notes |
|---|---|---|---|
| `cells inlet axial nx_in` | integer | `60` | Axial cells in the inlet straight. |
| `cells bend angular nx_bend` | integer | `80` | Angular cells around the 90° bend. |
| `cells outlet axial nx_out` | integer | `100` | Axial cells in the outlet. |
| `cells across channel ny` | integer | `40` | Cross-channel cells (sub-divided into N+1 radial bands). |
| `cells in z (nz)` | integer | `20` | `1` = quasi-2D + empty front/back; `>1` = full 3D, front/back merged into walls. |
| `auto layers` | string ("true"/"false") | `true` | If true: `n_layer_splits` computed Re-adaptively to hit y+_target at the first cell centre (Dresden convention). False: use the fixed value below. |
| `y+ target` | float | `40.0` | Target y+ for auto-layers. 40 = wall-function regime (robust+cheap). ~1 = wall-resolved (more expensive). |
| `wall-layer splits` | integer | `0` | Fixed `refineWallLayer` passes (used only when auto_layers=false). ~4 → y+<10 at Re=1e5, ny=40. |
| `layer-split thickness fraction` | float | `0.5` | Each pass halves the wall-adjacent cell at this fraction. |

### Run-time

| Field | Type | Default | Notes |
|---|---|---|---|
| `number of processors` | integer | `4` | MPI ranks. Capped at 128. |
| `contact email` | string | (auto-prefilled in portal) | Ignored by CLI. |
| `comments` | text_area | empty | Free-text notes. |

---

## Adding a new field

If you need a new design or fixed parameter, edit the workflow's
`interface.json` -- add a new 6-tuple entry with a stable position
number, type, and tooltip. Then teach the worker
(`py_<workflow>.py`) to read it. The CLI picks the change up
immediately (no portal rebuild needed for headless use).

Coordinate with Martin if the new field should land in the Streamlit
portal as well -- the portal redeploy is a separate step.
