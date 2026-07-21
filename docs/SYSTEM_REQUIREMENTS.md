# System Requirements -- streamlit-demo-workflows (headless mode)

Reference for the **non-Streamlit** use of the workflows: the Dresden
Bayes-loop driver, RunGui-next local integration, and CI smoke tests.
The Streamlit portal manages all of this on its own deployed VMs --
read this only if you want to run a workflow outside the portal.

## Operating system

| Platform | Status | Notes |
|---|---|---|
| **Linux x86_64** | **Primary target** | Tested: openSUSE Tumbleweed, Ubuntu 24.04, Rocky 9. All HPC clusters fall here. |
| **Windows 11 + blueCFD-Core** | Demo / test target | Workflow Python runs natively; OpenFOAM is sourced from blueCFD-Core's bundled `etc/bashrc`. Slower; not for production sweeps. |
| **macOS** | Untested | OpenFOAM via Docker is the only realistic path. Not on our roadmap. |

## OS-level dependencies

The workflows expect a working **Python 3.10+**, **OpenFOAM ESI v2406+**
and a handful of common build/visualisation utilities. The lists below
are tested on the distros mentioned. On clusters using `module load`,
substitute the equivalent module names; the package list is the same.

### Ubuntu / Debian

```bash
# System tools + Python + ParaView
sudo apt update
sudo apt install -y python3 python3-pip python3-venv \
                    gnuplot imagemagick paraview python3-paraview \
                    jq bc rsync git unzip xvfb

# Python libs (system-wide; alternatively use a venv -- see below)
sudo pip install --break-system-packages \
    matplotlib numpy scipy reportlab pillow networkx shapely \
    trimesh manifold3d mapbox-earcut pyvista vtk pyyaml
```

### openSUSE Leap 16 / Tumbleweed

```bash
# System tools + Python + ParaView
sudo zypper -n install python313 python313-pip python313-PyYAML \
                       paraview gnuplot ImageMagick jq bc rsync git unzip

# *** Caveat for Leap 16 (verified 2026-06-11): the shipped libexpat1 may
# be older than the one Python 3.13 was built against -- pip will then
# crash with `ImportError: /usr/lib64/python3.13/lib-dynload/
# pyexpat.cpython-313-x86_64-linux-gnu.so: undefined symbol:
# XML_SetAllocTrackerActivationThreshold`. Fix is one zypper update:
sudo zypper -n update libexpat1

# Python libs (system-wide; PEP-668 requires the --break-system-packages flag)
sudo pip install --break-system-packages \
    matplotlib numpy scipy reportlab pillow networkx shapely \
    trimesh manifold3d mapbox-earcut pyvista vtk
```

### Fedora / RHEL / Rocky 9

```bash
# System tools + Python + ParaView
sudo dnf install -y python3 python3-pip paraview python3-paraview \
                    gnuplot ImageMagick jq bc rsync git unzip xorg-x11-server-Xvfb

# Python libs
sudo pip install --break-system-packages \
    matplotlib numpy scipy reportlab pillow networkx shapely \
    trimesh manifold3d mapbox-earcut pyvista vtk pyyaml
```

### Venv alternative (recommended for multi-user / multi-project hosts)

If `--break-system-packages` is not your preference (or your distro blocks
it), drop the libs into an isolated env and point each workflow's
`workflow.yaml` `python_exec` at it:

```bash
python3 -m venv /opt/streamlit-workflows-venv
source /opt/streamlit-workflows-venv/bin/activate
pip install -r workflows/A2_brand_topology/requirements.txt
pip install -r workflows/A3_nozzle/requirements.txt
# repeat per workflow you use; or pip install the union of all requirements.txt
```

Then in each workflow's `workflow.yaml`:

```yaml
python_exec: /opt/streamlit-workflows-venv/bin/python
```

### Why each dep matters

| Package | Used by | Hard requirement? |
|---|---|---|
| `matplotlib` | every workflow (charts) | yes |
| `numpy`, `scipy` | A3 mesh generator, validation utils | A3: yes; others: optional but common |
| `pyyaml` | `tools/run_workflow.py` CLI (parses `workflow.yaml`) | only for stand-alone CLI use |
| `reportlab` | dean90Bend PDF report | only if you use dean90Bend |
| `pillow` | matplotlib backend | yes (pulled in transitively) |
| `trimesh`, `manifold3d`, `mapbox-earcut`, `shapely`, `networkx` | A3 STL boolean union | A3: yes |
| `pyvista`, `vtk` | not used at runtime, but handy for post-hoc case inspection | optional |
| `gnuplot` | RunGui `logfile_observer` (live residual charts) | optional (worker continues without it) |
| `imagemagick` | several workflows compose PNGs | optional |
| `xvfb` / `xvfb-run` | only needed if `pvbatch` runs on a host without an X server / VNC session | optional |
| `jq`, `bc` | defensive: some Allrun-style scripts parse JSON / do shell math | nice-to-have |

## Per-workflow Python dependencies

Each workflow ships its own `requirements.txt` -- a tighter list of only
what *that* worker really imports:

```bash
pip install -r workflows/A2_brand_topology/requirements.txt
pip install -r workflows/A3_nozzle/requirements.txt
```

If you only run A2 workflows you don't need the trimesh stack; the
A2 `requirements.txt` is `matplotlib` only. The OS-level table above
is the **union** for installations that may run any workflow.

Visualisation (`pv_*.py`) uses ParaView's bundled Python, NOT this pip
stack.

## OpenFOAM

The worker sources an OpenFOAM `etc/bashrc` lazily, at solver-spawn time.
Resolution order (see `workflows/A2_brand_topology/py_OF_utils.py`):

1. `STREAMLIT_OPENFOAM_BASHRC` env var (absolute path; admin override).
2. ESI install chain: `/opt/OpenFOAM/OpenFOAM-{v2512,v2506,v2406,v2312}/etc/bashrc`.
3. Foundation packaging: `/usr/lib/openfoam/openfoam{2512,12,11,10,9}/etc/bashrc`.
4. blueCFD-Core (Windows only): `C:/Program Files/blueCFD-Core-*/OpenFOAM-*/etc/bashrc`.

**Tested versions:** ESI v2512 (primary) and v2506. Foundation builds and
older ESI versions are best-effort -- dictionary syntax of the case
template is v2406+ compatible. If your cluster module-system needs a
custom path, set `STREAMLIT_OPENFOAM_BASHRC` in your job script.

### MPI

The worker calls `mpirun -np <N> -parallel`. Whatever your OpenFOAM
install provides (OpenMPI 4.x/5.x / MPICH / IntelMPI bundled with ESI)
works -- just make sure the `etc/bashrc` source step puts `mpirun` on
the `PATH`.

## CPU sizing

Auto-detected physical cores, capped at 128 (`MAX_AUTO_NPROCS` in
`py_OF_utils.py`). Override per-job via the workflow's "number of
processors" input, or globally via `STREAMLIT_NPROCS=<N>` env var.

For the Bayes loop, a typical Brand-Pos N-sweep case is ~1-5 minutes on
8 cores at the demo mesh resolution (cells across inlet = 16, nz=4).
Cluster jobs at production resolution scale by mesh^(4/3).

## Optional: RunGui logfile_observer (live residual charts)

The Worker spawns the rungui `logfile_observer` to render live residual
and functionObject PNGs into `progress/` during the solve. **Optional**:
if the observer isn't found, the worker logs a clear miss and continues
without live charts -- the final `summary.json` / per-case result PNGs
are produced regardless.

The resolver tries two layouts:

1. **rungui-classic** (single-file): `which("logfile_observer.py")` over
   `$PATH`, then `$RUNGUI_PATH` → `/opt/rungui-portal` → `/opt/rungui` →
   `~/rungui` → `/usr/local/share/rungui`. Found → spawned as
   `python3 <file> <log> <interval>`.
2. **rungui-next** (Python package, `__main__.py`): probes
   `<prefix>/logfile_observer/__main__.py` for `$RUNGUI_PATH` →
   `/opt/dhcae/rungui-next` → `/opt/rungui-next` → `~/dhcae/rungui-next`
   → `~/rungui-next`. Found → spawned as `python3 -m logfile_observer
   <log> <interval>` with `PYTHONPATH=<prefix>`.

`$RUNGUI_PATH` is checked in both modes -- set it for custom installs.

Requires `gnuplot` ≥ 5 on the host (rungui v9.3+ wants v6).

## Optional: ParaView (visualisation)

Required only if you want the `pv_*.py` PNG renders. Install ParaView 5.11
or newer with `pvbatch` on `PATH`. The PV 6.0+ `DecomposePolyhedra`
attribute change is already handled with a try/except. The portal's VMs
ship PV 5.11 or PV 6.0.1; both work.

Cluster compute-nodes that only run the solver don't need ParaView -- the
worker swallows pvbatch-missing errors and continues. The chart PNG
(`progress/brandTopology_dP_uniformity.png`) is produced by matplotlib
and does NOT need ParaView.

## Disk

Per-case the worker produces:

- `system/`, `constant/`, `0/`, `<latestTime>/` (~ 1-20 MB depending on mesh)
- `postProcessing/` with the surfaceFieldValue functionObjects (~ 50 kB)
- `logs/` with one log file per OF utility + `solver.log` (~ 500 kB)

A 13-case Brand-7.3.3 sweep at demo resolution stays under 200 MB total.

## Smoke test

After install, validate without OpenFOAM:

```bash
cd workflows/A2_brand_topology
python test_build_offline.py     # builds 8 cases dry-run, no solver
```

For an end-to-end smoke (requires OpenFOAM), any of the deployed
settings JSON files works as input -- e.g. one of the bundled
`examples/brand_kap_7_3_3_anzahl_sweep.json`:

```bash
python tools/run_workflow.py A2_brand_topology <path-to-settings.json> \
                             --name smoke-test --timeout 1800
```

See `tools/run_workflow.py --help` for all options.

A curated `examples/` set of per-chapter settings JSON files will be
shipped with the Dresden partner bundle (Phase 1).
