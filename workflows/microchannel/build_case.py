#!/usr/bin/env python3
"""Build a complete thermalShell OpenFOAM case from case_template/.

Port of the parametrization repo's main.py (build half only), driven by
params.CaseParams so py_microchannel.py can call it directly. No OpenFOAM
utility runs here: the template's Allrun / Allrun-parallel do extrudeMesh
(3 graded stages) + topoSet/createPatch, makeFaMesh and the solver.

Workflow (the caller has already passed both gates of
mesher/true_geometry.py):
    1. Mesh the 2D vertical-cut cross-section (mesher/CutMesh.py) and write
       it, scaled to metres, as constant/triSurface/mesh.obj.
    2. Write extrudeMeshDict.stage{1,2,3} (fine-coarse-fine Z grading).
    3. Substitute the physics/run parameters into the copied template's
       @PLACEHOLDER@ tokens.

Shape parameters are in mm, as in the parametrization repo.
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

import numpy as np

from params import CaseParams

HERE = Path(__file__).resolve().parent
TEMPLATE_DIR = HERE / "case_template"
# the mesher modules import each other flat (from true_geometry import ...)
sys.path.insert(0, str(HERE / "mesher"))

from geometry_figures import plot_full_pitch, plot_h_profile  # noqa: E402
from CutMesh import build_cut_mesh  # noqa: E402
from cut_geometry import cut_height  # noqa: E402
from Connect import Connect  # noqa: E402
from LineDistribution import LineDistribution  # noqa: E402
from mesh import Mesh  # noqa: E402
from core import take_gmsh_log  # noqa: E402

# Shape parameters are in mm; extrudeMeshDict has no blockMesh-style
# `scale`, so the surface handed to extrudeMesh is pre-scaled to metres.
BLOCKMESH_SCALE = 0.001

# Z layer count tracks l at 1.5 cells/mm (the reference blockMeshDict's 90
# layers over 60 mm), scaled by the refinement knob.
REFERENCE_LENGTH_Z = 60.0
REFERENCE_BLOCK_Z = 90.0

# Fine-coarse-fine Z grading, (length_fraction, cell_fraction, expansion):
# the reference's 0.30/0.40/0.30 splits as 3 chained extrudeMesh stages.
GRAD_Z_DIVISIONS = [(0.2, 0.30, 10.0), (0.6, 0.40, 1.0), (0.2, 0.30, 0.1)]

# Report slices (controlDict "slices", pv_slices.py): z positions as
# fractions of l. They land exactly on Z layer faces (stage boundaries at
# 0.2/0.8 l, whole cells in the uniform middle stage), where cuttingPlane is
# ill-defined, so the planes are nudged by SLICE_NUDGE * l (11 um at l = 110).
SLICE_Z_FRACS = (0.2, 0.4, 0.6, 0.8)
SLICE_NUDGE = 1e-4

PLACEHOLDER = re.compile(r"@([A-Z][A-Z0-9_]*)@")


def format_float(value: float) -> str:
    return f"{value:.9f}".rstrip("0").rstrip(".")


# ---------------------------------------------------------------------------
# template substitution
# ---------------------------------------------------------------------------
def template_values(p: CaseParams) -> dict[str, dict[str, str]]:
    """Per template file, the value of every @PLACEHOLDER@ it holds."""
    # s is the FULL foil thickness (geometry, channel count and gate use it
    # as is); each foil separates two channels, so this channel's shell
    # carries half of it.
    shell = format_float(0.5 * p.s * BLOCKMESH_SCALE)
    slice_z = {f"Z{i}": format_float((frac + SLICE_NUDGE) * p.l
                                     * BLOCKMESH_SCALE)
               for i, frac in enumerate(SLICE_Z_FRACS, start=1)}
    return {
        "0.orig/finite-area/h_outerShell": {"SHELL_THICKNESS": shell},
        "0.orig/finite-area/Ts_outerShell": {"T_INIT": format_float(p.T_init)},
        # hot gas enters at the inlet, cold at the outlet
        "0.orig/T": {"SHELL_THICKNESS": shell,
                     "T_HOT": format_float(p.T_hot),
                     "T_COLD": format_float(p.T_cold)},
        # the shared include of U and p_rgh
        "0.orig/flowOscillationParams": {"PERIOD": format(p.period, ".9g"),
                                         "MDOT_IN": format(p.mdot, ".9g")},
        # the report-image FOs ("slices", "shell") write every T/4
        "system/controlDict": {"END_TIME": format_float(p.end_time),
                               "MAX_CO": format_float(p.max_co),
                               "SLICE_INTERVAL": format(p.period / 4.0, ".9g"),
                               **slice_z},
        "system/fvSolution": {"N_CORRECTORS": str(p.n_correctors)},
        "system/decomposeParDict": {"N_PROC": str(p.n_proc)},
    }


def render(path: Path, values: dict[str, str]) -> None:
    """Replace every @NAME@ in `path` by values[NAME], in place. Fails on a
    placeholder without a value and on a value the file has no
    placeholder for, so the template and template_values() cannot drift."""
    text = path.read_text()
    found = set(PLACEHOLDER.findall(text))
    if found != set(values):
        raise ValueError(
            f"{path.name}: placeholders without a value "
            f"{sorted(found - set(values))}, values without a placeholder "
            f"{sorted(set(values) - found)}")
    path.write_text(PLACEHOLDER.sub(lambda m: values[m.group(1)], text))


# ---------------------------------------------------------------------------
# Z extrusion
# ---------------------------------------------------------------------------
def build_graded_extrude_dicts(case_dir: Path, lengthZ_mm: float,
                               n_z: int) -> None:
    """One extrudeMeshDict.stageN per GRAD_Z_DIVISIONS entry: the first
    extrudes the 2D surface, each later one the previous stage's
    "otherSide" patch, so the final boundary is still just
    sides/originalPatch/otherSide. The stage layer counts sum to n_z."""
    counts = LineDistribution.division_counts(
        [n_frac for _, n_frac, _ in GRAD_Z_DIVISIONS], n_z)
    for i, ((len_frac, _, ratio), n) in enumerate(zip(GRAD_Z_DIVISIONS,
                                                      counts)):
        n = int(n)
        thickness_m = (len_frac * lengthZ_mm) * BLOCKMESH_SCALE
        # blockMesh's expansionRatio is last/first cell over the division;
        # extrudeMesh's is the per-layer multiplier.
        per_layer_ratio = ratio ** (1.0 / (n - 1)) if n > 1 else 1.0

        if i == 0:
            source_block = (
                'constructFrom   surface;\n'
                'surface         "constant/triSurface/mesh.obj";\n'
                'flipNormals     true;\n'
            )
        else:
            source_block = (
                'constructFrom   mesh;\nsourceCase      ".";\n'
                'sourcePatches   (otherSide);\nexposedPatchName otherSide;\n'
                'flipNormals     false;\n'
            )
        text = (
            "FoamFile { version 2.0; format ascii; class dictionary; "
            "object extrudeMeshDict; }\n"
            f"{source_block}"
            "extrudeModel    linearDirection;\n"
            f"nLayers         {n};\n"
            f"expansionRatio  {per_layer_ratio!r};\n"
            "linearDirectionCoeffs\n{\n"
            "    direction   (0 0 1);\n"
            f"    thickness   {format_float(thickness_m)};\n"
            "}\n"
            "mergeFaces      false;\n"
            "mergeTol        0;\n"
        )
        (case_dir / "system" / f"extrudeMeshDict.stage{i + 1}").write_text(text)


def write_geometry_figures(target: Path, shape: dict, h_cut: float) -> None:
    """h_profile.png (true fluid vs the CFD vertical cut at h_cut) and
    full_pitch.png (the true geometry) into the case dir. Documentation
    only: a plotting failure is reported, never fatal to the build."""
    a, b, s, R = (float(shape[k]) for k in ("a", "b", "s", "R"))
    for name, plot in (
            ("h_profile.png", lambda p: plot_h_profile(a, b, s, R, p, h=h_cut)),
            ("full_pitch.png", lambda p: plot_full_pitch(a, b, s, R, p))):
        try:
            plot(str(target / name))
        except Exception as e:  # noqa: BLE001
            print("%s failed (non-fatal): %s" % (name, e))


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------
def build(target: Path, p: CaseParams) -> dict:
    """Write a runnable case for `p` into `target` (replaced if it exists).

    `p` must have passed both gates (check_parameter_values, then
    check_true_geometry). Returns the cut height and the 2D/Z cell counts."""
    shape = p.shape
    h_cut = cut_height(p.a, p.b)
    take_gmsh_log()                                  # drop any stale messages
    structured, unstructured = build_cut_mesh(shape, h_cut, p.refinement)
    vertices, connectivity = Connect().connectAllBlocks(structured, unstructured)

    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(TEMPLATE_DIR, target)
    (target / "log.gmsh").write_text("".join(f"{m}\n" for m in take_gmsh_log()))
    write_geometry_figures(target, shape, h_cut)

    # --- geometry: the 2D cross-section (in m) is extrudeMesh's surface ---
    tri_dir = target / "constant" / "triSurface"
    tri_dir.mkdir(parents=True, exist_ok=True)
    Mesh(np.asarray(vertices, dtype=float) * BLOCKMESH_SCALE,
         connectivity).write_obj(str(tri_dir / "mesh.obj"))

    # --- Z extrusion, fine-coarse-fine ---
    n_z = max(len(GRAD_Z_DIVISIONS),
              round(p.l * REFERENCE_BLOCK_Z / REFERENCE_LENGTH_Z
                    * p.refinement))
    build_graded_extrude_dicts(target, p.l, n_z)

    # --- physics / run parameters ---
    for rel, values in template_values(p).items():
        render(target / rel, values)

    for script in ("Allrun", "Allrun-parallel", "Allclean"):
        path = target / script
        if path.is_file():
            path.chmod(path.stat().st_mode | 0o111)

    return {"n_z": n_z, "cells_2d": len(connectivity), "h_cut_mm": h_cut}
