#!/usr/bin/env python3
"""Generate the A3 pressure-swirl nozzle mesh inputs into a case directory.

Drives the validated variant #6 generators (sketch3d.py O-grid + make_snap.py
snappy layer, PRISMATIC=True, box refinement) for ONE case. Geometry/resolution
come from params.py, which reads the A3_<param> env overrides at import time --
so this MUST run as a FRESH subprocess per case (the py_A3_nozzle worker sets the
env, then calls `python gen_a3_mesh.py <case_dir>`).

Writes into <case_dir>:
    system/blockMeshDict              (O-grid, named blocks)
    system/snappyHexMeshDict          (carve the n inlets, box-refine the chamber)
    system/surfaceFeatureExtractDict
    system/meshQualityDict
    constant/triSurface/nozzle.stl    (bore wall + inlets, regions wall/inlet)
    constant/triSurface/feat_inlet.stl

It does NOT write controlDict/fvSchemes/fvSolution -- those are the interFoam VOF
dicts shipped in case_template/ and must survive. sketch3d.write_blockmeshdict
would emit its own minimal ones, so the blockMeshDict is generated in a temp dir
and only that file is copied across.

Requires trimesh + manifold3d (for the unioned nozzle.stl), numpy.
"""
import os
import shutil
import sys
import tempfile

import sketch3d   # noqa: E402  computes the O-grid geometry from params (A3_* env)
import make_snap  # noqa: E402  snappy layer (PREREFINE_LEVEL=0 -> validated box mode)


def generate(case_dir: str) -> None:
    case_dir = os.path.abspath(case_dir)
    sysdir = os.path.join(case_dir, "system")
    tri = os.path.join(case_dir, "constant", "triSurface")
    os.makedirs(sysdir, exist_ok=True)
    os.makedirs(tri, exist_ok=True)

    # blockMeshDict only (sketch3d also drops minimal controlDict/fvSchemes/
    # fvSolution -> write to a throwaway dir, copy across just the blockMeshDict
    # so the interFoam VOF system dicts from case_template are preserved).
    tmp = tempfile.mkdtemp(prefix="a3bm_")
    try:
        sketch3d.write_blockmeshdict(tmp)
        shutil.copy(os.path.join(tmp, "system", "blockMeshDict"),
                    os.path.join(sysdir, "blockMeshDict"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # STLs (nozzle bore + inlets + ambient top-seal disk gemerged) + snappy /
    # feature / quality dicts. Mit AMBIENT ist die top-seal-Disk ins nozzle.stl
    # gemerged (boolean union, inner cone follows bore wall) -> kein separates
    # top.stl, kein topSeal-Patch in der polyMesh, kein createPatch noetig.
    make_snap.build_surfaces(tri)
    make_snap.write_sfe(sysdir)
    make_snap.write_meshquality(sysdir)
    make_snap.write_snappy(sysdir)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python gen_a3_mesh.py <case_dir>")
    generate(sys.argv[1])
    print("gen_a3_mesh OK -> %s" % os.path.abspath(sys.argv[1]))
