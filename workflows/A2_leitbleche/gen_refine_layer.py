#!/usr/bin/env python3
"""Patch list helper for ``refineWallLayer``.

The OpenFOAM utility ``refineWallLayer`` is CLI-driven (no dict file):

    refineWallLayer '(walls vane1_master ...)' thickness -overwrite

Each call splits the wall-adjacent cell of every listed patch in two, with
the new face at ``thickness * cell_extent`` from the wall. Multiple calls
accumulate refinement at the wall (each split halves the already-thin layer).

Running this post-blockMesh keeps the parametric mesh topology deterministic
(Bayes/ROM-safe) and avoids snappy's layer-collapse variability — refineWallLayer
performs a single deterministic edge split per patch face, no iteration.

``build_case.py`` substitutes the patch list and thickness into Allrun.
"""
from __future__ import annotations

import argparse

try:      # package form (see __init__.py); flat script dir falls back below
    from .gen_blockmesh import ElbowGeometry
except ImportError:
    from gen_blockmesh import ElbowGeometry


def patches_for(g: ElbowGeometry) -> list[str]:
    """Patches whose wall-adjacent cells should get refined.

    - ``walls`` is the outer/inner channel (+ frontAndBack merged in 3D).
    - ``vaneN_master`` / ``vaneN_slave`` are the createBaffles patch pair for vane N.
    """
    patches = ["walls"]
    for i in range(1, len(g.vane_radii) + 1):
        patches.append(f"vane{i}_master")
        patches.append(f"vane{i}_slave")
    return patches


def render_patch_arg(g: ElbowGeometry) -> str:
    """Render patches as a single OpenFOAM tuple string: ``(walls vane1_master ...)``.

    Used by build_case to substitute @LAYER_PATCHES@ in Allrun.
    """
    return "(" + " ".join(patches_for(g)) + ")"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-vanes", type=int, default=3)
    args = p.parse_args()

    g = ElbowGeometry(
        H=1.0, W=1.0, R=1.5, L_in=3.0, L_out=5.0,
        nx_in=1, nx_bend=1, nx_out=1, ny=1, nz=1,
        vane_radii=[1.0] * args.n_vanes,
    )
    print(render_patch_arg(g))


if __name__ == "__main__":
    main()
