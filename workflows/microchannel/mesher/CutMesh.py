r"""Cross-section mesher for the VERTICAL-CUT CFD approximation.

The approximation (the parametrization repo's study/05_domain_description;
README sec. 3.1): the horizontal ``y = h``
meets the true hypotenuse at ``Q``, and the vertical dropped from ``Q`` to
the floor is the domain's right-hand wall. Everything right of ``x = x_h``
is handed to the solid.

Points
------
::

    p3 (0, a-s) ..                     a-s = the fluid apex (`a` is the
      |           ``..                 OUTER envelope; the foil eats s)
    p6 .            ``.
      |   .             `.             the wall p2 -> p3 is the true
      |      .            \            profile: crest fillet, straight
      |        .           |           flank, and -- when the cut is low
    p4 +--------+ p5        | p2       enough -- part of the trough fillet
      |        |           |
    p0 +--------+-----------+ p1       p1 = (x_h, 0), p2 = (x_h, h)
               p7

``h`` (:func:`cut_geometry.cut_height`) is both the cut height and the
depth of the band along the wall. ``p5`` sits where the wall's own offset
``h`` crosses ``y = h``, so it lies on the band's inner edge exactly; it
is generally further than ``h`` from ``p2`` (``h/sin(theta)`` on a flank
cut). See :func:`cut_geometry.solve_cut_geometry`.

Blocks
------
=========  =================  ==========================================
block1     ``p0 p7 p5 p4``    the floor block, from the axis to p7/p5
block2     ``p7 p1 p2 p5``    wraps the p1 corner: floor band + cut band
block3     ``p2 p3 p6 p5``    the band along the wall
block4     ``p4 p5 p6``       the core: a structured block + a gmsh cap
=========  =================  ==========================================

blocks 1-3 are transfinite (TFI) via :class:`BlockMesh`; block4 is
:mod:`core`. p1 is a genuine 90-degree convex corner of the fluid, and
block2 turns it with both of its edges explicitly graded toward the
p1/p2 walls rather than by a layer stack converging into it.

Distribution (:mod:`distribution`)
----------------------------------
One across-band array is reused verbatim on every edge that spans the
band's depth -- ``p0->p4``, ``p7->p5``, ``p1->p2``, ``p3->p6`` -- so the
blocks' interior gridlines are continuous across their shared seams (TFI
blends same-index points; a shape mismatch shows up directly as shear).
The wall ``p2->p3`` is distributed piece by piece (fillets by swept angle,
the flank as the graded slack), and the inner edge ``p5->p6``
(:func:`cut_geometry.inner_band_edge`) is the wall's true offset, sampled
at the wall's own fractions so the cross-lines join like to like.

Edge orientation is load-bearing: ``left`` and ``right`` must run in the
SAME direction (lower -> upper), as must ``lower`` and ``upper``.
Reversing one of a pair twists the block and folds its cells -- the
symptom is a single block returning BOTH windings among its own cells.

Module layout: :mod:`cut_geometry` (points and band edges),
:mod:`distribution` (cell counts and fractions), :mod:`core` (block4),
this module (assembly).
"""
import numpy as np

from BlockMesh import BlockMesh
from LineDistribution import LineDistribution
from core import (CORE_CELL_SIZE_RATIO, as_tuples, build_structured_core,
                  gmsh_fill, plan_core)
from cut_geometry import inner_band_edge, solve_cut_geometry, wall_points
from distribution import (MIN_CELL_COUNT, Sizing, gp_fractions,
                          graded_fractions, merge_trough_sliver,
                          wall_fractions)


def build_cut_mesh(shape, h, refinement):
    """Build the 2D cross-section of the vertical-cut approximation.

    :param shape: ``a``, ``b``, ``s``, ``R`` [mm].
    :param h: the cut height (``cut_geometry.cut_height``) [mm].
    :param refinement: the global refinement knob (1 = reference mesh).
    :return: ``(structured_blocks, unstructured_patches)``, in the form
        ``Connect.connectAllBlocks`` expects.
    """
    geo = solve_cut_geometry(shape, h)
    x_h, p5_offset = geo["x_h"], geo["p5_offset"]
    p0, p1, p2, p3 = geo["p0"], geo["p1"], geo["p2"], geo["p3"]
    p4, p5, p6, p7 = geo["p4"], geo["p5"], geo["p6"], geo["p7"]
    sizing = Sizing.for_band(refinement, h)
    band_fr = sizing.band_fr

    # --- the floor: p0 -> p7 (and p4 -> p5, so block1's lines are vertical)
    # block2's corner edges span p5_offset, graded fine toward the p1/p2
    # wall. With a structured core, block1's floor takes the core's own
    # count n_x, as a GP whose first cell at p7 is block2's last -- the
    # floor carries on from the corner without a jump and coarsens toward
    # the axis; p4 -> p5 is then also the core block's first row.
    corner_fr = 1.0 - graded_fractions(sizing.n_band, p5_offset)[::-1]
    plan = plan_core(geo, sizing.refinement)
    floor_len = x_h - p5_offset
    if plan is not None:
        h_p7 = float(corner_fr[1] - corner_fr[0]) * p5_offset
        floor_fr = 1.0 - gp_fractions(h_p7, floor_len, plan["n_x"])[::-1]
    else:
        n_floor = max(MIN_CELL_COUNT,
                      int(round(sizing.nx_per_mm * floor_len)))
        floor_fr = 1.0 - graded_fractions(n_floor, floor_len)[::-1]

    # --- block1: p0 p7 p5 p4 ---------------------------------------------
    e_p0p7 = LineDistribution.apply_fractions(p0, p7, floor_fr)
    e_p4p5 = LineDistribution.apply_fractions(p4, p5, floor_fr)
    e_p0p4 = LineDistribution.apply_fractions(p0, p4, band_fr)
    e_p7p5 = LineDistribution.apply_fractions(p7, p5, band_fr)
    block1 = BlockMesh()
    block1.transfinite(boundary=[e_p0p7, e_p4p5, e_p0p4, e_p7p5])

    # --- block2: p7 p1 p2 p5 -- wraps the p1 corner -----------------------
    e_p7p1 = LineDistribution.apply_fractions(p7, p1, corner_fr)
    e_p5p2 = LineDistribution.apply_fractions(p5, p2, corner_fr)
    e_p1p2 = LineDistribution.apply_fractions(p1, p2, band_fr)
    block2 = BlockMesh()
    block2.transfinite(boundary=[e_p7p1, e_p5p2, e_p7p5, e_p1p2])

    # --- block3: p2 p3 p6 p5 -- the band along the wall -------------------
    wall = geo["wall"]
    spans = merge_trough_sliver(
        [dict(length=p.length, kappa=p.kappa) for p in wall],
        sizing.band_last_cell)
    wall_fr = wall_fractions(
        spans, sum(p.length for p in wall), sizing.ref_len, width=geo["b"],
        block_rows=(None if plan is None
                    else (plan["n_y_min"], plan["n_y_max"])))
    outer = wall_points(wall, wall_fr)
    outer[0], outer[-1] = np.asarray(p2), np.asarray(p3)
    inner = inner_band_edge(geo, wall_fr, h)

    # The cross edges run WALL -> INNER. The p2 end IS block2's e_p5p2
    # reversed -- the same point sequence, so the seam is shared verbatim.
    e_p2p5 = e_p5p2[::-1]
    e_p3p6 = LineDistribution.apply_fractions(p3, p6, band_fr)
    block3 = BlockMesh()
    block3.transfinite(boundary=[as_tuples(outer), as_tuples(inner),
                                 e_p2p5, e_p3p6])

    # --- block4: p4 p5 p6 -- the core -------------------------------------
    # The crest arc's start is an exact wall station (every piece junction
    # is one); after a vanishing-flank merge it is the station where the
    # two arcs meet, half a flank below T_c -- still the nearest.
    core = None
    if plan is not None:
        L = np.array([sp["length"] for sp in spans], dtype=float)
        j_tc = int(np.argmin(np.abs(wall_fr - L[:-1].sum() / L.sum())))
        core = build_structured_core(list(e_p4p5), inner, floor_fr,
                                     plan["c"], j_tc)
    if core is None:
        # No room for the structured block: gmsh fills the whole core,
        # p4 -> p5 (block1), p5 -> p6 (block3), p6 -> p4 (the axis).
        cell = CORE_CELL_SIZE_RATIO * sizing.ref_len
        n_axis = max(MIN_CELL_COUNT, int(round((p6[1] - p4[1]) / cell)))
        axis = [(0.0, float(y))
                for y in np.linspace(p6[1], p4[1], n_axis + 1)]
        loop = list(e_p4p5)[:-1] + as_tuples(inner)[:-1] + axis[:-1]
        core = gmsh_fill(loop, [cell] * len(loop))

    return [block1, block2, block3], [core]

