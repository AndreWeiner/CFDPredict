#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parametric blockMeshDict generator for A2 Brand-topology variant.

Variant of the A2_leitbleche elbow with parallel-translated bend arcs
(instead of concentric ones). Geometrically closer to the original Brand
2020 figures with full-height baffles separating constant-width channels
through the bend.

Coordinate system identical to A2-mother (gen_blockmesh.py):
    Inlet flows in +x from x=-L_in to x=0.
    Bend center at (0, R); pre-bend channel at y in [-H/2, +H/2] (centerline y=0).
    Post-bend channel at x in [R-H/2, R+H/2], y in [some_y, R+L_out].
    Inner wall (closer to bend center (0,R)) = y=+H/2 pre-bend, x=R-H/2 post-bend.
    Outer wall (further from bend center)    = y=-H/2 pre-bend, x=R+H/2 post-bend.

Pathline numbering (5 main vertices per pathline + 1 arc-midpoint per arc):
    Pathline n has vertices p_{5n}..p_{5n+4}:
        pos 0 = Inlet (x = -L_in)
        pos 1 = Bend arc start (just before the bend)
        pos 2 = Bend arc end
        pos 3 = End of Trail block (y = arc_end_y + L_trail)
        pos 4 = Outlet (y = R + L_out)
    Bend arc midpoint = p_{5n+1}_{5n+2}.

Pathlines numbered outer (s=1) -> inner (s=0); inner_arc_fractions insert
additional pathlines between them, sorted descending.

If L_trail == 0, the Trail block collapses (pos 2 == pos 3) and is skipped.

Constraint: r_bend < R - H/2 (geometrically valid; at r_bend = R - H/2 the
inner-arc-end and outer-arc-start can coincide for R = 3H/2, causing block
degeneration).
"""
from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence


# ----------------------------------------------------------------------
# Geometry data class
# ----------------------------------------------------------------------
@dataclass
class BrandTopology:
    H: float = 1.0
    W: float = 1.0
    R: float = 1.5
    L_in: float = 3.0
    L_trail: float = 0.0
    L_out: float = 5.0
    r_bend: float | None = None     # None -> R - H/2 (tight fit)
    inner_arc_fractions: list[float] = field(default_factory=list)
    N_inlet: int = 16
    nz: int = 16
    scale: float = 1.0  # multiplier applied at writing time via 'scale'

    # ------------------------------------------------------------------
    def _r(self) -> float:
        """Resolved bend arc radius."""
        return (self.R - self.H / 2.0) if self.r_bend is None else float(self.r_bend)

    def validate(self) -> None:
        if self.H <= 0:        raise ValueError("H must be > 0")
        if self.W <= 0:        raise ValueError("W must be > 0")
        if self.R <= self.H / 2.0:
            raise ValueError(f"R must be > H/2 (got R={self.R}, H/2={self.H/2})")
        if self.L_in <= 0:     raise ValueError("L_in must be > 0")
        if self.L_out <= 0:    raise ValueError("L_out must be > 0")
        if self.L_trail < 0:   raise ValueError("L_trail must be >= 0")
        if self.L_trail > self.L_out:
            raise ValueError(f"L_trail ({self.L_trail}) must be <= L_out ({self.L_out})")
        r = self._r()
        if not (0 < r <= self.R - self.H / 2.0 + 1e-12):
            raise ValueError(
                f"r_bend ({r}) must satisfy 0 < r_bend <= R-H/2 ({self.R-self.H/2:.6f})"
            )
        for f in self.inner_arc_fractions:
            if not (0.0 < f < 1.0):
                raise ValueError(
                    f"inner_arc_fractions must be strictly in (0,1); got {f}")
        if self.N_inlet < 2:   raise ValueError("N_inlet must be >= 2")
        if self.nz < 1:        raise ValueError("nz must be >= 1")

    def derived(self) -> dict:
        """Return a dict of derived quantities (arc lengths, pathlines, ...)."""
        self.validate()
        r = self._r()
        s_values = sorted([1.0] + list(self.inner_arc_fractions) + [0.0],
                          reverse=True)  # outer (s=1) -> inner (s=0)
        return {
            "r_bend": r,
            "arc_length": r * math.pi / 2.0,
            "pathline_s_values": s_values,
            "N_pathlines": len(s_values),
            "N_lanes": len(s_values) - 1,
            "cell_size": self.H / self.N_inlet,
            "has_trail": self.L_trail > 0.0,
        }


# ----------------------------------------------------------------------
# Point generation in A2 coordinate system
# ----------------------------------------------------------------------
def build_points(geo: BrandTopology):
    """Return (pts, arcs_edges, vertex_offset_z, N_pathlines).

    pts:        dict[str, (x, y)] -- 2D coordinates (m), z=0 plane.
    arcs_edges: list[(start_name, mid_on_arc_name, end_name)].
    vertex_offset_z: integer; +offset for the z=W layer at export.

    For each pathline n at fraction s_n (outer=1 at y=-H/2, inner=0 at y=+H/2):
      pos 0 (Inlet):     (-L_in,                                (0.5-s)*H)
      pos 1 (Arc start): (R + (s-0.5)*H - r,                    (0.5-s)*H)
      pos 2 (Arc end):   (R + (s-0.5)*H,                        (0.5-s)*H + r)
      pos 3 (Trail end): (R + (s-0.5)*H,                        arc_end_y + L_trail)
      pos 4 (Outlet):    (R + (s-0.5)*H,                        R + L_out)

    The post-bend wall x = R + (s-0.5)*H interpolates between R-H/2 (inner, s=0)
    and R+H/2 (outer, s=1). The pre-bend wall y = (0.5-s)*H interpolates between
    +H/2 (inner) and -H/2 (outer). Centerline of channel: y=0 pre-bend, x=R post-bend.
    """
    d = geo.derived()
    r = d["r_bend"]
    s_values = d["pathline_s_values"]
    H = geo.H
    R = geo.R
    L_in = geo.L_in
    L_trail = geo.L_trail
    L_out = geo.L_out

    pts = {}
    arcs_edges = []

    for n, s in enumerate(s_values):
        x_post_wall = R + (s - 0.5) * H               # post-bend wall x (in [R-H/2, R+H/2])
        y_pre_wall  = (0.5 - s) * H                   # pre-bend wall y  (in [-H/2, +H/2])
        x_arc_start = x_post_wall - r                 # where the bend arc starts
        y_arc_end   = y_pre_wall + r                  # where the bend arc ends
        x_arc_center = x_arc_start                    # center on +y side of start
        y_arc_center = y_arc_end
        y_outlet = R + L_out
        y_trail_end = y_arc_end + L_trail

        base = 5 * n
        pts[f"p{base + 0}"] = (-L_in,         y_pre_wall)
        pts[f"p{base + 1}"] = (x_arc_start,   y_pre_wall)
        pts[f"p{base + 2}"] = (x_post_wall,   y_arc_end)
        pts[f"p{base + 3}"] = (x_post_wall,   y_trail_end)
        pts[f"p{base + 4}"] = (x_post_wall,   y_outlet)

        # Arc sweeps from start (angle 3pi/2 from center, tangent +x) to end
        # (angle 0 from center, tangent +y), CCW. Midpoint at angle 7pi/4 (315 deg).
        mid_angle = 7.0 * math.pi / 4.0
        x_mid = x_arc_center + r * math.cos(mid_angle)
        y_mid = y_arc_center + r * math.sin(mid_angle)
        mid_key = f"p_{base + 1}_{base + 2}"
        pts[mid_key] = (x_mid, y_mid)
        arcs_edges.append((f"p{base + 1}", mid_key, f"p{base + 2}"))

    vertex_offset_z = 5 * len(s_values)
    return pts, arcs_edges, vertex_offset_z, len(s_values)


# ----------------------------------------------------------------------
# Block + patch generation (identical pattern to brand_classic sketch.py)
# ----------------------------------------------------------------------
SEGMENT_NAMES = ["Feed", "Bend", "Trail", "Run"]  # pos 0->1, 1->2, 2->3, 3->4


def build_topology(geo: BrandTopology):
    """Return (pts, arcs_edges, blocks, patches, vertex_offset_z, derived)."""
    pts, arcs_edges, voz, N_pathlines = build_points(geo)
    d = geo.derived()

    cell_size = d["cell_size"]
    has_trail = d["has_trail"]
    arc_length = d["arc_length"]
    s_values = d["pathline_s_values"]

    # cells per lane in width direction
    n_width = []
    for k in range(N_pathlines - 1):
        s_diff = s_values[k] - s_values[k + 1]   # always > 0 (descending order)
        n_width.append(int(math.ceil(geo.N_inlet * s_diff)))

    # cells per longitudinal segment (same across all lanes; from longest pathline)
    def seg_len(n: int, pos: int) -> float:
        if pos == 1:
            return arc_length
        a = pts[f"p{5*n + pos}"]
        b = pts[f"p{5*n + pos + 1}"]
        return math.hypot(b[0] - a[0], b[1] - a[1])

    # Determine which segments are active. If L_trail == 0, the Trail block
    # (segment pos 2->3) collapses -- skip it and let Run start at pos 2.
    active_segments = [0, 1]                       # Feed, Bend always present
    if has_trail:
        active_segments.append(2)                  # Trail
    active_segments.append(3)                      # Run

    # n_long per ACTIVE segment
    n_long = {}
    for pos in active_segments:
        max_len = max(seg_len(n, pos) for n in range(N_pathlines))
        n_long[pos] = int(math.ceil(max_len / cell_size))

    # When trail is collapsed, the Run block goes from pos 2 directly to pos 4.
    # We must adjust segment representation accordingly.

    # n_z
    n_z = max(1, int(math.ceil(geo.W / cell_size)))

    blocks = []
    for lane_idx in range(N_pathlines - 1):
        n_lo = lane_idx
        n_hi = lane_idx + 1

        # Helper: a single hex block name and vertex tuple.
        # CCW order (viewed from +z) so that blockMesh's right-hand rule yields
        # outward face normals. In A2 frame, outer wall (n_lo, s=1) is at LOWER y
        # pre-bend (y=-H/2) and HIGHER x post-bend (x=R+H/2), so swapping start
        # with inner (n_hi) gives the correct CCW orientation:
        #     v0 = inner_pos_start (upper-left in pre-bend view)
        #     v1 = outer_pos_start (lower-left)
        #     v2 = outer_pos_end   (lower-right)
        #     v3 = inner_pos_end   (upper-right)
        def add_block(pos_start: int, pos_end: int, seg_name: str, ny_cells: int):
            v0 = f"p{5*n_hi + pos_start}"
            v1 = f"p{5*n_lo + pos_start}"
            v2 = f"p{5*n_lo + pos_end}"
            v3 = f"p{5*n_hi + pos_end}"
            blocks.append({
                "name":    f"{seg_name}_Lane{lane_idx + 1}",
                "verts":   (v0, v1, v2, v3),
                "cells":   (n_width[lane_idx], ny_cells, n_z),
                "grading": (1.0, 1.0, 1.0),
            })

        # Feed (pos 0 -> 1)
        add_block(0, 1, "Feed", n_long[0])
        # Bend (pos 1 -> 2)
        add_block(1, 2, "Bend", n_long[1])
        if has_trail:
            # Trail (pos 2 -> 3), then Run (pos 3 -> 4)
            add_block(2, 3, "Trail", n_long[2])
            add_block(3, 4, "Run", n_long[3])
        else:
            # Skip Trail; Run goes pos 2 -> 4 directly. We need its own
            # ny_cells based on the actual Run length:
            run_lens = []
            for n in range(N_pathlines):
                a = pts[f"p{5*n + 2}"]
                b = pts[f"p{5*n + 4}"]
                run_lens.append(math.hypot(b[0] - a[0], b[1] - a[1]))
            ny_run = int(math.ceil(max(run_lens) / cell_size))
            add_block(2, 4, "Run", ny_run)

    # Patches: inlet (at x = -L_in) and outlet (at y = R + L_out)
    # In A2 frame: outer wall pre-bend at y = -H/2, inner pre-bend at y = +H/2.
    #              outer wall post-bend at x = R+H/2, inner post-bend at x = R-H/2.
    # Inlet  plane normal = -x. Edge in -y direction gives -x normal (right-hand rule
    #        with +z extrusion). At pos 0: inner (y=+H/2) -> outer (y=-H/2) is -y. ✓
    # Outlet plane normal = +y. Edge in -x direction gives +y normal.
    #        At pos 4: outer (x=R+H/2) -> inner (x=R-H/2) is -x. ✓
    patches = {"patch inlet": [], "patch outlet": []}
    for lane_idx in range(N_pathlines - 1):
        n_lo = lane_idx         # outer (higher s, y=-H/2 side)
        n_hi = lane_idx + 1     # inner (lower  s, y=+H/2 side)
        patches["patch inlet"].append((f"p{5*n_hi + 0}", f"p{5*n_lo + 0}"))
        patches["patch outlet"].append((f"p{5*n_lo + 4}", f"p{5*n_hi + 4}"))

    return {
        "pts": pts,
        "arcs_edges": arcs_edges,
        "blocks": blocks,
        "patches": patches,
        "vertex_offset_z": voz,
        "N_pathlines": N_pathlines,
        "n_width": n_width,
        "n_long": n_long,
        "n_z": n_z,
        "derived": d,
    }


# ----------------------------------------------------------------------
# blockMeshDict writer
# ----------------------------------------------------------------------
_HEADER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
|  \\\\      /  F ield        | OpenFOAM                                        |
|   \\\\    /   O peration    | CFDPredict A2 -- Leitbleche / brand_topology   |
|    \\\\  /    A nd          | generated by gen_blockmesh.py                  |
|     \\\\/     M anipulation |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def write_block_mesh_dict(geo: BrandTopology, path: str | os.PathLike) -> str:
    """Write blockMeshDict to <path>. Returns absolute path."""
    topo = build_topology(geo)
    pts = topo["pts"]
    arcs_edges = topo["arcs_edges"]
    blocks = topo["blocks"]
    patches = topo["patches"]
    voz = topo["vertex_offset_z"]
    W = geo.W

    # Sort main points by integer index
    main_keys = sorted((k for k in pts if not k.startswith("p_")),
                       key=lambda k: int(k[1:]))
    pt_index = {k: i for i, k in enumerate(main_keys)}
    offset = voz

    out: list[str] = []
    out.append(_HEADER)
    out.append("scale   1.0;   // coordinates already in meters")
    out.append("")

    out.append("vertices")
    out.append("(")
    for k in main_keys:
        x, y = pts[k]
        out.append(f"    ( {x:14.6e} {y:14.6e} {0.0:14.6e} )  // "
                   f"{pt_index[k]:3d}  {k}")
    for k in main_keys:
        x, y = pts[k]
        out.append(f"    ( {x:14.6e} {y:14.6e} {W:14.6e} )  // "
                   f"{pt_index[k] + offset:3d}  {k}+offset")
    out.append(");")
    out.append("")

    out.append("blocks")
    out.append("(")
    for blk in blocks:
        v0, v1, v2, v3 = blk["verts"]
        nx, ny, nz_ = blk["cells"]
        gx, gy, gz = blk["grading"]
        ix = [pt_index[v0], pt_index[v1], pt_index[v2], pt_index[v3]]
        all8 = ix + [i + offset for i in ix]
        v_str = " ".join(str(i) for i in all8)
        out.append(
            f"    hex ({v_str}) {blk['name']} "
            f"({nx} {ny} {nz_}) simpleGrading ({gx} {gy} {gz})"
        )
    out.append(");")
    out.append("")

    out.append("edges")
    out.append("(")
    for layer in (0, 1):
        z_val = 0.0 if layer == 0 else W
        idx_off = 0 if layer == 0 else offset
        for (s_name, mid_name, e_name) in arcs_edges:
            i = pt_index[s_name] + idx_off
            j = pt_index[e_name] + idx_off
            mx, my = pts[mid_name]
            out.append(f"    arc {i} {j} ({mx:14.6e} {my:14.6e} {z_val:14.6e})")
    out.append(");")
    out.append("")

    out.append("boundary")
    out.append("(")
    for patch_key, edge_list in patches.items():
        ptype, pname = patch_key.split(" ", 1)
        out.append(f"    {pname}")
        out.append("    {")
        out.append(f"        type {ptype};")
        out.append("        faces")
        out.append("        (")
        for (a, b) in edge_list:
            ia = pt_index[a]
            ib = pt_index[b]
            out.append(f"            ( {ia:3d} {ib:3d} "
                       f"{ib + offset:3d} {ia + offset:3d} )")
        out.append("        );")
        out.append("    }")
    # In quasi-2D (nz=1) we need an explicit empty front/back; for 3D the
    # defaultPatch=walls catches it. For now we let blockMesh emit them via
    # defaultPatch (works for both, becomes wall in 3D and we patch it to
    # empty post-hoc in 2D mode -- TODO when 2D is needed).
    out.append(");")
    out.append("")
    out.append("defaultPatch")
    out.append("{")
    out.append("    name walls;")
    out.append("    type wall;")
    out.append("}")
    out.append("")
    out.append("mergePatchPairs")
    out.append("(")
    out.append(");")
    out.append("")
    out.append("// "
               + "*" * 73 + " //")

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(out), encoding="utf-8", newline="\n")
    return str(p.resolve())


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def _from_params_json(path: str) -> BrandTopology:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    p = data["parameters"]
    def v(name): return p[name]["value"]
    return BrandTopology(
        H=v("H"), W=v("W"), R=v("R"),
        L_in=v("L_in"), L_trail=v("L_trail"), L_out=v("L_out"),
        r_bend=v("r_bend"),
        inner_arc_fractions=list(v("inner_arc_fractions")),
        N_inlet=v("N_inlet"), nz=v("nz"),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--params", default="params.json",
                    help="Pfad zu params.json (Default: ./params.json)")
    ap.add_argument("--out", default=None,
                    help="Ziel-blockMeshDict (Default: system/blockMeshDict relativ zu cwd)")
    args = ap.parse_args()

    geo = _from_params_json(args.params)
    out_path = args.out or "system/blockMeshDict"
    written = write_block_mesh_dict(geo, out_path)
    topo = build_topology(geo)
    d = topo["derived"]
    print(f"r_bend             = {d['r_bend']:.6f} m  "
          f"(R-H/2 = {geo.R - geo.H/2:.6f})")
    print(f"arc_length         = {d['arc_length']:.6f} m")
    print(f"N_pathlines        = {d['N_pathlines']}")
    print(f"N_lanes            = {d['N_lanes']}")
    print(f"cell_size          = {d['cell_size']:.6f} m  "
          f"(N_inlet = {geo.N_inlet})")
    print(f"has_trail          = {d['has_trail']}  (L_trail = {geo.L_trail})")
    print(f"n_width            = {topo['n_width']}")
    print(f"n_long             = {topo['n_long']}")
    print(f"n_z                = {topo['n_z']}  (W = {geo.W})")
    print(f"blocks ({len(topo['blocks'])}):")
    for blk in topo["blocks"]:
        print(f"  {blk['name']:16s} verts={blk['verts']}  cells={blk['cells']}")
    print(f"patches:")
    for pk, edges in topo["patches"].items():
        print(f"  {pk:14s} ({len(edges)} face(s))")
    print(f"---")
    print(f"geschrieben: {written}")


if __name__ == "__main__":
    main()
