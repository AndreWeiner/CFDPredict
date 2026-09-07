# -*- coding: utf-8 -*-
"""
A3 SwirlNozzleInjector -- O-Grid + Quadrat-Aussengitter (2D, xy-Ebene).

Iterative blockMesh-Skizze. Stil/Helper analog examples/brand_classic/sketch.py.

Topologie (von innen nach aussen)
---------------------------------
    Schale 0  Quadrat            p0  p1  p2  p3    Ecken (+-a, +-a)
    Schale 1  Kreis 1            p4  p5  p6  p7    Radius R1, Diagonalrichtungen
    Schale 2  Kreis 2            p8  p9  p10 p11   Radius R2, Diagonalrichtungen
    Schale 3  Quadrat (zentral)  p12 p13 p14 p15   Ecken (+-c, +-c)
    Aussen    4x4-Gitterpunkte   p16 .. p27        Perimeter des 3x3-Rasters

Bloecke
-------
    Core               1   Quadrat-Kern (Schale 0)
    Ring1_{S,E,N,W}    4   Quadrat -> Kreis 1   (arc auf Kreis 1)
    Ring2_{S,E,N,W}    4   Kreis 1 -> Kreis 2   (arc auf Kreis 2)
    Trans_{S,E,N,W}    4   Kreis 2 -> Quadrat   (arc innen = Kreis 2, aussen gerade)
    Out_{SW,S,SE,E,NE,N,NW,W}  8   3x3-Raster, identische Quadrate (2c x 2c)
    --------------------------------------------------------------
    Summe             21 Bloecke

Die 8 Aussen-Quadrate sind CCW ab links-unten nummeriert
(SW, S, SE, E, NE, N, NW, W), passend zur Eckpunkt-Konvention.

Einheiten: mm  ->  scale 0.001 beim blockMeshDict-Export.
"""

import math
import os
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Helper: Kreis durch 3 Punkte + Bogen-Sampling  (aus brand_classic)
# ----------------------------------------------------------------------
def circle_from_3pts(p1, p2, p3):
    ax, ay = p1
    bx, by = p2
    cx, cy = p3
    d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-12:
        raise ValueError("Drei Punkte kollinear -- kein Kreis.")
    ux = ((ax**2 + ay**2) * (by - cy)
          + (bx**2 + by**2) * (cy - ay)
          + (cx**2 + cy**2) * (ay - by)) / d
    uy = ((ax**2 + ay**2) * (cx - bx)
          + (bx**2 + by**2) * (ax - cx)
          + (cx**2 + cy**2) * (bx - ax)) / d
    r = math.hypot(ax - ux, ay - uy)
    return (ux, uy), r


def arc_samples(p_start, p_mid, p_end, n=60):
    (cx, cy), r = circle_from_3pts(p_start, p_mid, p_end)
    a_s = math.atan2(p_start[1] - cy, p_start[0] - cx)
    a_e = math.atan2(p_end[1]   - cy, p_end[0]   - cx)
    a_m = math.atan2(p_mid[1]   - cy, p_mid[0]   - cx)

    def ccw_dist(a, b):
        return (b - a) % (2.0 * math.pi)
    if ccw_dist(a_s, a_m) <= ccw_dist(a_s, a_e):
        end_angle = a_s + ccw_dist(a_s, a_e)
    else:
        end_angle = a_s - ccw_dist(a_e, a_s)
    angles = np.linspace(a_s, end_angle, n)
    return [(cx + r * math.cos(t), cy + r * math.sin(t)) for t in angles]


# ----------------------------------------------------------------------
# 1) Parameter -- aus zentralem params.py (Single Source of Truth)
# ----------------------------------------------------------------------
# Alle Querschnitts-Masze sind dort als Verhaeltnisse zu Zhang-Do abgeleitet
# (a, t1=Uebergang, t2=BL/Film-Band, t3=Buffer, N_core, n_r1, n_r2 = BL-Sizing).
from params import a, t1, t2, t3, N_core, n_r1, n_r2
import params as P
try:
    from params import TOPOLOGY, R3 as _R3_abs
except ImportError:                  # backward compat (pre-v2 params)
    TOPOLOGY = "v1"
    _R3_abs = None

# z-Extrusion (3D) -- Platzhalter, echter z-Aufbau in sketch_side.py.
dz = 2.0 * a


# ----------------------------------------------------------------------
# 2) Abgeleitete Geometrie
# ----------------------------------------------------------------------
R0 = a * math.sqrt(2.0)       # Radius der Quadrat-Ecke (Diagonale)
R1 = R0 + t1                  # Kreis 1
R2 = R1 + t2                  # Kreis 2
c  = R2 + t3                  # halbe Seitenlaenge des zentralen Aussen-Quadrats
                              # (Kantenmitte bei x=c liegt um t3 ausserhalb R2)

# Diagonal-Einheitsrichtungen der 4 Eckpunkte (CCW ab links-unten):
CORNER_DIRS = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
_INV_SQRT2 = 1.0 / math.sqrt(2.0)

# Kardinal-Winkel (deg) fuer den Bogen-Mittelpunkt zwischen Ecke i und i+1:
ARC_MID_ANGLE_DEG = [270.0, 0.0, 90.0, 180.0]


# ----------------------------------------------------------------------
# 3) Punkte
# ----------------------------------------------------------------------
pts = {}
ring_level = {}     # p-Name -> Schalen-Index 0..4  (fuer Zellzahl-Logik)


def _add_corner_shell(base, level, radius_or_halfside, is_square):
    for i, (dx, dy) in enumerate(CORNER_DIRS):
        if is_square:
            pts[f"p{base + i}"] = (dx * radius_or_halfside, dy * radius_or_halfside)
        else:
            pts[f"p{base + i}"] = (dx * _INV_SQRT2 * radius_or_halfside,
                                   dy * _INV_SQRT2 * radius_or_halfside)
        ring_level[f"p{base + i}"] = level


_add_corner_shell(0,  0, a,  is_square=True)    # p0..p3   Quadrat
_add_corner_shell(4,  1, R1, is_square=False)   # p4..p7   Kreis 1
_add_corner_shell(8,  2, R2, is_square=False)   # p8..p11  Kreis 2
_add_corner_shell(12, 3, c,  is_square=True)    # p12..p15 zentrales Quadrat

# Aussen-Perimeter (4x4-Gitter, ohne die inneren 4 = Schale 3) CCW ab (-3c,-3c)
_outer_coords = [
    (-3*c, -3*c), (-c, -3*c), (c, -3*c), (3*c, -3*c),   # p16..p19  unten
    (3*c, -c), (3*c, c), (3*c, 3*c),                     # p20..p22  rechts
    (c, 3*c), (-c, 3*c), (-3*c, 3*c),                    # p23..p25  oben
    (-3*c, c), (-3*c, -c),                               # p26..p27  links
]
for i, xy in enumerate(_outer_coords):
    pts[f"p{16 + i}"] = xy
    ring_level[f"p{16 + i}"] = 4

# TOPOLOGY=v2: EINE zusaetzliche Schale ueber dem 3x3-Aussen-Raster:
#   Schale 5  Ring3 (p28..p31)  -- Frustum (am exit_top R=R3=10mm,
#                                  am amb_far R=R3+L_amb*tan(half_angle)=62mm).
# 2D-Position in pts[] ist die "Reference" am R3_max=R3+L_amb*tan(half_angle).
# sketch3d.py skaliert ring_level=5 Frustum-abhaengig pro Station.
#
# Urspruengliche 2-Ring-Idee (Ring3 konstant + Ring4 Frustum) verworfen, weil
# Ring4-Block-Bottom an exit_top 0-Volume bekommt (alle 4 Vertices auf R=R3).
# Ein gemeinsamer Frustum-Ring vermeidet das.
if TOPOLOGY == "v2" and _R3_abs is not None:
    _R3_max = _R3_abs + (P.z_amb_far - P.z_exit_top) * math.tan(
        math.radians(getattr(P, "cone_half_angle_deg", 60.0)))
    _add_corner_shell(28, 5, _R3_max, is_square=False)   # p28..p31  Ring3 max-corners
    MAIN_KEYS_COUNT = 32
else:
    MAIN_KEYS_COUNT = 28
    _R3_max = None
vertex_offset_z = MAIN_KEYS_COUNT


# ----------------------------------------------------------------------
# 4) Boegen (arc-Kanten) -- nur Kreis 1 und Kreis 2
# ----------------------------------------------------------------------
arcs_edges = []     # (start_name, mid_key, end_name)
arc_lookup = {}     # (a,b) und (b,a) -> mid_key


def _add_ring_arcs(base, R):
    for i in range(4):
        s = f"p{base + i}"
        e = f"p{base + (i + 1) % 4}"
        ang = math.radians(ARC_MID_ANGLE_DEG[i])
        mid_key = f"p_{base + i}_{base + (i + 1) % 4}"
        pts[mid_key] = (R * math.cos(ang), R * math.sin(ang))
        arcs_edges.append((s, mid_key, e))
        arc_lookup[(s, e)] = mid_key
        arc_lookup[(e, s)] = mid_key


_add_ring_arcs(4, R1)   # Kreis 1
_add_ring_arcs(8, R2)   # Kreis 2
if TOPOLOGY == "v2" and _R3_max is not None:
    _add_ring_arcs(28, _R3_max)   # Ring3 (Ref-Position R3_max; Frustum in sketch3d)


# ----------------------------------------------------------------------
# 5) Bloecke
# ----------------------------------------------------------------------
def _signed_area(verts):
    xy = [pts[v] for v in verts]
    s = 0.0
    for i in range(len(xy)):
        x1, y1 = xy[i]
        x2, y2 = xy[(i + 1) % len(xy)]
        s += x1 * y2 - x2 * y1
    return 0.5 * s


def _ccw(verts):
    return tuple(verts) if _signed_area(verts) > 0 else tuple(reversed(verts))


cell_size = (2.0 * a) / N_core
# n_r1 (Uebergang) und n_r2 (BL-Band, uniform) kommen aus params (BL-Sizing).
# Nur n_r3 (Buffer Kreis2->Quadrat) und n_z lokal:
n_r3 = max(1, round(t3 / cell_size))
n_z = max(1, int(math.ceil(dz / cell_size)))


def _ring_cells(verts, n_circ, n_rad):
    """(n1, n2, nz): Kante v0->v1 circumferentiell, wenn gleiche Schale."""
    v0, v1 = verts[0], verts[1]
    if ring_level[v0] == ring_level[v1]:
        return (n_circ, n_rad, n_z)
    return (n_rad, n_circ, n_z)


# TOPOLOGY=v2 "segments"-Tag pro Block: gibt an in welchen 3D-z-Segmenten
# (Swirl, Contraction, Orifice, Expansion, Exit, Ambient1, Ambient2) dieser
# 2D-Block existiert. sketch3d.py liest das und ueberspringt Bloecke ausserhalb
# ihrer aktiven Segmente. Bei TOPOLOGY=v1 (default) wird "ALL" gesetzt -> alte
# uniform Topologie.
SEG_ALL          = "ALL"          # alle z-Segmente (Inner-Schale)
SEG_CHAMBER      = "CHAMBER"      # nur Drallkammer + Kontraktion
SEG_EXIT_AMBIENT = "EXIT_AMBIENT" # ab Expansion bis amb_far
SEG_AMBIENT      = "AMBIENT"      # nur Ambient1 + Ambient2 (Ring4-Frustum)


def _segtag_for(name):
    if TOPOLOGY != "v2":
        return SEG_ALL
    if name.startswith("Trans_") or name.startswith("Out_"):
        return SEG_CHAMBER
    if name.startswith("Ring3_"):
        return SEG_EXIT_AMBIENT
    if name.startswith("Ring4_"):
        return SEG_AMBIENT
    return SEG_ALL


blocks = []

# Kern
blocks.append({"name": "Core", "verts": _ccw(("p0", "p1", "p2", "p3")),
               "cells": (N_core, N_core, n_z), "grading": (1.0, 1.0, 1.0),
               "segments": _segtag_for("Core")})

SIDE_NAMES = ["S", "E", "N", "W"]


def _add_ring(inner_base, outer_base, label, n_rad):
    for i in range(4):
        ia = f"p{inner_base + i}"
        ib = f"p{inner_base + (i + 1) % 4}"
        oa = f"p{outer_base + i}"
        ob = f"p{outer_base + (i + 1) % 4}"
        verts = _ccw((ia, ib, ob, oa))
        nm = f"{label}_{SIDE_NAMES[i]}"
        blocks.append({"name": nm, "verts": verts,
                       "cells": _ring_cells(verts, N_core, n_rad),
                       "grading": (1.0, 1.0, 1.0),
                       "segments": _segtag_for(nm)})


_add_ring(0, 4, "Ring1", n_r1)    # Quadrat -> Kreis 1
_add_ring(4, 8, "Ring2", n_r2)    # Kreis 1 -> Kreis 2
_add_ring(8, 12, "Trans", n_r3)   # Kreis 2 -> zentrales Quadrat (TOPOLOGY v2: nur Chamber)

# Aussen-Quadrate (3x3-Raster), CCW ab links-unten -- TOPOLOGY v2 nur Chamber
OUTER_BLOCKS = [
    ("Out_SW", ("p16", "p17", "p12", "p27")),
    ("Out_S",  ("p17", "p18", "p13", "p12")),
    ("Out_SE", ("p18", "p19", "p20", "p13")),
    ("Out_E",  ("p13", "p20", "p21", "p14")),
    ("Out_NE", ("p14", "p21", "p22", "p23")),
    ("Out_N",  ("p15", "p14", "p23", "p24")),
    ("Out_NW", ("p26", "p15", "p24", "p25")),
    ("Out_W",  ("p27", "p12", "p15", "p26")),
]
for name, verts in OUTER_BLOCKS:
    blocks.append({"name": name, "verts": _ccw(verts),
                   "cells": (N_core, N_core, n_z), "grading": (1.0, 1.0, 1.0),
                   "segments": _segtag_for(name)})

# TOPOLOGY=v2: 1 zusaetzlicher Frustum-Ring (Kreis2 -> Ring3)
if TOPOLOGY == "v2" and _R3_max is not None:
    # n_r5: Zellen radial im Ring3-Block (zwischen Kreis 2 und Ring3-Position).
    # Ring3-Block hat im 2D radiale Dicke (R3_max - R2). Im 3D wird Ring3-
    # Aussenkante per Frustum-scale auf R3(z) reduziert (zwischen R3 am
    # exit_top und R3_max am amb_far).
    # 16 cells radial -> @amb_far Cell-Groesse (R3_max-R2)/16 ~3.5mm, am
    # exit_top scale-reduziert auf (R3-R2)/16 ~0.5mm -- fein an der Lamelle.
    n_r5 = 16
    _add_ring(8, 28, "Ring3", n_r5)    # Kreis 2 -> Ring3 (Frustum)


# ----------------------------------------------------------------------
# 6) Patches
# ----------------------------------------------------------------------
# Aussen-Perimeter des 3x3-Rasters = Domain-Rand (provisorisch wall).
# Der Kreis (Schale 2) ist jetzt INTERNE Schnittstelle, kein Patch mehr.
patches = {
    "wall outerWall": [
        ("p16", "p17"), ("p17", "p18"), ("p18", "p19"),   # unten
        ("p19", "p20"), ("p20", "p21"), ("p21", "p22"),   # rechts
        ("p22", "p23"), ("p23", "p24"), ("p24", "p25"),   # oben
        ("p25", "p26"), ("p26", "p27"), ("p27", "p16"),   # links
    ],
}


# ----------------------------------------------------------------------
# 7) Skizze
# ----------------------------------------------------------------------
def _draw_edge(ax, a_name, b_name, **kw):
    mid = arc_lookup.get((a_name, b_name))
    if mid is not None:
        s = arc_samples(pts[a_name], pts[mid], pts[b_name], n=60)
        ax.plot([p[0] for p in s], [p[1] for p in s], **kw)
    else:
        xa, ya = pts[a_name]
        xb, yb = pts[b_name]
        ax.plot([xa, xb], [ya, yb], **kw)


def draw(fname="sketch.png", title="A3 SwirlNozzle O-Grid + Aussengitter (2D)",
         lim=None, label_fs=6.8):
    fig, ax = plt.subplots(figsize=(8.5, 8.5))

    # Alle Block-Kanten zeichnen (arc oder gerade); dedup ueber Kanten-Set
    seen = set()
    for blk in blocks:
        v = blk["verts"]
        for i in range(4):
            e = (v[i], v[(i + 1) % 4])
            key = tuple(sorted(e))
            if key in seen:
                continue
            seen.add(key)
            is_arc = (e in arc_lookup)
            _draw_edge(ax, e[0], e[1],
                       color="tab:blue" if is_arc else "0.55",
                       linewidth=1.8 if is_arc else 1.1, zorder=2)

    _lim = lim if lim is not None else 3 * c * 1.12

    # Aussen-Wand farbig
    for (a_name, b_name) in patches["wall outerWall"]:
        _draw_edge(ax, a_name, b_name, color="tab:red", linewidth=3.5,
                   zorder=3, solid_capstyle="round")
    if _lim >= 3 * c:
        ax.text(0, 3 * c * 1.03, "wall: outerWall", ha="center", va="bottom",
                color="tab:red", fontsize=12, fontweight="bold")

    # Block-Labels
    seg_colors = {"Core": "#332288", "Ring1": "#117733", "Ring2": "#882255",
                  "Trans": "#ddaa33", "Out": "#44aa99"}
    for blk in blocks:
        xs = [pts[v][0] for v in blk["verts"]]
        ys = [pts[v][1] for v in blk["verts"]]
        cx, cy = sum(xs) / 4.0, sum(ys) / 4.0
        if abs(cx) > _lim or abs(cy) > _lim:
            continue
        seg = blk["name"].split("_")[0]
        col = seg_colors.get(seg, "0.2")
        ax.text(cx, cy, f"{blk['name']}\n{blk['cells']}", ha="center",
                va="center", fontsize=label_fs, color=col, zorder=5,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=col,
                          lw=0.7, alpha=0.85))

    # Hauptpunkte
    for name, (x, y) in pts.items():
        if name.startswith("p_"):
            continue
        if abs(x) > _lim or abs(y) > _lim:
            continue
        ax.plot(x, y, "ko", markersize=4.5, zorder=6)
        ax.annotate(name, (x, y), textcoords="offset points",
                    xytext=(5, 5), fontsize=8.5, zorder=7)

    ax.set_aspect("equal")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title(title)
    if lim is None:
        lim = 3 * c * 1.12
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    fig.tight_layout()
    fig.savefig(fname, dpi=120)
    print(f"geschrieben: {fname}")


# ----------------------------------------------------------------------
# 8) Konsolen-Report
# ----------------------------------------------------------------------
if __name__ == "__main__":
    print(f"a  = {a} mm   (Kern-Quadrat {2*a} x {2*a} mm)")
    print(f"R0 (Ecke) = {R0:.4f} mm")
    print(f"R1 = {R1:.4f} mm  (t1 = {t1})")
    print(f"R2 = {R2:.4f} mm  (t2 = {t2})")
    print(f"c  = {c:.4f} mm  (zentr. Quadrat-Halbseite; t3 = {t3} an Kantenmitte)")
    print(f"   Eckspalt Kreis2->Quadrat = c*sqrt2 - R2 = {c*math.sqrt(2)-R2:.4f} mm")
    print(f"   Aussen-Domain = {6*c:.3f} x {6*c:.3f} mm  (3x3-Raster, Kachel 2c={2*c:.3f})")
    print(f"cell_size = {cell_size:.4f} mm   N_core = {N_core}")
    print(f"n_r1={n_r1} n_r2={n_r2} n_r3={n_r3} n_z={n_z} (dz={dz} Platzhalter)")
    print(f"-> {len(blocks)} Bloecke:")
    for blk in blocks:
        print(f"     {blk['name']:9s} verts={blk['verts']} cells={blk['cells']}")
    total_2d = sum(b["cells"][0] * b["cells"][1] for b in blocks)
    print(f"-> Zellen je z-Lage: {total_2d}  (x n_z = {total_2d*n_z} gesamt)")
    print(f"-> arcs: {len(arcs_edges)}  | outerWall: "
          f"{len(patches['wall outerWall'])} Kanten")
    draw(title=f"A3 SwirlNozzle O-Grid + 3x3-Aussengitter -- "
               f"a={a}, t1={t1}, t2={t2}, t3={t3}, N_core={N_core}")
    draw(fname="sketch_zoom.png", lim=c * 1.35, label_fs=8.5,
         title=f"A3 SwirlNozzle -- Zoom Kern-O-Grid (Core/Ring1/Ring2/Trans), "
               f"R1={R1:.3f} R2={R2:.3f} c={c:.3f} mm")
