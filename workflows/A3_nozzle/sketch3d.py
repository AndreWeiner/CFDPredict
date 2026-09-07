# -*- coding: utf-8 -*-
"""
A3 SwirlNozzleInjector -- 3D-Extrusion (koerperangepasst, Variante B).

Stapelt die 2D-O-Grid-Topologie (sketch.py) ueber die z-Stationen aus
params.py und SKALIERT jeden Querschnitt mit dem Bohrungsradius:

    Station            z            scale s = R(z)/(Do/2)
    swirl_bot   z_swirl_bot         Ds/Do  (=2)
    contr_bot   z_contr_bot         Ds/Do  (=2)
    orif_bot    0                   1
    orif_top    Lo                  1
    exp_top     Lo+Lk               R_exit/(Do/2)  (=2)
    exit_top    Lo+Lk+L_exit        R_exit/(Do/2)  (=2)

-> 21 2D-Bloecke x 5 z-Segmente = 105 Hex-Bloecke.
Kontraktion (contr_bot->orif_bot) und Expansion (orif_top->exp_top) sind
Frusta (lineare Skalenaenderung -> schraege Seitenwaende).

Patches (provisorisch, snappy verfeinert spaeter):
    outlet         Top-Flaeche (Exit, +z)
    chamberBottom  Boden Drallkammer (-z)  [wall; Zulauf via snappy]
    outerWall      Aussen-Perimeter (3x3-Raster), alle Segmente  [wall]
Innerer Kreis = interne Schnittstelle.

Einheiten: mm -> scale 0.001 im blockMeshDict.
"""

import math
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import params as P
import sketch as S   # 2D-Topologie: S.pts, S.blocks, S.arcs_edges, S.patches


# ----------------------------------------------------------------------
# 1) 2D-Basis
# ----------------------------------------------------------------------
MAIN = sorted((k for k in S.pts if not k.startswith("p_")), key=lambda k: int(k[1:]))
N2D = len(MAIN)                       # 28 Hauptpunkte
idx2d = {k: i for i, k in enumerate(MAIN)}
xy = {k: S.pts[k] for k in MAIN}      # 2D-Koordinaten (Orifice-Skala, s=1)

# ----------------------------------------------------------------------
# 2) z-Stationen (z, scale)
# ----------------------------------------------------------------------
R_orif = P.R_wall
# Chamber O-grid scale. Body-fitted would be Ds/Do = 2.0 (Kreis 2 lands exactly
# on the chamber wall Rs -> the tangential inlets then cut through the curved
# arc + BL band, which snaps to slivers/high skew). Reducing it keeps the inner
# O-grid (incl. p12-15) ~inside Rs so the chamber wall + inlets are carved by
# snappy through the SQUARE outer raster (Cartesian) instead -> much less skew.
# 1.4: Kreis 2 ~3.5mm, Trans corners ~5.6mm (just past Rs=5); fully-inside would
# need ~1.25, the trade keeps more chamber O-grid resolution.
S_CHAMBER = 1.4
# PRISMATIC: constant cross-section over z (all scales = 1.0) -> the O-grid
# "circle segments go straight up" with no frusta. snappy then carves the whole
# nozzle (chamber Rs, contraction, orifice Rw, expansion, exit) out of the
# uniform prism. The orifice still matches Kreis 2 (=Rw) arc-exact for the
# air-core/film; everything wider is cut from the square raster. Uniform cells
# -> cleaner refinement transitions (no tiny cells from frusta + refinement).
# False -> the body-fitted variant (S_CHAMBER chamber, R_exit expansion).
PRISMATIC = True
_HAS_EXIT = P.L_exit > 1e-6           # L_exit=0 -> kein Exit-Geradstueck (paper-konform)

if PRISMATIC:
    # Chamber + contraction + orifice prismatic (scale 1, circle segments straight)
    # -> uniform cells, clean tangential carve. Expansion/exit KEEP the body-fitted
    # widening (scale R_exit/Do/2): the core blocks widen toward the outlet -- this
    # gave a nice spray-exit region, not made prismatic there.
    _se = P.R_exit / R_orif
    STATIONS = [
        ("swirl_bot", P.z_swirl_bot, 1.0),
        ("contr_bot", P.z_contr_bot, 1.0),
        ("orif_bot",  0.0,           1.0),
        ("orif_top",  P.z_orif_top,  1.0),
        ("exp_top",   P.z_exp_top,   _se),
    ]
    if _HAS_EXIT:
        STATIONS.append(("exit_top", P.z_exit_top, _se))
else:
    _se = P.R_exit / R_orif
    STATIONS = [
        ("swirl_bot", P.z_swirl_bot, S_CHAMBER),
        ("contr_bot", P.z_contr_bot, S_CHAMBER),
        ("orif_bot",  0.0,           1.0),
        ("orif_top",  P.z_orif_top,  1.0),
        ("exp_top",   P.z_exp_top,   _se),
    ]
    if _HAS_EXIT:
        STATIONS.append(("exit_top", P.z_exit_top, _se))
CASE_DIR = "case_a3_prism" if PRISMATIC else "case_a3"
SEG_NAMES = ["Swirl", "Contraction", "Orifice", "Expansion"]
if _HAS_EXIT:
    SEG_NAMES.append("Exit")
N_NOZZLE_SEG = len(SEG_NAMES)        # 4 oder 5 nozzle segments before ambient

# AMBIENT extension (Task 3): zwei zusaetzliche z-Stationen unter dem Nozzle-Exit
# fuer die Sprueh-Kegelwinkel-Auswertung. Die O-Grid-Frusta oeffnen sich frustum-
# artig in den Spray-Konus (Scale linear vom Exit-Scale zu R_amb/R_orif). Das
# Aussen-Raster dieser Segmente + die Top-Flaeche der amb_far-Station landen im
# 'outlet'-Patch (atmosphaerische pressureInletOutletVelocity-BC, identisch zur
# Nozzle-Exit-BC) -> case_template/0/* bleibt unveraendert.
if P.AMBIENT:
    _s_exit  = STATIONS[-1][2]
    if getattr(P, "TOPOLOGY", "v1") == "v2":
        # v2: Inner-O-Grid faechert mit demselben Cone-Half-Angle wie Ring3 auf
        # -- damit folgt der Inner-O-Grid (Core+Ring1+Ring2) der Spray-Lamelle.
        # Ring3 bleibt dabei radial konstant in der Dicke (R3-R_exit = 5mm).
        # Vorher prismatisch: am chamberBottom-Uebergang R_inner=R_exit=5mm fix,
        # die Lamelle musste vom feinen Inner-O-Grid in den groben Ring3 umsteigen
        # -> Druck-Spike (Martin's ParaView 6.13ms). Jetzt waechst Inner-O-Grid
        # mit der Lamelle mit -> Druck-Spike loest besser auf.
        # SHRINK 0.9 (Martin's Vorschlag 2026-06-02): an den Ambient-Boden-
        # Stationen (amb_mid, amb_far) den Inner-O-Grid-Aussenrand leicht nach
        # innen versetzen (~10%), damit Ring3 mehr Pufferzone bekommt und der
        # spitze Block-Vertex am Frustum-Knick entschaerft wird (siehe Skizze
        # mit orange/gruen Pfeil).
        _half = math.radians(getattr(P, "cone_half_angle_deg", 60.0))
        _tan_h = math.tan(_half)
        _shrink = getattr(P, "inner_amb_shrink", 0.9)
        _R_inner_m = (P.R_exit + (P.z_amb_mid - P.z_exit_top) * _tan_h) * _shrink
        _R_inner_f = (P.R_exit + (P.z_amb_far - P.z_exit_top) * _tan_h) * _shrink
        _s_amb_m = _R_inner_m / P.R_wall
        _s_amb_f = _R_inner_f / P.R_wall
    else:
        _s_amb_f = P.R_amb / R_orif
        _s_amb_m = _s_exit + (_s_amb_f - _s_exit) * P.amb_mid_frac_s
    STATIONS.append(("amb_mid", P.z_amb_mid, _s_amb_m))
    STATIONS.append(("amb_far", P.z_amb_far, _s_amb_f))
    SEG_NAMES = SEG_NAMES + ["Ambient1", "Ambient2"]
NST = len(STATIONS)

# axiale Aufloesung. Orifice und Expansion bekommen 1.5x mehr z-Cells als die
# default h_ax-Berechnung -- das lange Rohr und der sich aufweitende Outflow
# brauchen feinere axiale Aufloesung, damit die Wasser/Luft-Front nicht ueber
# zu wenige Cells streckt (Spray-Detonation @ 7.97ms ohne Refinement).
axial_aspect = 4.0
h_ax = P.h_base * axial_aspect
_z_mult = {"Orifice": 1.5, "Expansion": 2.33, "Ambient1": 2.0, "Contraction": 2.0}
nz_seg = [max(1, round(_z_mult.get(SEG_NAMES[j], 1.0) *
                       abs(STATIONS[j + 1][1] - STATIONS[j][1]) / h_ax))
          for j in range(NST - 1)]


# ----------------------------------------------------------------------
# 3) 3D-Vertices:  V[j*N2D + i]
# ----------------------------------------------------------------------
def vid(station_j, pname):
    return station_j * N2D + idx2d[pname]


# FLIP: 180 deg about the Y axis -> (x,y,z) -> (-x, y, -z). Turns the nozzle so
# the exit faces -z (fluid sprays in -z, chamber on top), matching the CFDPredict
# convention. A proper rotation (det +1) -> hex volumes stay positive, no winding
# fix. (x is mirrored too -> swirl handedness flips; cosmetic for symmetric BCs.)
def _flip(x, y, z):
    return (-x, y, -z)


# TOPOLOGY=v2: ring_level=5 (Ring3) wird Frustum-skaliert:
#   - 2D-Position pts[] in sketch.py setzt Ring3 corners auf R3_max
#     (= R3 + L_amb*tan(half_angle)).
#   - Pro Station wird Ring3-Radius linear interpoliert:
#       R3(z) = R3                       fuer z <= z_exit_top
#               R3 + (z-z_exit_top)*tan(half_angle)  fuer z > z_exit_top (Ambient)
#   - sketch3d-scale fuer Ring3 = R3(z) / R3_max
TOPOLOGY_V2 = (getattr(S, "TOPOLOGY", "v1") == "v2")
_R3_ABS  = getattr(S, "_R3_abs", None) or 0.0
_R3_MAX  = getattr(S, "_R3_max", None)
_TAN_HC  = math.tan(math.radians(getattr(P, "cone_half_angle_deg", 60.0)))


def _vertex_scale(k, j, station_scale, station_z_orig):
    """Per-vertex scale: ring_level <=4 = uniform station_scale (wie v1).
    Ring3 (level 5): Frustum mit half-angle. Im 2D liegt der pt auf R3_max
    (= R3 + L_amb*tan(half_angle)). Pro Station: R3(z) = R3 fuer z <= z_exit_top,
    sonst R3 + (z-z_exit_top)*tan. scale = R3(z) / R3_max."""
    if not TOPOLOGY_V2:
        return station_scale
    lvl = S.ring_level.get(k, 0)
    if lvl <= 4:
        return station_scale
    if lvl == 5 and _R3_MAX is not None:
        # Ring3 Frustum: R3 am exit_top, waechst downstream.
        delta = max(0.0, station_z_orig - P.z_exit_top)
        R3_at = _R3_ABS + delta * _TAN_HC
        return R3_at / _R3_MAX
    return station_scale


verts3d = []
for j, (_nm, z, s) in enumerate(STATIONS):
    for k in MAIN:
        x, y = xy[k]
        s_v = _vertex_scale(k, j, s, z)
        verts3d.append(_flip(x * s_v, y * s_v, z))


# ----------------------------------------------------------------------
# 4) 3D-Hex-Bloecke:  je 2D-Block x je z-Segment
# ----------------------------------------------------------------------
def hex_orient(v8):
    """Spatprodukt an v0:  ((v1-v0) x (v3-v0)) . (v4-v0).
    > 0  => blockMesh-konforme Orientierung (positives Zellvolumen).
    Betrag ~ lokales Volumenmasz."""
    p = [np.array(verts3d[i], dtype=float) for i in v8]
    return float(np.dot(np.cross(p[1] - p[0], p[3] - p[0]), p[4] - p[0]))


# TOPOLOGY=v2: Welche Segments existieren fuer einen Block-segments-Tag?
# Bei v1 ist segments="ALL" -> immer aktiv, kompatibel mit alter Logik.
def _seg_active(seg_tag, seg_name):
    if seg_tag == "ALL":
        return True
    if seg_tag == "CHAMBER":
        return seg_name in ("Swirl", "Contraction")
    if seg_tag == "EXIT_AMBIENT":
        # ab Expansion + Exit (falls vorhanden) + Ambient1/2
        return seg_name in ("Expansion", "Exit", "Ambient1", "Ambient2")
    if seg_tag == "AMBIENT":
        return seg_name in ("Ambient1", "Ambient2")
    return False


blocks3d = []
for j in range(NST - 1):
    seg_name = SEG_NAMES[j]
    for b in S.blocks:
        seg_tag = b.get("segments", "ALL")
        if not _seg_active(seg_tag, seg_name):
            continue
        v0, v1, v2, v3 = b["verts"]
        v8 = [vid(j, v0), vid(j, v1), vid(j, v2), vid(j, v3),
              vid(j + 1, v0), vid(j + 1, v1), vid(j + 1, v2), vid(j + 1, v3)]
        n1, n2, _ = b["cells"]
        blocks3d.append({
            "name": f"{b['name']}_{seg_name}",
            "v8": v8,
            "cells": (n1, n2, nz_seg[j]),
            "seg": seg_name,
        })


# ----------------------------------------------------------------------
# 5) Arc-Kanten je Station (skalierter Kreis 1 / Kreis 2)
# ----------------------------------------------------------------------
# TOPOLOGY=v2: arcs fuer Ring3 (p28..p31) nur an Stations wo Ring3-Block
# tatsaechlich existiert (= EXIT_AMBIENT segments). Sonst meldet blockMesh
# "Curved edge does not correspond to a block edge" Warnings und einen
# fatalen Topology-Konsistenz-Check (Patch-Face-Edge passt nicht zur
# Block-Face-Edge weil arc den geraden Block-Edge biegt).
def _arc_active_at_station(midkey, j):
    if not TOPOLOGY_V2:
        return True
    # Ring3 arcs (p28..p31): nur in EXIT_AMBIENT-Segmenten ODER an Stations
    # die diese Segmente begrenzen.
    if "28_" in midkey or "29_" in midkey or "30_" in midkey or "31_" in midkey:
        seg_below = SEG_NAMES[j - 1] if j > 0 else None
        seg_above = SEG_NAMES[j] if j < NST - 1 else None
        exit_amb = ("Expansion", "Exit", "Ambient1", "Ambient2")
        return (seg_below in exit_amb) or (seg_above in exit_amb)
    return True


arcs3d = []   # (i_start, i_end, (mx,my,mz))
for j, (_nm, z, s) in enumerate(STATIONS):
    for (sname, midkey, ename) in S.arcs_edges:
        if not _arc_active_at_station(midkey, j):
            continue
        mx, my = S.pts[midkey]
        # ring_level vom Endpunkt nutzen (arcs liegen auf einer Schale).
        s_v = _vertex_scale(sname, j, s, z)
        arcs3d.append((vid(j, sname), vid(j, ename), _flip(mx * s_v, my * s_v, z)))


# ----------------------------------------------------------------------
# 6) Patches
# ----------------------------------------------------------------------
# outlet: Top-Flaeche aller 2D-Bloecke an der LETZTEN Station -- nur Bloecke
# deren segments-tag das letzte 3D-Segment enthaelt (in v2 ist das "Ambient2";
# Trans/Out gibt's dort nicht).
outlet_faces = []
jlast = NST - 1
last_seg = SEG_NAMES[jlast - 1]
for b in S.blocks:
    if not _seg_active(b.get("segments", "ALL"), last_seg):
        continue
    v0, v1, v2, v3 = b["verts"]
    outlet_faces.append((vid(jlast, v0), vid(jlast, v1),
                         vid(jlast, v2), vid(jlast, v3)))

# chamberBottom: Boden-Flaeche an Station 0 (Normale -z -> umgekehrte Reihenfolge)
# Nur Bloecke deren segments-tag das erste 3D-Segment enthaelt (in v2 "Swirl").
bottom_faces = []
first_seg = SEG_NAMES[0]
for b in S.blocks:
    if not _seg_active(b.get("segments", "ALL"), first_seg):
        continue
    v0, v1, v2, v3 = b["verts"]
    bottom_faces.append((vid(0, v0), vid(0, v3), vid(0, v2), vid(0, v1)))

# outerWall: Aussen-Perimeter, **NUR Nozzle-Segmente** (j < N_NOZZLE_SEG).
# Mit AMBIENT: die Ambient-Segmente-Aussenseiten werden weiter unten in den
# outlet-Patch gemerged (atmosphaerische BC).
#
# TOPOLOGY=v2: Out_*-Bloecke existieren nur in Swirl+Contraction. Im Orifice/
# Expansion/Exit/Ambient gibt es keinen Outer-Quadrat-Raster. Stattdessen:
#   - In Expansion/Ambient: Ring3-Block's Aussenkante (Ring3-Kreis) ist
#     der "ehemalige outerWall" -- bei AMBIENT=True wird der von snappy
#     ggf. weggeschnitten / gehoert sonst zum atmosphere/outlet-Patch.
#   - In Ambient1+2: Ring4-Block's Aussenkante = Frustum-Mantel (atmosphere).
outer_edges = S.patches["wall outerWall"]
outerwall_faces = []                # nozzle: wall (snappy schneidet hier)
atmosphere_faces = []               # ambient: patch (-> outlet)
# L_exit Aussenrand-Faces. Bei P.AMBIENT=False ist das L_exit-Geradstueck eine
# artificialle Wand-Verlaengerung der Duese (in Zhang's Setup nicht vorhanden);
# Coanda-Zirkulation an dieser Wand triggerte in v25 den limitVelocity dauerhaft
# (z=-25mm Spike-Region). Ohne Ambient bekommt der Mantel jetzt einen eigenen
# atmosphaerischen Patch. Mit Ambient bleibt es wall (echte Material-Aussenwand
# der Duese).
lexit_side_faces = []

if TOPOLOGY_V2:
    # v2: outer_edges (p16-p27) sind nur in Swirl+Contraction-Bloecken referenziert
    chamber_segs = ("Swirl", "Contraction")
    for j in range(NST - 1):
        if SEG_NAMES[j] not in chamber_segs:
            continue
        for (a, b) in outer_edges:
            face = (vid(j, a), vid(j, b), vid(j + 1, b), vid(j + 1, a))
            outerwall_faces.append(face)
    # Ring3-Aussenkante = Frustum-Mantel des einzigen v2-Aussenrings.
    # In Expansion (Nozzle): wall-aequivalent (snappy schneidet hier den
    # echten Nozzle-Aussenrand). In Ambient1/2: atmosphere -> outlet.
    ring3_edges = [(f"p{28+i}", f"p{28+(i+1)%4}") for i in range(4)]
    # Ring2-Aussenkante (p8..p11): im ORIFICE ist Ring2 die outermost shell
    # (kein Ring3 dort) -> outerwall-patch noetig. In Expansion+Ambient ist
    # Ring2-Aussenkante SHARED mit Ring3-Innenkante (internal face).
    ring2_edges = [(f"p{8+i}", f"p{8+(i+1)%4}") for i in range(4)]
    for j in range(NST - 1):
        sn = SEG_NAMES[j]
        if sn == "Orifice":
            # Ring2-Aussenrand = Bohrungswand (snappy fittet hier)
            for (a, b) in ring2_edges:
                face = (vid(j, a), vid(j, b), vid(j + 1, b), vid(j + 1, a))
                outerwall_faces.append(face)
        elif sn == "Expansion":
            # Ring3-Aussenrand = Bohrungs-/Expansion-Wand-Aequivalent
            for (a, b) in ring3_edges:
                face = (vid(j, a), vid(j, b), vid(j + 1, b), vid(j + 1, a))
                outerwall_faces.append(face)
        elif sn == "Exit":
            # L_exit-Geradstueck: ohne Ambient = atmosphaerischer Patch
            # (sonst artificialle Wand, siehe lexit_side_faces-Kommentar).
            target = outerwall_faces if P.AMBIENT else lexit_side_faces
            for (a, b) in ring3_edges:
                face = (vid(j, a), vid(j, b), vid(j + 1, b), vid(j + 1, a))
                target.append(face)
        elif sn in ("Ambient1", "Ambient2"):
            for (a, b) in ring3_edges:
                face = (vid(j, a), vid(j, b), vid(j + 1, b), vid(j + 1, a))
                atmosphere_faces.append(face)

    # Ring3-Block-Bottom-Face an station orif_top (Expansion-Start): an dieser
    # Station gibt es keinen Ring3-Block darunter (Orifice-Segment hat kein
    # Ring3) -> Patch (chamberBottom-aequivalent).
    j_orif_top = None
    for j, (nm, _z, _s) in enumerate(STATIONS):
        if nm == "orif_top":
            j_orif_top = j
            break
    if j_orif_top is not None:
        for b in S.blocks:
            if not b["name"].startswith("Ring3_"):
                continue
            v0, v1, v2, v3 = b["verts"]
            bottom_faces.append((vid(j_orif_top, v0), vid(j_orif_top, v3),
                                 vid(j_orif_top, v2), vid(j_orif_top, v1)))
else:
    for j in range(NST - 1):
        sn = SEG_NAMES[j]
        is_nozzle = (j < N_NOZZLE_SEG)
        if sn == "Exit" and not P.AMBIENT:
            target = lexit_side_faces
        elif is_nozzle:
            target = outerwall_faces
        else:
            target = atmosphere_faces
        for (a, b) in outer_edges:
            face = (vid(j, a), vid(j, b), vid(j + 1, b), vid(j + 1, a))
            target.append(face)

# chamberBottom_top (TOPOLOGY=v2): die TOP-Faces der Trans/Out-Bloecke an
# Station "Orif_bot" (= Ende Contraction). Diese werden mit chamberBottom
# gemerged (gleicher noSlip-wall-Patch, funktional).
if TOPOLOGY_V2:
    j_orif_bot = None
    for j, (nm, _z, _s) in enumerate(STATIONS):
        if nm == "orif_bot":
            j_orif_bot = j
            break
    if j_orif_bot is not None:
        # Trans+Out blocks endeten an j_orif_bot-1 (Contraction-Segment).
        # Ihre TOP-Faces (= blockMesh face 5, order v0..v3 an Station j+1) an
        # Station j_orif_bot werden bottom_faces zugefuegt.
        for b in S.blocks:
            if b.get("segments") != "CHAMBER":
                continue
            v0, v1, v2, v3 = b["verts"]
            bottom_faces.append((vid(j_orif_bot, v0), vid(j_orif_bot, v1),
                                 vid(j_orif_bot, v2), vid(j_orif_bot, v3)))


# ----------------------------------------------------------------------
# 7) blockMeshDict-Export
# ----------------------------------------------------------------------
_HEADER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
|  \\\\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox           |
|   \\\\    /   O peration    | Version:  v2512                                 |
|    \\\\  /    A nd          | Web:      www.OpenFOAM.com                      |
|     \\\\/     M anipulation |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def write_blockmeshdict(case_path):
    sysdir = os.path.join(case_path, "system")
    os.makedirs(sysdir, exist_ok=True)
    L = []
    L.append("scale   0.001;   // mm -> m")
    L.append("")
    L.append("vertices")
    L.append("(")
    for n, (x, y, z) in enumerate(verts3d):
        st = n // N2D
        pn = MAIN[n % N2D]
        L.append(f"    ( {x:13.6e} {y:13.6e} {z:13.6e} )  // {n:3d} {STATIONS[st][0]}:{pn}")
    L.append(");")
    L.append("")
    L.append("blocks")
    L.append("(")
    # z-Grading 1.5 in Ambient1/Ambient2: feine Cells am exit_top (Lamellen-
    # Austritt am Lampenschirm-Disk-Knick), groebere weiter weg. Vorher
    # uniform -> langgestreckte Cells direkt unter dem Lampenschirm, die
    # Lamelle sah dort einen Geometrie-Knick + lange Cells (Martin's ParaView
    # 2026-06-02, v5 alpha=1 vs Mesh-Stretching).
    _z_grad = {"Ambient1": 4.0, "Ambient2": 2.0}
    for blk in blocks3d:
        vs = " ".join(str(i) for i in blk["v8"])
        nx, ny, nz = blk["cells"]
        gz = _z_grad.get(blk.get("seg", ""), 1.0)
        L.append(f"    hex ({vs}) {blk['name']} ({nx} {ny} {nz}) "
                 f"simpleGrading (1 1 {gz})")
    L.append(");")
    L.append("")
    L.append("edges")
    L.append("(")
    for (i, jj, (mx, my, mz)) in arcs3d:
        L.append(f"    arc {i} {jj} ({mx:13.6e} {my:13.6e} {mz:13.6e})")
    L.append(");")
    L.append("")
    L.append("boundary")
    L.append("(")

    def _patch(name, ptype, faces):
        L.append(f"    {name}")
        L.append("    {")
        L.append(f"        type {ptype};")
        L.append("        faces")
        L.append("        (")
        for f in faces:
            L.append("            (" + " ".join(str(i) for i in f) + ")")
        L.append("        );")
        L.append("    }")

    # outlet = nozzle-exit/ambient-far Top-Flaeche. Frustum-Mantel im Ambient
    # geht in einen separaten "atmosphere"-Patch (s.u.) -- vorher in outlet
    # gemerged, das gab BC-Probleme: outlet ist stiff fixedValue 0 fuer das
    # bottom (sauber non-reflecting), Frustum-Mantel braucht atmosphere-typ
    # BC (alpha=inletOutlet, p_rgh=prghPressure, U=pressureInletOutletVelocity)
    # damit Luft fuer Air-Core nachstroemen kann ohne dass das stiff p=0 die
    # Druck-Pulsation an der Aussenseite des Sprays clamped.
    # outlet =
    # Raster-Faces der Ambient-Segmente. Beide tragen die gleiche atmosphaerische
    # pressureInletOutletVelocity/totalPressure-BC.
    _patch("outlet", "patch", outlet_faces)
    if atmosphere_faces:
        _patch("atmosphere", "patch", atmosphere_faces)
    if lexit_side_faces:
        _patch("lexit_side", "patch", lexit_side_faces)
    _patch("chamberBottom", "wall", bottom_faces)
    _patch("outerWall", "wall", outerwall_faces)
    L.append(");")
    L.append("")
    L.append("mergePatchPairs ( );")
    body = "\n".join(L)
    path = os.path.join(sysdir, "blockMeshDict")
    with open(path, "w", newline="\n") as f:
        f.write(_HEADER.format(obj="blockMeshDict"))
        f.write("\n" + body + "\n// " + "*" * 70 + " //\n")
    print(f"geschrieben: {path}")

    # minimale controlDict, damit blockMesh/checkMesh direkt laufen
    cdict = (_HEADER.format(obj="controlDict") + "\n"
             "application     blockMesh;\n"
             "startFrom       startTime;\nstartTime       0;\n"
             "stopAt          endTime;\nendTime         1;\n"
             "deltaT          1;\nwriteControl    timeStep;\nwriteInterval   1;\n"
             "writeFormat     ascii;\nwritePrecision  7;\nrunTimeModifiable true;\n"
             "// " + "*" * 70 + " //\n")
    with open(os.path.join(sysdir, "controlDict"), "w", newline="\n") as f:
        f.write(cdict)

    fvsch = (_HEADER.format(obj="fvSchemes") + "\n"
             "ddtSchemes { default steadyState; }\n"
             "gradSchemes { default Gauss linear; }\n"
             "divSchemes { default none; }\n"
             "laplacianSchemes { default Gauss linear corrected; }\n"
             "interpolationSchemes { default linear; }\n"
             "snGradSchemes { default corrected; }\n"
             "// " + "*" * 70 + " //\n")
    with open(os.path.join(sysdir, "fvSchemes"), "w", newline="\n") as f:
        f.write(fvsch)
    fvsol = (_HEADER.format(obj="fvSolution") + "\n"
             "solvers {}\n"
             "// " + "*" * 70 + " //\n")
    with open(os.path.join(sysdir, "fvSolution"), "w", newline="\n") as f:
        f.write(fvsol)
    os.makedirs(os.path.join(case_path, "constant"), exist_ok=True)
    print(f"geschrieben: controlDict, fvSchemes, fvSolution")


# ----------------------------------------------------------------------
# 8) 3D-Preview (Silhouette: Aussen-Perimeter + Bohrungskreis je Station)
# ----------------------------------------------------------------------
def draw3d(fname="sketch3d.png"):
    fig = plt.figure(figsize=(8.5, 9.5))
    ax = fig.add_subplot(111, projection="3d")

    # Bohrungskreis (Kreis2 p8..p11) + Aussen-Perimeter je Station
    outer_names = ["p16", "p17", "p18", "p19", "p20", "p21", "p22", "p23",
                   "p24", "p25", "p26", "p27", "p16"]
    for j, (nm, z, s) in enumerate(STATIONS):
        # Bohrung (Kreis 2) als Polygon der 4 Bogenmittelpunkte + Ecken (gesampelt)
        th = np.linspace(0, 2 * np.pi, 81)
        ax.plot(P.R_wall * s * np.cos(th), P.R_wall * s * np.sin(th),
                zs=z, color="tab:blue", lw=1.6)
        # Aussen-Perimeter
        ox = [xy[n][0] * s for n in outer_names]
        oy = [xy[n][1] * s for n in outer_names]
        ax.plot(ox, oy, zs=z, color="tab:red", lw=1.6)

    # vertikale Kanten der Aussen-Ecken + Bohrungs-Eckpunkte
    for nm in ["p16", "p19", "p22", "p25", "p8", "p9", "p10", "p11"]:
        xs = [xy[nm][0] * s for (_n, _z, s) in STATIONS]
        ys = [xy[nm][1] * s for (_n, _z, s) in STATIONS]
        zs = [z for (_n, z, _s) in STATIONS]
        col = "tab:blue" if nm in ("p8", "p9", "p10", "p11") else "tab:red"
        ax.plot(xs, ys, zs, color=col, lw=1.0, alpha=0.7)

    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    ax.set_title("A3 SwirlNozzle 3D (koerperangepasst, B)\n"
                 "rot=Aussenrand, blau=Bohrungswand")
    ax.view_init(elev=14, azim=-60)
    try:
        ax.set_box_aspect((1, 1, 3))
    except Exception:
        pass
    fig.tight_layout()
    fig.savefig(fname, dpi=120)
    print(f"geschrieben: {fname}")


# ----------------------------------------------------------------------
# 9) Report + Sanity
# ----------------------------------------------------------------------
if __name__ == "__main__":
    print(f"Stationen: {[(n, round(z,2), round(s,3)) for n,z,s in STATIONS]}")
    print(f"nz je Segment {dict(zip(SEG_NAMES, nz_seg))}  (h_ax={h_ax:.3f} mm)")
    print(f"Vertices: {len(verts3d)}   Bloecke: {len(blocks3d)}   Arcs: {len(arcs3d)}")
    total = sum(b['cells'][0] * b['cells'][1] * b['cells'][2] for b in blocks3d)
    print(f"Zellen gesamt: {total}")
    # Orientierungs-Check (Spatprodukt > 0 = blockMesh-konform)
    orients = [hex_orient(b["v8"]) for b in blocks3d]
    bad = [blocks3d[i]["name"] for i, v in enumerate(orients) if v <= 0]
    print(f"min Spatprodukt={min(orients):.4e}  invertiert/0: {len(bad)}")
    if bad:
        print("  !! invertierte Bloecke:", bad[:8])
    print(f"Patches: outlet={len(outlet_faces) + len(atmosphere_faces)} "
          f"(top={len(outlet_faces)} +ambientSide={len(atmosphere_faces)})  "
          f"chamberBottom={len(bottom_faces)}  outerWall={len(outerwall_faces)}  "
          f"lexit_side={len(lexit_side_faces)}  AMBIENT={P.AMBIENT}")
    draw3d()
    write_blockmeshdict(CASE_DIR)
    print(f"case dir: {CASE_DIR}  (PRISMATIC={PRISMATIC})")
