# -*- coding: utf-8 -*-
"""A3 swirl nozzle -- snappyHexMesh layer on top of the body-fitted O-grid.

The blockMesh O-grid (sketch3d.py -> case_a3) is the structured background: the
bore (chamber -> contraction -> orifice -> expansion -> exit) already lies
arc-exact on the O-grid circles, wrapped in the 3x3 outer raster (collar) at
every z. snappyHexMesh has to:

    1. delete the collar everywhere (everything outside the bore wall), keeping
       the bore interior + the 4 tangential inlet tubes;
    2. open the chamber wall at each tangential penetration;
    3. snap the bore wall (level 0, minimal -- it already matches the O-grid).

Surfaces (metres; the blockMesh uses scale 0.001 so the mesh is in metres):
    nozzle.stl      ONE union(bore body, 4 inlets) with the top (outlet) and
                    bottom (chamberBottom) caps removed so the bore stays open
                    to those blockMesh patches. 2 named solids:
                      wall   bore lateral wall + inlet feed pipe laterals
                      inlet  the 4 outer feed end-caps (velocity-inlet patch)
    feat_inlet.stl  watertight union(chamber cylinder, inlets) -> the
                    tangential penetration curves as feature edges (the bore
                    junctions are NOT in this surface, so the body-fitted O-grid
                    frusta stay untouched -- minimal snappy work, by design).

Run from this directory (imports params.py):
    python make_snap.py            # writes case_a3/constant/triSurface + system dicts
"""
import math
import os
import sys

import numpy as np
import trimesh

import params as P

_ENGINE = "manifold"
M = 0.001                      # mm -> m (blockMesh scale)
SECTIONS = 180                 # azimuthal facets of the bore body of revolution

# --- Top-Seal (nur mit AMBIENT) -------------------------------------------
# Annulare duenne Scheibe direkt stromab des Nozzle-Exits (im original Frame:
# z >= z_exit_top; im geflippten Frame: z <= -z_exit_top). Schliesst den Aussen-
# Raster-Collar (R zwischen R_exit und R_outer*scale_exit) gegen die Ambient-Cells
# ab, damit snappy den Collar korrekt als "ausserhalb der Bore" identifiziert und
# loescht. Ohne diese Scheibe verbinden die Ambient-Cells den Collar lateral zur
# Bore-Exit-Oeffnung, snappy haelt die Collar-Cells dann faelschlich.
# Patch: gemerged in den existierenden `wall`-Patch (gleiche noSlip-Physik wie
# die Bore-Innenwand) -> keine Aenderung im case_template/0/* noetig.
# WICHTIG: das blockMesh-Outer-Raster ist QUADRATISCH (3x3-Raster). Der Disk-
# Aussenradius muss die DIAGONALE der Quadrat-Ecken abdecken (= sqrt(2)*R_outer
# *scale), sonst leaken die 4 Eckspalten den Collar zum Ambient. Verifiziert
# Martin am 2026-06-01.
import math as _math
TOP_SEAL_THICK_MM    = 0.5      # axiale Dicke der Scheibe
TOP_SEAL_MARGIN_MM   = 2.0      # wie weit ueber die Quadrat-Diagonal-Ecke hinaus
TOP_SEAL_INNER_OFFSET_MM = 0.1  # R_in = R_exit + offset; sub-cell radialer Abstand
                                # zur Bore-Wand-STL, damit der Disk-Inner-Cylinder
                                # NICHT an der Bore-Exit-Edge (R=R_exit, z=z_exit_top)
                                # koinzidiert -> snappy snap erzeugt sonst Slivers
                                # (Min volume ~1e-13 m^3 -> FPE in RNG-k-epsilon).

# ---- chamber refinement strategy ------------------------------------------
# PREREFINE_LEVEL == 0  -> DEFAULT / validated endform (variant #6): 4 searchableBox
#   refinementRegions tile the chamber outer raster to a uniform absolute level 2.
#   The box extents are COMPUTED from the chamber topology (P.c .. P.Rs+2*Di, z
#   swirl_bot..orif_bot), not hand-tuned -> snappy-native, owns its hexRef8 history,
#   checkMesh skew 2.58 / 0 illegal.
# PREREFINE_LEVEL > 0  -> INVESTIGATED, NOT recommended. Pre-refine the named
#   chamber outer-raster cellZones (Out_*_Swirl/_Contraction) via topoSet
#   (zoneToCell) + refineMesh (the Allmesh runs it PREREFINE_LEVEL times). Elegant
#   on paper (topology-exact, cannot bleed into the inner O-grid) but does NOT
#   compose with snappy: standalone refineMesh writes no cell/point-level history,
#   so refineMesh->snappy crashes hexRef8 at L=2 (2-level jump -> danglingCellRefine)
#   and only degrades to skew 5.44 at L=1; snappy->refineMesh refines snap-distorted
#   cells -> skew 21.7. Kept behind the flag so the negative result is reproducible.
PREREFINE_LEVEL = 0
# cellZones to pre-refine: outer square raster over the chamber + contraction
# (where the inlets sit). Orifice/Expansion/Exit Out_* are pure collar snappy
# deletes -> not listed. Inner O-grid blocks are never listed.
PREREFINE_SEGMENTS = ("Swirl", "Contraction")

# FLIP: 180 deg about Y -> (x,y,z) -> (-x, y, -z). Matches sketch3d.py so the
# nozzle exit faces -z (fluid sprays in -z). Proper rotation -> STL stays valid
# and in the same frame as the flipped blockMesh background.
_R_FLIP = np.array([[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]], float)


def _flip_pt(p):
    p = np.asarray(p, float)
    return np.array([-p[0], p[1], -p[2]])


# ----------------------------------------------------------------------
# geometry (metres)
# ----------------------------------------------------------------------
def _bore_solid():
    """Closed body of revolution of the bore (chamber..exit), axis +z, metres."""
    Rs, Rw, Re = P.Rs, P.R_wall, P.R_exit
    prof = [
        [0.0, P.z_swirl_bot],   # axis, chamber bottom
        [Rs,  P.z_swirl_bot],   # chamber bottom rim
        [Rs,  P.z_contr_bot],   # chamber top (contraction start)
        [Rw,  0.0],             # orifice entrance (contraction end)
        [Rw,  P.z_orif_top],    # orifice top
        [Re,  P.z_exp_top],     # expansion end
    ]
    if P.L_exit > 1e-6:
        prof.append([Re, P.z_exit_top])     # exit top (gerades L_exit-Stueck)
    prof.append([0.0, P.z_exit_top])        # close across the top to the axis
    return trimesh.creation.revolve(np.array(prof) * M, sections=SECTIONS)


def _chamber_cylinder():
    """Capped chamber cylinder (Rs over the chamber z-band), metres."""
    h = (P.z_contr_bot - P.z_swirl_bot) * M
    cyl = trimesh.creation.cylinder(radius=P.Rs * M, height=h, sections=SECTIONS)
    cyl.apply_translation([0.0, 0.0, 0.5 * (P.z_swirl_bot + P.z_contr_bot) * M])
    return cyl


def _inlet_geom():
    """(inner_pt, outer_pt, dir) per inlet, metres, in the chamber band."""
    r_in = 0.5 * P.Di * M
    offset = P.Rs * M - r_in                       # tangent to the chamber wall
    feed = 5.0 * M                                 # feed pipe length
    z_in = 0.5 * (P.z_swirl_bot + P.z_contr_bot) * M   # chamber centre
    axes, geoms = [], []
    for i in range(P.n_inlet):
        ang = 2.0 * math.pi * i / P.n_inlet
        ca, sa = math.cos(ang), math.sin(ang)
        inner = np.array([-sa * offset, ca * offset, z_in])
        outer = np.array([ca * feed - sa * offset, sa * feed + ca * offset, z_in])
        d = outer - inner
        axes.append((inner, outer, d / (np.linalg.norm(d) + 1e-30)))
        cyl = trimesh.creation.cylinder(radius=r_in, height=feed,
                                        sections=max(24, SECTIONS // 4))
        cyl.apply_transform(trimesh.transformations.rotation_matrix(math.pi / 2, [0, 1, 0]))
        cyl.apply_translation([0.5 * feed, offset, z_in])
        cyl.apply_transform(trimesh.transformations.rotation_matrix(ang, [0, 0, 1]))
        geoms.append(cyl)
    return axes, geoms, r_in


def _bore_r_at(z):
    """Bore-Wall-Radius an z (original frame, mm) entlang des Expansion-Profils.
    Zwischen z_orif_top (R=R_wall) und z_exp_top (R=R_exit) linear interpoliert."""
    if z >= P.z_exp_top:
        return P.R_exit
    if z <= P.z_orif_top:
        return P.R_wall
    frac = (z - P.z_orif_top) / (P.z_exp_top - P.z_orif_top)
    return P.R_wall + frac * (P.R_exit - P.R_wall)


def _top_seal_solid():
    """Annulare Scheibe als (outer cylinder) MINUS (inner frustum), wobei der
    inner frustum EXAKT entlang der Bore-Wall im z-Range des Disks verlaeuft.
    Wird in `build_surfaces` per Boolean Union in den Bore-Body integriert ->
    eine einzige durchgaengige nozzle.stl, kein Spalt zwischen Disk und Bore.

    z-Range (original frame): z_low = z_exit_top - thick ... z_high = z_exit_top.
    Inner radii: R_in_low = bore_r_at(z_low), R_in_high = bore_r_at(z_high)
    (= R_exit bei z_exit_top mit L_exit=0). Outer: konstanter R_out.

    Beide Sub-Solids werden via revolve(profile) mit AXIS-SCHLUSS gebaut, sonst
    liefert revolve mit annulaerem Profil Torus-Topologie (euler=0), gilt als
    nicht-watertight in trimesh -> Boolean union crasht.
    """
    scale_exit = P.R_exit / P.R_wall
    R_out = P.R_outer * scale_exit * _math.sqrt(2.0) + TOP_SEAL_MARGIN_MM
    # Disk-Position WIE URSPRUENGLICH (z=[z_exit_top - thick, z_exit_top]),
    # damit die Nozzle nicht effektiv laenger wird. Disk-Innenwand R=R_exit
    # (kein Offset) -- der Trick gegen Boolean-Verschmelzung ist dass die
    # Disk NICHT mehr in boolean.union(bore, ...) eingebaut wird, sondern als
    # separate STL `top_seal.stl` exportiert + von snappy als zweite
    # refinementSurface eingebunden wird (siehe build_surfaces).
    z_high = P.z_exit_top
    z_low  = z_high - TOP_SEAL_THICK_MM
    R_in_low  = _bore_r_at(z_low)
    R_in_high = _bore_r_at(z_high)
    # Outer body: rectangular profile from axis to R_out -> revolved cylinder.
    outer = trimesh.creation.revolve(np.array([
        [0.0,   z_low],
        [R_out, z_low],
        [R_out, z_high],
        [0.0,   z_high],
    ]) * M, sections=SECTIONS)
    # Inner cutter: triangle-like profile from axis to bore-wall-frustum.
    # Slightly larger in z to ensure clean cut at top + bottom (avoid coplanar faces).
    eps = 1e-6
    inner = trimesh.creation.revolve(np.array([
        [0.0,        z_low  - eps],
        [R_in_low,   z_low  - eps],
        [R_in_high,  z_high + eps],
        [0.0,        z_high + eps],
    ]) * M, sections=SECTIONS)
    return trimesh.boolean.difference([outer, inner], engine=_ENGINE)


def _classify_inlet_faces(mesh, axes, r_in):
    fc, fn = mesh.triangles_center, mesh.face_normals
    mask = np.zeros(len(mesh.faces), dtype=bool)
    for _inner, outer, d in axes:
        mask |= (np.linalg.norm(fc - outer, axis=1) < 1.3 * r_in) & ((fn @ d) > 0.7)
    return mask


def _write_multisolid_stl(path, mesh, inlet_mask):
    v, fn = mesh.vertices, mesh.face_normals

    def emit(fh, name, ids):
        fh.write(f"solid {name}\n")
        for fi in ids:
            n = fn[fi]
            fh.write(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n    outer loop\n")
            for vi in mesh.faces[fi]:
                p = v[vi]
                fh.write(f"      vertex {p[0]:.8e} {p[1]:.8e} {p[2]:.8e}\n")
            fh.write("    endloop\n  endfacet\n")
        fh.write(f"endsolid {name}\n")

    with open(path, "w", newline="\n") as fh:
        emit(fh, "wall", np.where(~inlet_mask)[0])
        emit(fh, "inlet", np.where(inlet_mask)[0])


# ----------------------------------------------------------------------
# build
# ----------------------------------------------------------------------
def build_surfaces(tri_dir):
    os.makedirs(tri_dir, exist_ok=True)
    bore = _bore_solid()
    axes, inlets, r_in = _inlet_geom()

    # nozzle.stl: bore + inlets via boolean.union (OHNE Disk).
    # Die Disk wird separat als top_seal.stl exportiert (siehe unten). Damit
    # bleibt die Bohrungs-Wand R=R_exit, z=[20,25] vollstaendig in nozzle.stl
    # erhalten -- keine Boolean-Verschmelzung mit der Disk-Innenwand
    # (Martin 2026-06-07 abends).
    surfaces = [bore, *inlets]
    union = trimesh.boolean.union(surfaces, engine=_ENGINE)
    if not union.is_watertight:
        raise RuntimeError(f"nozzle union not watertight (euler={union.euler_number})")

    # open the bore at the outlet (z=z_exit_top) and chamberBottom (z=z_swirl_bot).
    fc = union.triangles_center
    rr = np.linalg.norm(fc[:, :2], axis=1)
    z_top, z_bot = P.z_exit_top * M, P.z_swirl_bot * M
    cap = (((fc[:, 2] > z_top - 1e-6) & (rr < P.R_exit * M)) |
           ((fc[:, 2] < z_bot + 1e-6) & (rr < P.Rs * M)))
    open_union = trimesh.Trimesh(vertices=union.vertices,
                                 faces=union.faces[~cap], process=False)

    inlet_mask = _classify_inlet_faces(open_union, axes, r_in)   # before flip (per-face)
    open_union.apply_transform(_R_FLIP)
    _write_multisolid_stl(os.path.join(tri_dir, "nozzle.stl"), open_union, inlet_mask)

    feat = trimesh.boolean.union([_chamber_cylinder(), *inlets], engine=_ENGINE)
    if not feat.is_watertight:
        raise RuntimeError(f"feat union not watertight (euler={feat.euler_number})")
    feat.apply_transform(_R_FLIP)
    feat.export(os.path.join(tri_dir, "feat_inlet.stl"))

    # AMBIENT: top-seal Disk als SEPARATE STL (nicht in nozzle.stl), damit die
    # Disk-Innenwand R=R_exit als eigenstaendige Surface erhalten bleibt und
    # snappy sie gemeinsam mit der bore-Wand auf dieselben Cell-Faces snappen
    # kann (kein 0.1mm-Spalt mehr, keine Boolean-Verschmelzung mehr).
    if P.AMBIENT:
        seal = _top_seal_solid()
        if not seal.is_watertight:
            raise RuntimeError(f"top seal not watertight (euler={seal.euler_number})")
        seal.apply_transform(_R_FLIP)
        seal.export(os.path.join(tri_dir, "top_seal.stl"))

    print(f"nozzle.stl: {len(open_union.faces)} faces "
          f"(wall {int((~inlet_mask).sum())}, inlet {int(inlet_mask.sum())}); "
          f"dropped {int(cap.sum())} cap faces")
    print(f"feat_inlet.stl: {len(feat.faces)} faces; "
          f"chamber z=[{P.z_swirl_bot:.1f},{P.z_contr_bot:.1f}]mm Rs={P.Rs}mm, "
          f"{P.n_inlet} inlets")
    if P.AMBIENT:
        scale_exit = P.R_exit / P.R_wall
        R_out_mm = P.R_outer * scale_exit * _math.sqrt(2.0) + TOP_SEAL_MARGIN_MM
        print(f"top_seal.stl SEPARATE (NOT merged): annular ring R[{P.R_exit:.2f},"
              f"{R_out_mm:.2f}]mm z=[{P.z_exit_top - TOP_SEAL_THICK_MM:.2f},"
              f"{P.z_exit_top:.2f}]mm (orig.); inner wall co-located with bore "
              f"wall R={P.R_exit}mm -- snappy snaps both to same mesh face; "
              f"covers square-corner R_corner="
              f"{P.R_outer * scale_exit * _math.sqrt(2.0):.2f}mm")
    return axes, r_in


_HDR = ("/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "FoamFile {{ version 2.0; format ascii; class dictionary; object {obj}; }}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")


def write_sfe(sysdir):
    """surfaceFeatureExtractDict: nur feat_inlet (Tangentialeinlass-Penetrationen).
    Die top-seal Disk-Kanten sind nicht mehr separat — die Annulus ist Teil der
    nozzle.stl, ihre Kanten (z_high/Disk-Top und z_low/Disk-Bottom) sind dort
    natuerliche Knick-Kanten und werden ggf. via feat_inlet's includedAngle
    abgedeckt, wenn die Geometry-Wand das hergibt. Falls Disk-Edges schlecht
    snappen, kann man hier ein zusaetzliches `nozzle.stl`-Feature-Extract
    nachschieben."""
    with open(os.path.join(sysdir, "surfaceFeatureExtractDict"), "w", newline="\n") as f:
        f.write(_HDR.format(obj="surfaceFeatureExtractDict"))
        f.write('\nfeat_inlet.stl\n{\n'
                '    extractionMethod    extractFromSurface;\n'
                '    extractFromSurfaceCoeffs { includedAngle 150; }\n'
                '    subsetFeatures { nonManifoldEdges no; openEdges no; }\n'
                '    writeObj no;\n}\n')


def write_meshquality(sysdir):
    with open(os.path.join(sysdir, "meshQualityDict"), "w", newline="\n") as f:
        f.write(_HDR.format(obj="meshQualityDict"))
        f.write('\n#includeEtc "caseDicts/meshQualityDict"\n\n'
                'maxNonOrtho         65;\nmaxBoundarySkewness 20;\n'
                'maxInternalSkewness 4;\nminTetQuality       1e-30;\n'
                'minVol              1e-30;\nminDeterminant      0.001;\n'
                'nSmoothScale        4;\nerrorReduction      0.75;\n')


def write_toposet(sysdir):
    """cellSet 'chamberRefine' = the chamber outer-raster cellZones (zoneToCell
    regex). blockMesh makes one cellZone per named block (sketch3d.py); refineMesh
    maps zones to child cells, so this same dict re-selects after each refine pass."""
    pats = " ".join(f'"Out_.*_{seg}"' for seg in PREREFINE_SEGMENTS)
    text = (_HDR.format(obj="topoSetDict") + "\n"
            "actions\n(\n"
            "    {\n"
            "        name    chamberRefine;\n"
            "        type    cellSet;\n"
            "        action  new;\n"
            "        source  zoneToCell;\n"
            f"        zones   ( {pats} );\n"
            "    }\n"
            ");\n"
            "// " + "*" * 60 + " //\n")
    with open(os.path.join(sysdir, "topoSetDict"), "w", newline="\n") as f:
        f.write(text)


def write_refinemesh(sysdir):
    """refineMeshDict: isotropic hex split (2x2x2) of cellSet 'chamberRefine'."""
    text = (_HDR.format(obj="refineMeshDict") + "\n"
            "set             chamberRefine;\n"
            "coordinateSystem global;\n"
            "globalCoeffs    { tan1 (1 0 0); tan2 (0 1 0); }\n"
            "directions      ( tan1 tan2 normal );\n"
            "useHexTopology  true;\n"
            "geometricCut    false;\n"
            "writeMesh       false;\n"
            "// " + "*" * 60 + " //\n")
    with open(os.path.join(sysdir, "refineMeshDict"), "w", newline="\n") as f:
        f.write(text)


def write_snappy(sysdir):
    z_mid = 0.5 * (P.z_swirl_bot + P.z_contr_bot) * M
    Di_m = P.Di * M
    regions, zone_refs = [], []
    if PREREFINE_LEVEL > 0:
        # Pre-refine endform: the chamber outer raster is already refined to level
        # PREREFINE_LEVEL by topoSet+refineMesh (run before snappy). No region
        # boxes -> snappy only carves the inlets and snaps the bore wall. wall_level
        # 0: the bore already lies arc-exact on the O-grid. The inlet jet + its
        # penetration edge are held at ABSOLUTE level 2 (~18 cells across Di,
        # matching the validated box variant) regardless of the chamber level:
        # relative add = 2 - PREREFINE_LEVEL (>= 0).
        wall_level = 0
        inlet_level = max(0, 2 - PREREFINE_LEVEL)
        feat_level = max(0, 2 - PREREFINE_LEVEL)
    else:
        # legacy box mode (validated variant #6): 4 BOXES tile the chamber square
        # outer raster (outside +-c) over chamber+contraction -> uniform absolute
        # level 2; inlet surface + features level 2; round O-grid core stays 0.
        wall_level, inlet_level, feat_level = 0, 2, 2
        c_m = P.c * M
        out = P.Rs * M + 2.0 * Di_m             # cover the chamber wall + inlet reach
        z0 = -P.z_orif_bot * M - 0.25 * Di_m    # down to the orifice transition (+ margin)
        z1 = -P.z_swirl_bot * M + 0.3 * Di_m
        _boxes = {
            "E": ((c_m, -out, z0), (out, out, z1)),
            "W": ((-out, -out, z0), (-c_m, out, z1)),
            "N": ((-c_m, c_m, z0), (c_m, out, z1)),
            "S": ((-c_m, -out, z0), (c_m, -c_m, z1)),
        }
        for nm, (mn, mx) in _boxes.items():
            regions.append(f'    chamberBox{nm} {{ type searchableBox; '
                           f'min ({mn[0]:.6g} {mn[1]:.6g} {mn[2]:.6g}); '
                           f'max ({mx[0]:.6g} {mx[1]:.6g} {mx[2]:.6g}); }}')
            zone_refs.append(f'        chamberBox{nm} {{ mode inside; levels ((1e15 2)); }}')

    # AMBIENT (Task 3 / Martin 2026-06-07): top-seal Disk ist eine SEPARATE
    # STL (top_seal.stl), nicht via boolean.union in nozzle.stl gemerged.
    # Damit bleibt die Bohrungs-Wand R=R_exit in nozzle.stl als eigenstaendige
    # Surface erhalten (z=[20, 25]mm). top_seal.stl liefert die Disk-Annulus-
    # Surfaces (Outer/Inner/Top/Bottom). Snappy snapt die co-located bore-Wand
    # und Disk-Innenwand bei R=R_exit auf dieselben Mesh-Faces -> kein Spalt,
    # keine Sliver, Disk-Annulus bleibt vom Bore-Fluid isoliert.
    #
    # KEIN sprayCone refinementRegion: AMR (adaptive mesh refinement) macht
    # die Spray-Lamelle/Tropfen-Aufloesung spaeter zur Laufzeit automatisch
    # (Martin 2026-06-07 abends).
    ambient_geom_block = ""
    ambient_surf_block = ""
    ambient_region_block = ""
    if P.AMBIENT:
        ambient_geom_block = (
            "\n    top_seal.stl\n"
            "    {\n"
            "        type triSurfaceMesh;\n"
            "        name top_seal;\n"
            "    }"
        )
        ambient_surf_block = (
            f"\n        top_seal\n"
            f"        {{\n"
            f"            level ({wall_level} {wall_level});\n"
            f"            patchInfo {{ type wall; name top_seal; }}\n"
            f"        }}"
        )

    nl = chr(10)
    text = f"""{_HDR.format(obj="snappyHexMeshDict")}
castellatedMesh true;
snap            true;
addLayers       false;

geometry
{{
    nozzle.stl
    {{
        type triSurfaceMesh;
        name nozzle;
        regions {{ wall {{ name wall; }} inlet {{ name inlet; }} }}
    }}{ambient_geom_block}
{nl.join(regions)}
}}

castellatedMeshControls
{{
    maxLocalCells       4000000;
    maxGlobalCells      8000000;
    minRefinementCells  5;
    nCellsBetweenLevels 4;
    allowFreeStandingZoneFaces true;
    resolveFeatureAngle 40;

    features ( {{ file "feat_inlet.eMesh"; level {feat_level}; }} );

    refinementSurfaces
    {{
        nozzle
        {{
            level ({wall_level} {wall_level});
            regions
            {{
                wall  {{ level ({wall_level} {wall_level}); patchInfo {{ type wall;  name wall; }} }}
                inlet {{ level ({inlet_level} {inlet_level}); patchInfo {{ type patch; name inlet; }} }}
            }}
        }}{ambient_surf_block}
    }}

    refinementRegions
    {{
{nl.join(zone_refs)}{ambient_region_block}
    }}

    locationInMesh (0 0 {(-z_mid):.6g});
}}

snapControls
{{
    nSmoothPatch 3; tolerance 2.0; nSolveIter 50; nRelaxIter 5;
    nFeatureSnapIter 15; implicitFeatureSnap false; explicitFeatureSnap true;
    multiRegionFeatureSnap false;
}}

addLayersControls
{{
    relativeSizes true; expansionRatio 1.2; finalLayerThickness 0.4;
    minThickness 0.05; layers {{}}; nGrow 0; featureAngle 130; nRelaxIter 5;
    nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10;
    maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50;
}}

meshQualityControls {{ #include "meshQualityDict" }}

writeFlags ( scalarLevels );
mergeTolerance 1e-6;
// ************************************************************************* //
"""
    with open(os.path.join(sysdir, "snappyHexMeshDict"), "w", newline="\n") as f:
        f.write(text)


def write_allmesh(case):
    """Reproducible meshing pipeline following the logs/ convention. Runs
    PREREFINE_LEVEL passes of topoSet+refineMesh before snappy; A3_PREREFINE env
    overrides the pass count (0 -> legacy box mode, no pre-refine)."""
    allmesh = ('#!/bin/sh\n'
               'cd "${0%/*}" || exit 1\n'
               'rm -rf logs constant/polyMesh && mkdir -p logs\n'
               'blockMesh > logs/blockMesh.log 2>&1 || exit 1\n'
               f'PRE=${{A3_PREREFINE:-{PREREFINE_LEVEL}}}\n'
               'i=1\n'
               'while [ "$i" -le "$PRE" ]; do\n'
               '    topoSet -dict system/topoSetDict > logs/topoSet.$i.log 2>&1 || exit 1\n'
               '    refineMesh -dict system/refineMeshDict -overwrite > logs/refineMesh.$i.log 2>&1 || exit 1\n'
               '    i=$((i + 1))\n'
               'done\n'
               'surfaceFeatureExtract > logs/surfaceFeatureExtract.log 2>&1 || exit 1\n'
               'snappyHexMesh -overwrite > logs/snappyHexMesh.log 2>&1 || exit 1\n'
               'checkMesh -allGeometry -allTopology > logs/checkMesh.log 2>&1\n'
               'echo "meshing done (PRE=$PRE); see logs/checkMesh.log"\n')
    allclean = ('#!/bin/sh\n'
                'cd "${0%/*}" || exit 1\n'
                'rm -rf logs constant/polyMesh constant/extendedFeatureEdgeMesh\n'
                'rm -f constant/triSurface/*.eMesh\n'
                'rm -rf [0-9]* processor* VTK postProcessing\n'
                'echo cleaned\n')
    with open(os.path.join(case, "Allmesh"), "w", newline="\n") as f:
        f.write(allmesh)
    with open(os.path.join(case, "Allclean"), "w", newline="\n") as f:
        f.write(allclean)


if __name__ == "__main__":
    _flags = [a for a in sys.argv[1:] if a.startswith("--")]
    _pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    casename = _pos[0] if _pos else "case_a3"
    dicts_only = "--dicts-only" in _flags   # reuse existing STLs (geometry unchanged)
    case = os.path.join(os.path.dirname(os.path.abspath(__file__)), casename)
    if not dicts_only:
        axes, r_in = build_surfaces(os.path.join(case, "constant", "triSurface"))
    sysdir = os.path.join(case, "system")
    os.makedirs(sysdir, exist_ok=True)
    write_sfe(sysdir)
    write_meshquality(sysdir)
    if PREREFINE_LEVEL > 0:
        write_toposet(sysdir)
        write_refinemesh(sysdir)
    write_snappy(sysdir)
    write_allmesh(case)
    _extra = " +topoSet+refineMesh" if PREREFINE_LEVEL > 0 else ""
    _stls = "" if dicts_only else " + STLs"
    print(f"wrote snappy/sfe/meshQuality{_extra} + Allmesh{_stls} into {case} "
          f"(PREREFINE_LEVEL={PREREFINE_LEVEL})")
