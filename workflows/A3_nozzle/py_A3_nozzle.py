#!/usr/bin/env python3
"""A3_nozzle worker -- Zhang 2023 pressure-swirl nozzle (CFDPredict Testfall 3).

Per run: build the validated variant #6 hybrid mesh (structured O-grid bore +
snappyHexMesh-carved tangential inlets, prismatic chamber + body-fitted exit,
-z spray) via gen_a3_mesh.py, then run interFoam (VOF, RNG k-epsilon, water/air)
from the case_template, and extract the two CFDPredict objectives:

    film      liquid film thickness at the orifice exit  [mm]  (annular, from
              areaIntegrate(alpha.water) on the filmPlane sampledSurface)
    air_core  swirl-section diameter on the same plane = 2*R_core [mm]
    U_theta   circumferential (swirl) velocity mean/max from swirlU surfaces [m/s]
    vordruck  supply pressure = areaAverage(p) at the inlet patch  [MPa]
    cone      spray cone-angle = mean(2*atan(r_alpha=0.25 / |z-z_exit|)) over
              5 sampledSurface-planes 5..25 mm downstream of the exit  [deg]
              -- only with enable_ambient=1 (2 zusaetzliche z-Stationen).

Worker contract (WORKFLOW_AUTHORING.md section 4):
    python py_A3_nozzle.py <series_dir>
  - <series_dir>/interface.json holds the injected 6-tuple schema (values in the
    default slot, inject_schema: true).
  - progress/<step>_info.txt = live log; result PNG into progress/.
  - honour <series_dir>/command_kill (exit clean + touch command_finished).
  - touch <series_dir>/command_finished as the very last step (always).

Heavy VOF pipeline -> target host is cfdtools (128 cores). Requires OpenFOAM
v2512, plus trimesh + manifold3d for gen_a3_mesh.py. The build/extract logic is
offline-testable (do_run=False); the solver section runs only where OF exists.
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import foam_dictionary
import py_OF_utils
from py_OF_utils import (touch, run_OF_utility, run_solver_copy_progress,
                         clean_processor_directories)
try:
    from py_OF_utils import run_pvbatch_utility
except ImportError:
    run_pvbatch_utility = None  # pvbatch not available -- skipped at runtime

HERE = Path(__file__).resolve().parent
OF_BASHRC = os.environ.get("STREAMLIT_OPENFOAM_BASHRC",
                           "/opt/OpenFOAM/OpenFOAM-v2512/etc/bashrc")

# Zhang Fig.1b design parameters exposed in the form -> A3_<name> env overrides
# consumed by the vendored params.py. (float unless listed in _INT_PARAMS.)
GEOM_PARAMS = ("Di", "Ds", "Ls", "alpha", "Do", "Lo", "Lk", "n_inlet",
               "L_exit", "R_exit")
MESH_PARAMS = ("N_Do", "bl_cells")
AMBIENT_PARAMS = ("L_amb", "R_amb")     # only forwarded when enable_ambient=1
_INT_PARAMS = {"n_inlet", "N_Do", "bl_cells"}

# Form-Label -> internal param name. Zhang DOE-A/B/C/D get a one-letter
# prefix in the UI; other params drop the workflow prefix entirely. R_exit
# is NOT in the form anymore -- the user picks the expansion half-angle
# theta and the worker computes R_exit = Do/2 + Lk * tan(theta/2).
GEOM_LABELS = {
    "Di":      "Di [mm]",
    "Ds":      "D - Ds [mm]",
    "Ls":      "Ls [mm]",
    "alpha":   "C - alpha [deg]",
    "Do":      "B - Do [mm]",
    "Lo":      "Lo [mm]",
    "Lk":      "Lk [mm]",
    "n_inlet": "n_inlet",
    "L_exit":  "L_exit [mm]",
    # "R_exit" excluded -- derived from "A - theta [deg]" below.
}
THETA_LABEL = "A - theta [deg]"

# Cone-angle sampledSurface-planes: 5 stueck, 5/10/15/20/25 mm stromab des
# Nozzle-Exits. Die absolute z-Position haengt von der Geometrie ab (z_exit =
# -(Lo+Lk+L_exit)) -> der Worker patcht controlDict zur Bauzeit und uebergibt
# z_exit an _read_cone_alpha bei der Auswertung.
CONE_OFFSETS_MM = (5, 10, 15, 20, 25)
CONE_PLANE_NAMES = ("coneZ05", "coneZ10", "coneZ15", "coneZ20", "coneZ25")

# Zhang 2023 L16 orthogonal-array factor levels (paper sec. 4.3). A = expansion
# angle beta (OUTLET, deg) -> drives R_exit; B = straight-section diameter Do (mm);
# C = contraction angle alpha (chamber->straight, deg); D = swirl-section Ds (mm).
_DOE_A = {1: 0.0, 2: 10.0, 3: 20.0, 4: 30.0}   # expansion angle beta [deg]
_DOE_B = {1: 4.0, 2: 5.0, 3: 6.0, 4: 7.0}      # Do [mm]
_DOE_C = {1: 30.0, 2: 45.0, 3: 60.0, 4: 75.0}  # contraction angle alpha [deg]
_DOE_D = {1: 9.0, 2: 10.0, 3: 11.0, 4: 12.0}   # Ds [mm]


def parse_doe_code(code: str):
    """Decode a Zhang L16 label 'A_iB_jC_kD_l' (also 'A4B1C2D4') into geometry.
    Returns {beta, Do, alpha, Ds, label} or None if the string is not a DOE code.
    Raises ValueError on an out-of-range level."""
    m = re.search(r"A_?(\d)\s*B_?(\d)\s*C_?(\d)\s*D_?(\d)", code.strip(), re.I)
    if not m:
        return None
    a, b, c, d = (int(g) for g in m.groups())
    for lv, tab, nm in ((a, _DOE_A, "A"), (b, _DOE_B, "B"), (c, _DOE_C, "C"), (d, _DOE_D, "D")):
        if lv not in tab:
            raise ValueError("DOE level %s%d out of range 1..4" % (nm, lv))
    return {"beta": _DOE_A[a], "Do": _DOE_B[b], "alpha": _DOE_C[c], "Ds": _DOE_D[d],
            "label": "A%dB%dC%dD%d" % (a, b, c, d)}


# ---------------------------------------------------------------------------
# input handling (same convention as the zetaBrand worker)
# ---------------------------------------------------------------------------
def check_args() -> str:
    if len(sys.argv) != 2:
        sys.exit("usage: python %s <series_dir>" % os.path.basename(sys.argv[0]))
    return sys.argv[1]


def _flatten_with_all_defaults(schema: dict) -> dict:
    """Walk EVERY subdict-branch (not just the selected one) and return a flat
    dict of defaults. This lets the worker rely on every field being present
    even when the user kept a parent subdict on the non-default branch (e.g.
    "Show details" = "compact" hides the "advanced" tree -- but the worker
    still needs the Geometry/Mesh/Ambient/Flow defaults from there).
    """
    out = {}
    for key, entry in schema.items():
        if isinstance(entry, list) and len(entry) >= 5 and isinstance(entry[4], dict):
            out[key] = entry[1]
            for _branch_key, branch_entry in entry[4].items():
                sub_schema = branch_entry[1] if isinstance(branch_entry, list) and len(branch_entry) >= 2 else None
                if isinstance(sub_schema, dict):
                    out.update(_flatten_with_all_defaults(sub_schema))
        else:
            out[key] = entry[1] if isinstance(entry, list) and len(entry) >= 2 else entry
    return out


def unroll_interface(series_path: str) -> dict:
    with open(os.path.join(series_path, "interface.json"), encoding="utf-8") as f:
        raw = json.load(f)
    # Merge order: schema defaults (covers ALL branches) <- active selection
    # (only the user's currently-selected subdict, but with user-edited values).
    flat_defaults = _flatten_with_all_defaults(raw)
    flat_active = foam_dictionary.unroll_dict(raw)
    flat = {**flat_defaults, **flat_active}
    with open(os.path.join(series_path, "interfaceUnrolled2.json"), "w",
              encoding="utf-8") as f:
        f.write(json.dumps(flat))
    return flat


def update_progress(progress_path: str, info: list[str], step: str = "1") -> None:
    try:
        with open(os.path.join(progress_path, "%s_info.txt" % step), "w",
                  encoding="utf-8") as f:
            f.write("\n".join(info) + "\n")
    except Exception as e:  # noqa: BLE001
        print("update_progress failed (non-fatal): %s" % e)


def killed(series_path: str) -> bool:
    return os.path.isfile(os.path.join(series_path, "command_kill"))


def _write_overview(case_dir, progress_path: str, *,
                    project: str, doe_code: str, geom: dict, mesh: dict,
                    flow: dict, toggles: dict, ambient: dict | None) -> None:
    """Render the run's overview (table + Zhang-Fig.1b sketch) as PNG + HTML.

    PNG is also mirrored to progress/ so the streamlit live view shows it.
    Non-fatal on any failure -- the case build is what matters.
    """
    import shutil
    params = {
        "project": project, "doe_code": doe_code,
        "geom": geom, "mesh": mesh, "flow": flow,
        "toggles": toggles, "ambient": ambient,
    }
    try:
        case_path = Path(case_dir)
        (case_path / ".overview_params.json").write_text(
            json.dumps(params), encoding="utf-8")
        # import locally to keep the worker importable without matplotlib at test time
        import importlib.util as _ilu
        spec = _ilu.spec_from_file_location("_a3_make_overview",
                                            case_path / "make_overview.py")
        mod = _ilu.module_from_spec(spec)
        spec.loader.exec_module(mod)
        info_md = HERE / "info.md"
        mod.render_overview_png(params, case_path / "overview.png")
        mod.render_overview_html(params, info_md if info_md.exists() else None,
                                  case_path / "overview.html")
        # Mirror PNG to progress/ so the streamlit live view picks it up
        shutil.copy2(case_path / "overview.png",
                     os.path.join(progress_path, "0_overview.png"))
    except Exception as e:  # noqa: BLE001
        print("overview generation failed (non-fatal): %s" % e)


# ---------------------------------------------------------------------------
# case build
# ---------------------------------------------------------------------------
def build_case_dir(case_dir: Path, *, geom: dict, mesh: dict, ambient: dict | None,
                   u_inlet: float, end_time: float, n_proc: int) -> None:
    """Copy the interFoam VOF case_template, generate the variant-#6 mesh into it
    (fresh subprocess so params.py picks up the A3_* env), then patch the run
    controls (inlet velocity, endTime, film-plane z, decomposition).
    ambient=None -> classic #6; ambient={L_amb,R_amb} -> 2 z-stations downstream
    of the exit (frustum-aufgeweitet) fuer den Sprueh-Kegelwinkel."""
    if case_dir.exists():
        shutil.rmtree(case_dir)
    shutil.copytree(HERE / "case_template", case_dir)

    env = dict(os.environ)
    for k in GEOM_PARAMS:
        env["A3_" + k] = repr(int(geom[k])) if k in _INT_PARAMS else repr(float(geom[k]))
    for k in MESH_PARAMS:
        env["A3_" + k] = str(int(mesh[k]))
    if ambient:
        env["A3_AMBIENT"] = "1"
        for k in AMBIENT_PARAMS:
            env["A3_" + k] = repr(float(ambient[k]))
        if "spray_halfangle" in ambient:
            env["A3_SPRAY_HALFANGLE_DEG"] = repr(float(ambient["spray_halfangle"]))
        if "cone_half_angle_deg" in ambient:
            env["A3_cone_half_angle_deg"] = repr(float(ambient["cone_half_angle_deg"]))
        if "inner_amb_shrink" in ambient:
            env["A3_inner_amb_shrink"] = repr(float(ambient["inner_amb_shrink"]))
        topo = str(ambient.get("topology", "auto")).strip().lower()
        if topo in ("v1", "v2"):
            env["A3_TOPOLOGY"] = topo
        # topo == "auto" -> params.py setzt TOPOLOGY=v2 wegen A3_AMBIENT=1
    r = subprocess.run([sys.executable, str(HERE / "gen_a3_mesh.py"), str(case_dir)],
                       cwd=str(HERE), env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("gen_a3_mesh failed:\n%s\n%s" % (r.stdout, r.stderr))

    # inlet velocity (surfaceNormalFixedValue refValue, inward = negative)
    u_path = case_dir / "0" / "U"
    txt = u_path.read_text(encoding="utf-8")
    txt = re.sub(r"refValue\s+uniform\s+-?[0-9.eE+-]+;",
                 "refValue        uniform -%g;" % u_inlet, txt)
    u_path.write_text(txt, encoding="utf-8", newline="\n")

    # controlDict: endTime + film-plane z + cone-plane z. Film/swirlU sind auf der
    # Orifice-Ebene (~0.95*Lo). Cone-Planes liegen 5..25 mm stromabwaerts vom Exit
    # (= -(Lo+Lk+L_exit)). Alle Plane-Punkte werden NAME-anchored gepatcht (sonst
    # ueberschreibt ein globaler `point (0 0 z)`-Replace die anderen).
    z_film  = -0.95 * float(geom["Lo"]) * 1e-3
    z_exit  = -(float(geom["Lo"]) + float(geom["Lk"]) + float(geom["L_exit"])) * 1e-3
    cd = case_dir / "system" / "controlDict"
    txt = cd.read_text(encoding="utf-8")
    txt = re.sub(r"endTime\s+[0-9.eE+-]+;", "endTime         %g;" % end_time, txt)

    def _patch_plane_z(text, region_name, z_val):
        """Replace the FIRST `point ( 0 0 X )` after the named anchor block.
        count=1 + lazy `.*?` -> nur DAS erste point hinter dem Anker wird ersetzt,
        auch wenn der pointAndNormalDict-Block geschachtelt drin steckt (filmPlane
        -> sampledSurfaceDict -> pointAndNormalDict)."""
        return re.sub(
            r"(%s\s*\{.*?)point\s+\([^)]*\)" % region_name,
            r"\1point ( 0 0 %g )" % z_val, text, count=1, flags=re.S)

    # filmPlane FO (Anker: FO-Name) und swirlU.orificePlane (Anker: surface-Name).
    txt = _patch_plane_z(txt, "filmPlane",     z_film)
    txt = _patch_plane_z(txt, "orificePlane",  z_film)
    # Cone-Planes je nach Offset (5..25 mm stromab vom Exit).
    for off_mm, nm in zip(CONE_OFFSETS_MM, CONE_PLANE_NAMES):
        txt = _patch_plane_z(txt, nm, z_exit - off_mm * 1e-3)
    cd.write_text(txt, encoding="utf-8", newline="\n")

    # decomposeParDict ranks
    dpd = case_dir / "system" / "decomposeParDict"
    if dpd.is_file():
        txt = dpd.read_text(encoding="utf-8")
        txt = re.sub(r"numberOfSubdomains\s+\d+\s*;",
                     "numberOfSubdomains  %d;" % n_proc, txt)
        dpd.write_text(txt, encoding="utf-8", newline="\n")


# ---------------------------------------------------------------------------
# mesh + solve
# ---------------------------------------------------------------------------
def _apply_bl_refinement(case_dir: Path, geom: dict, passes: int) -> None:
    """Apply 2x refineWallLayer on the bore wall (wall_straight) + Lampenschirm
    (wall_expansion), with an unrefined upper buffer (wall_top) at the Drallkammer
    transition (z~0) to avoid the 48-negative-volume crash documented in Lesson 038-v1.

    Generates two dicts programmatically (topoSetDict_bl_split + createPatchDict_bl_split),
    runs the OF pipeline, then patches 0/* boundary fields to add wall_top/wall_straight/
    wall_expansion entries cloned from the existing outerWall BC.

    Geometry-aware z-cutoffs:
      - z_buffer = -0.0007 m (0.7 mm below chamber bottom z=0, 2-3 cell-layers safety)
      - z_lampenschirm = -(Lo - 0.1mm) (just above orifice end; only relevant if R_exit > Do/2)
    """
    if passes <= 0:
        return
    Lo_m = float(geom["Lo"]) * 1e-3
    z_buf = -0.0007
    z_lamp = -(Lo_m - 0.0001)

    (case_dir / "system" / "topoSetDict_bl_split").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }\n"
        "actions\n"
        "(\n"
        "    { name wall_top_faces; type faceSet; action new;\n"
        "      source patchToFace; sourceInfo { name outerWall; } }\n"
        "    { name wall_top_faces; type faceSet; action subset;\n"
        "      source boxToFace; sourceInfo { box (-0.020 -0.020 %g) (0.020 0.020 0.050); } }\n"
        "    { name wall_straight_faces; type faceSet; action new;\n"
        "      source patchToFace; sourceInfo { name outerWall; } }\n"
        "    { name wall_straight_faces; type faceSet; action subset;\n"
        "      source boxToFace; sourceInfo { box (-0.020 -0.020 -0.060) (0.020 0.020 %g); } }\n"
        "    { name wall_expansion_faces; type faceSet; action new;\n"
        "      source patchToFace; sourceInfo { name wall; } }\n"
        "    { name wall_expansion_faces; type faceSet; action subset;\n"
        "      source boxToFace; sourceInfo { box (-0.020 -0.020 -0.060) (0.020 0.020 %g); } }\n"
        ");\n" % (z_buf, z_buf, z_lamp),
        encoding="utf-8", newline="\n")

    (case_dir / "system" / "createPatchDict_bl_split").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object createPatchDict; }\n"
        "pointSync false;\n"
        "patches\n"
        "(\n"
        "    { name wall_top;       patchInfo { type wall; } constructFrom set; set wall_top_faces; }\n"
        "    { name wall_straight;  patchInfo { type wall; } constructFrom set; set wall_straight_faces; }\n"
        "    { name wall_expansion; patchInfo { type wall; } constructFrom set; set wall_expansion_faces; }\n"
        ");\n",
        encoding="utf-8", newline="\n")

    cdir = str(case_dir)
    run_OF_utility(cdir, "topoSet", params=["-dict", "system/topoSetDict_bl_split"])
    run_OF_utility(cdir, "createPatch", params=["-dict", "system/createPatchDict_bl_split", "-overwrite"])
    for _ in range(int(passes)):
        run_OF_utility(cdir, "refineWallLayer",
                       params=["-overwrite", "(wall_straight wall_expansion)", "0.5"])

    # Remove snappy leftovers: 0/cellLevel + 0/pointLevel carry the cell/point
    # count from BEFORE refineWallLayer, but the mesh has grown -- the next
    # decomposePar will FATAL with "Size N is not equal to the expected length M
    # file: 0/cellLevel/internalField" otherwise. Allrun_cfdtools.sh already
    # cleans these; the worker path must do the same.
    for stale in ("cellLevel", "pointLevel"):
        f = case_dir / "0" / stale
        if f.exists():
            f.unlink()

    # Patch 0/* boundary fields so that the mesh-patches and field-patches match.
    #
    # The BL-split creates up to three new mesh patches: wall_top, wall_straight,
    # wall_expansion. Whether they originate from an `outerWall` source (ambient-
    # variant) or from `wall` (classical innenströmung) doesn't matter -- they all
    # need a BC in 0/* fields.
    #
    # case_template/0/* historically also carries a `lexit_side` patch that is no
    # longer in the modern geometry. If we leave it in, decomposePar finds a field
    # patch that has no mesh counterpart, and the parallel solver later FATALs with
    # "cannot find file processor*/0/p_rgh" because field decomposition silently
    # bails out. So: strip lexit_side first, then ensure the three split-patches
    # exist (cloned from the wall BC).
    def _strip_block(text: str, patch_name: str) -> str:
        """Remove '<indent>patch_name { ... }' (single or multi-line) plus any
        trailing inline comment on the same line. Idempotent."""
        # Single-line form first: '    name { type X; value Y; }   // comment'
        text = re.sub(
            r"^[ \t]*" + re.escape(patch_name) + r"\s*\{[^}]*\}[ \t]*(?://[^\n]*)?\n?",
            "", text, flags=re.M)
        # Multi-line form: '    name\n    {\n  ...  }   // comment\n'
        text = re.sub(
            r"^[ \t]*" + re.escape(patch_name) + r"\s*\n[ \t]*\{[^}]*\}[ \t]*(?://[^\n]*)?\n?",
            "", text, flags=re.M)
        return text

    def _extract_wall_bc(text: str):
        """Return (indent, bc_body_with_braces, single_line_bool) for the `wall`
        patch BC. None if not found."""
        # Single-line: '    wall  { type X; ... }'
        m = re.search(r"^([ \t]*)wall(\s*\{[^}]*\})\s*$", text, flags=re.M)
        if m:
            return (m.group(1), m.group(2), True)
        # Multi-line: '    wall\n    {\n  ...  }'
        m = re.search(r"^([ \t]*)wall\s*\n[ \t]*(\{[^}]*\})", text, flags=re.M | re.S)
        if m:
            return (m.group(1), m.group(2), False)
        return None

    for fn in ("U", "k", "epsilon", "nut", "p_rgh", "alpha.water"):
        fp = case_dir / "0" / fn
        if not fp.exists():
            continue
        s = fp.read_text(encoding="utf-8")

        # 1) Drop legacy `lexit_side` (not in the modern mesh -> decomposePar mismatch).
        s = _strip_block(s, "lexit_side")

        # 2) Already patched on a prior call? Skip (idempotent).
        if "wall_expansion" in s:
            fp.write_text(s, encoding="utf-8", newline="\n")
            continue

        # 3) Ambient-variant: outerWall exists -> replace with the three new patches.
        m = re.search(r"^([ \t]*)outerWall(\s*\{[^}]*\})\s*$", s, flags=re.M)
        if m:
            indent, bc = m.group(1), m.group(2)
            block = (indent + "wall_top" + bc + "\n" +
                     indent + "wall_straight" + bc + "\n" +
                     indent + "wall_expansion" + bc)
            s = re.sub(r"^[ \t]*outerWall\s*\{[^}]*\}\s*$", block, s, count=1, flags=re.M)
            fp.write_text(s, encoding="utf-8", newline="\n")
            continue
        m = re.search(r"^([ \t]*)outerWall\s*\n[ \t]*\{([^}]*)\}", s, flags=re.M | re.S)
        if m:
            indent, body = m.group(1), m.group(2)
            block = (indent + "wall_top\n" + indent + "{" + body + "}\n" +
                     indent + "wall_straight\n" + indent + "{" + body + "}\n" +
                     indent + "wall_expansion\n" + indent + "{" + body + "}")
            s = re.sub(r"^[ \t]*outerWall\s*\n[ \t]*\{[^}]*\}", block, s, count=1, flags=re.M | re.S)
            fp.write_text(s, encoding="utf-8", newline="\n")
            continue

        # 4) Classical innenströmung: outerWall doesn't exist. Clone the `wall` BC
        #    and insert wall_top + wall_straight + wall_expansion after it.
        ext = _extract_wall_bc(s)
        if ext is None:
            # nothing to do (no wall patch found -- unusual)
            fp.write_text(s, encoding="utf-8", newline="\n")
            continue
        indent, bc, single_line = ext
        if single_line:
            block = ("\n" + indent + "wall_top" + bc +
                     "\n" + indent + "wall_straight" + bc +
                     "\n" + indent + "wall_expansion" + bc)
            s = re.sub(r"^([ \t]*wall\s*\{[^}]*\})\s*$", r"\1" + block,
                       s, count=1, flags=re.M)
        else:
            body = bc.strip()[1:-1]  # strip the surrounding { }
            block = ("\n" + indent + "wall_top\n" + indent + "{" + body + "}" +
                     "\n" + indent + "wall_straight\n" + indent + "{" + body + "}" +
                     "\n" + indent + "wall_expansion\n" + indent + "{" + body + "}")
            s = re.sub(r"^([ \t]*wall\s*\n[ \t]*\{[^}]*\})", r"\1" + block,
                       s, count=1, flags=re.M | re.S)
        fp.write_text(s, encoding="utf-8", newline="\n")

    # Re-checkMesh after refinement so the log captures the post-refine quality.
    run_OF_utility(cdir, "checkMesh", params=["-allGeometry", "-allTopology"])


def _apply_zhang_outlet_bcs(case_dir: Path) -> None:
    """Switch outlet BCs to Zhang-paper-conform: p_rgh outlet = fixedValue 0,
    alpha.water outlet = inletOutlet 0. Replaces the DHCAE non-reflecting variant
    (advective for p_rgh, zeroGradient for alpha.water -- Martins WIP).

    Lesson 029: k_avg drops -90%, dt grows +44% (fixedValue is Dirichlet without
    own time derivative) -- numerically cheaper + Zhang-aligned.
    """
    p_path = case_dir / "0" / "p_rgh"
    if p_path.exists():
        s = p_path.read_text(encoding="utf-8")
        s = re.sub(
            r"^(\s*)outlet\s*\{[^}]*\}",
            r"\1outlet         { type fixedValue; value uniform 0; }",
            s, count=1, flags=re.M | re.S)
        p_path.write_text(s, encoding="utf-8", newline="\n")
    a_path = case_dir / "0" / "alpha.water"
    if a_path.exists():
        s = a_path.read_text(encoding="utf-8")
        s = re.sub(
            r"^(\s*)outlet\s*\{[^}]*\}",
            r"\1outlet\n\1{\n\1    type            inletOutlet;\n"
            r"\1    inletValue      uniform 0;\n\1    value           uniform 0;\n\1}",
            s, count=1, flags=re.M | re.S)
        a_path.write_text(s, encoding="utf-8", newline="\n")


def _strip_prefill(case_dir: Path) -> None:
    """Disable setFields prefill: replace regions(...) with regions ();. alpha=0
    everywhere (Lesson 028). Aircore builds up naturally from inlet ramp.

    Greedy `.*` (not `.*?`) is required: the region entry nests its own
    parenthesised `fieldValues (...)` list, so a non-greedy match stops at
    that INNER closing `);` and truncates mid-dict, leaving the rest of the
    original region (closing `}` + outer `);`) as orphaned text -- a FOAM
    FATAL IO ERROR at the next `setFields` call ("Unexpected '}'"). Greedy
    DOTALL backtracks from the end of the string, landing on the true final
    `);` of the top-level `regions(...)` list instead. Caught 2026-09
    (Dresden-export smoke test): the default `Disable prefill (0/1) = 1`
    means this ran, silently, on nearly every A3_nozzle case to date --
    `run_OF_utility()` doesn't check setFields' exit code, so the crash was
    swallowed and alpha.water just kept its 0/-template default (uniform 0),
    which happens to equal the intended "no prefill" outcome. Correctness of
    prior results is very likely unaffected; the log noise and confusing
    "success" despite a FATAL error are what this fixes."""
    sfd = case_dir / "system" / "setFieldsDict"
    if not sfd.is_file():
        return
    s = sfd.read_text(encoding="utf-8")
    s = re.sub(r"regions\s*\(.*\)\s*;", "regions ();", s, count=1, flags=re.S)
    sfd.write_text(s, encoding="utf-8", newline="\n")


def _run_pv_script(case_dir: Path, script_name: str) -> None:
    """Invoke a pvbatch script in the case_dir. Non-fatal if pvbatch unavailable."""
    if run_pvbatch_utility is None:
        return
    script = case_dir / script_name
    if not script.is_file():
        print("WARNING: pvbatch script %s missing, skipping" % script_name)
        return
    try:
        run_pvbatch_utility(str(case_dir),
                            params=[str(script), "--force-offscreen-rendering"])
    except Exception as e:  # noqa: BLE001
        print("WARNING: pvbatch %s failed (non-fatal): %s" % (script_name, e))


def _copy_pv_pngs(series_path: str, case_dir: Path, filenames: list) -> None:
    """Copy PNGs from case_dir to progress/ with case-name prefix."""
    progress = os.path.join(series_path, "progress")
    for fn in filenames:
        src = case_dir / fn
        if not src.is_file():
            continue
        dst = os.path.join(progress, "%s_%s" % (case_dir.name, fn))
        try:
            shutil.copy2(src, dst)
        except OSError as e:
            print("WARNING: copy %s -> %s failed: %s" % (src, dst, e))


# pvbatch script lists, split by data-dependency:
#   MESH-only: need only blockMesh + snappyHexMesh + refineWallLayer output;
#              rendered BEFORE decomposePar so the user sees the final mesh
#              even if the solver fails to start (Lesson 2026-06-08).
#   FIGURE: need solver fields (U, p_rgh, alpha.water); rendered AFTER the
#           solver finishes (reconstructPar -latestTime).
_PV_MESH_SCRIPTS = [
    ("pv_mesh_side.py",   "pv_mesh_side.png"),
    ("pv_mesh_zoom_x.py", "pv_mesh_zoom_x.png"),
    ("pv_mesh_zoom_z.py", "pv_mesh_zoom_z.png"),
]
_PV_FIGURE_SCRIPTS = [
    ("pv_figure4.py", "pv_figure4.png"),
    ("pv_figure7.py", "pv_figure7.png"),
]


def _ensure_foam_anchor(case_dir: Path) -> None:
    """ParaView's OpenFOAMReader needs a <casename>.foam anchor file."""
    foam_anchor = case_dir / (case_dir.name + ".foam")
    if not foam_anchor.exists():
        foam_anchor.touch()


def _render_pv_group(case_dir: Path, scripts: list) -> None:
    """Run a list of pvbatch scripts and copy their PNGs to progress/."""
    _ensure_foam_anchor(case_dir)
    for script, _png in scripts:
        _run_pv_script(case_dir, script)
    _copy_pv_pngs(str(case_dir.parent), case_dir, [png for _, png in scripts])


def _run_pv_mesh_renders(case_dir: Path) -> None:
    """Three mesh-view PNGs (side / zoom_x / zoom_z) -- after the mesh is final."""
    _render_pv_group(case_dir, _PV_MESH_SCRIPTS)


def _run_pv_figure_renders(case_dir: Path) -> None:
    """Two solver-result PNGs (Zhang Fig.4 / Fig.7) -- after reconstructPar."""
    _render_pv_group(case_dir, _PV_FIGURE_SCRIPTS)


def mesh_and_solve(series_path: str, case_dir: Path, n_proc: int,
                   geom: dict | None = None, bl_passes: int = 2,
                   zhang_bcs: bool = True, disable_prefill: bool = True) -> None:
    cdir = str(case_dir)
    run_OF_utility(cdir, "blockMesh")
    run_OF_utility(cdir, "surfaceFeatureExtract")
    run_OF_utility(cdir, "snappyHexMesh", params=["-overwrite"])
    # createPatch nur wenn der dict existiert -- bei der jetzigen Geometrie (Disk
    # ins nozzle.stl gemerged) wird er nicht mehr erzeugt, bleibt aber kompatibel
    # zu aelteren Cases mit getrennter top.stl.
    if (case_dir / "system" / "createPatchDict").is_file():
        run_OF_utility(cdir, "createPatch", params=["-overwrite"])
    run_OF_utility(cdir, "checkMesh", params=["-allGeometry", "-allTopology"])
    # BL refinement (Lesson 038-v2): outerWall -> wall_top (Puffer) + wall_straight (refine)
    # + wall -> wall + wall_expansion (refine). 2x refineWallLayer = 4x naeher zur Wand.
    if bl_passes > 0 and geom is not None:
        _apply_bl_refinement(case_dir, geom=geom, passes=bl_passes)
    # Zhang-konforme outlet-BCs (Lesson 029) + no-prefill (Lesson 028).
    if zhang_bcs:
        _apply_zhang_outlet_bcs(case_dir)
    if disable_prefill:
        _strip_prefill(case_dir)
    # setFields: pre-fill bore body with water (alpha=1) to remove the
    # cold-start velocity spike that drives kEpsilon eps-production into
    # detonation. Serial -> runs before decomposePar so the seeded
    # alpha.water gets distributed to all processor*/0/. Conditional on
    # the dict so legacy cases without setFieldsDict still build.
    if (case_dir / "system" / "setFieldsDict").is_file():
        run_OF_utility(cdir, "setFields")
    # Mesh-Vorschau-PNGs JETZT (Mesh ist final, vor decomposePar). Damit hat
    # der User die Bilder selbst dann, wenn mpirun-Start spaeter scheitert
    # (z.B. konkurrierende 128-Rank-Jobs auf einem produktiven HPC-Server).
    try:
        _run_pv_mesh_renders(case_dir)
    except Exception as e:  # noqa: BLE001
        print("pv mesh renders failed (non-fatal): %s" % e)
    run_OF_utility(cdir, "decomposePar", params=["-force"])
    run_solver_copy_progress(series_path, ["residuals.png"], cdir, "interFoam",
                             log_file="logs/solver.log", n_proc=str(n_proc),
                             start_observer=True)
    # Solver-Erfolgs-Check: bei interFoam-Init-Crash (SIGFPE/SIGSEGV nach
    # "Reading field p_rgh") meldet mpirun mitunter rc=0 und solver.log
    # endet abrupt ohne FATAL-Block -- die Pipeline wuerde sonst still
    # weiterlaufen, summary.json zeigt alle Metriken = null, command_finished
    # wird gesetzt, und der Nutzer denkt "Erfolg". Hier proben wir, ob
    # tatsaechlich Zeitschritt-Output entstanden ist.
    if not _solver_wrote_timesteps(case_dir):
        raise RuntimeError(
            "Solver produced no time-step output. interFoam wahrscheinlich "
            "beim Init gecrasht (SIGFPE/SIGSEGV) -- pruefe logs/solver.log "
            "(typisch nach 'Reading field p_rgh' ohne FATAL-Block "
            "abgebrochen). Mesh + decomposePar liefen fehlerfrei, Crash "
            "liegt im Solver selbst (haeufige Ursachen: BC-Mismatch in "
            "0/, init-Felder mit NaN, ambient-Topologie-Instabilitaet, "
            "konkurrierender 128-Rank-Job auf demselben Server).")
    run_OF_utility(cdir, "reconstructPar", params=["-latestTime"])
    clean_processor_directories(cdir)


def _solver_wrote_timesteps(case_dir: Path) -> bool:
    """True if the solver wrote at least one t>0 time directory.

    Looks for time dirs both in the top-level case (post-reconstruct or
    serial run) and inside processor0/ (still decomposed). Anything with
    a numeric name strictly > 0 counts as "solver advanced past init".
    """
    candidates = [case_dir]
    proc0 = case_dir / "processor0"
    if proc0.is_dir():
        candidates.append(proc0)
    for parent in candidates:
        for child in parent.iterdir():
            if not child.is_dir():
                continue
            try:
                if float(child.name) > 0.0:
                    return True
            except ValueError:
                continue
    return False


# ---------------------------------------------------------------------------
# objective extraction
# ---------------------------------------------------------------------------
def _fo_series(case_dir: Path, fo_name: str) -> list[tuple[float, float]]:
    """All (time, value) rows from postProcessing/<fo>/<startTime>/*.dat
    (newest startTime dir; concatenate columns 1,2). Empty list if missing."""
    base = case_dir / "postProcessing" / fo_name
    if not base.is_dir():
        return []
    starts = sorted([d for d in base.iterdir() if d.is_dir()],
                    key=lambda d: float(d.name))
    rows: list[tuple[float, float]] = []
    for d in starts:
        for dat in sorted(d.glob("*.dat")):
            for ln in dat.read_text(encoding="utf-8", errors="ignore").splitlines():
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                parts = ln.split()
                try:
                    rows.append((float(parts[0]), float(parts[1])))
                except (ValueError, IndexError):
                    pass
    return rows


def _read_cone_alpha(case_dir: Path, z_exit_m: float) -> dict:
    """Sprueh-Kegelwinkel aus den 5 coneAlpha-Planes (latest writeTime).
    Pro Plane: faces mit alpha.water > 0.25 raussuchen, max-Radius r = sqrt(x^2+y^2)
    bestimmen, Cone-Angle = 2*atan(r/|z - z_exit|) [deg]. Mean+per-Plane in dict.
    z_exit_m = -(Lo+Lk+L_exit)*1e-3 (geflippter Frame), wird vom Aufrufer aus der
    aktuellen Geometrie berechnet; die Plane-Position ergibt sich aus dem
    Offset-Tisch CONE_OFFSETS_MM.
    Robust gegen: alpha-Spalte fehlt, kein nasses Face (-> Plane wird ausgelassen),
    nicht aktive Ambient (-> base dir fehlt; return {}).
    """
    base = case_dir / "postProcessing" / "coneAlpha"
    if not base.is_dir():
        return {}
    times = sorted([d for d in base.iterdir() if d.is_dir()],
                   key=lambda d: float(d.name))
    if not times:
        return {}
    latest = times[-1]
    per_plane: dict[str, dict] = {}
    angles_deg: list[float] = []
    try:
        import numpy as _np
        for plane_name, off_mm in zip(CONE_PLANE_NAMES, CONE_OFFSETS_MM):
            z_m = z_exit_m - off_mm * 1e-3
            # raw-Output: coneAlpha.coneZ05_alpha.water.raw oder coneZ05.raw,
            # je nach OF-Variante. Wir nehmen das erste matchende file.
            raws = (list(latest.glob("%s_alpha.water.raw" % plane_name)) +
                    list(latest.glob("%s.raw" % plane_name)))
            if not raws:
                continue
            d = _np.loadtxt(raws[0], comments="#")
            if d.ndim != 2 or d.shape[1] < 4:
                continue
            x, y, alpha = d[:, 0], d[:, 1], d[:, 3]
            wet = alpha > 0.25
            if not wet.any():
                per_plane[plane_name] = {"z_m": z_m, "n_wet": 0,
                                         "r_max_mm": None, "cone_deg": None}
                continue
            r = _np.hypot(x[wet], y[wet])
            r_max = float(r.max())
            arm = abs(z_m - z_exit_m)
            cone_deg = math.degrees(2.0 * math.atan(r_max / arm)) if arm > 0 else None
            if cone_deg is not None:
                angles_deg.append(cone_deg)
            per_plane[plane_name] = {"z_m": z_m, "n_wet": int(wet.sum()),
                                     "r_max_mm": r_max * 1e3, "cone_deg": cone_deg}
    except Exception as e:  # noqa: BLE001
        print("coneAlpha read failed (non-fatal): %s" % e)
        return {}
    mean_deg = sum(angles_deg) / len(angles_deg) if angles_deg else None
    return {"mean_deg": mean_deg, "planes": per_plane}


def _read_swirl_u(case_dir: Path) -> tuple[float | None, float | None]:
    """(mean, max) circumferential velocity |U_theta| [m/s] on the orifice plane,
    from the swirlU `surfaces` raw output (latest write). U_theta about the z spray
    axis = (-Ux*y + Uy*x)/r. Returns (None, None) if the raw surface is absent."""
    base = case_dir / "postProcessing" / "swirlU"
    if not base.is_dir():
        return None, None
    times = sorted([d for d in base.iterdir() if d.is_dir()], key=lambda d: float(d.name))
    raws = sorted(times[-1].glob("*.raw")) if times else []
    if not raws:
        return None, None
    try:
        import numpy as _np
        d = _np.loadtxt(raws[0], comments="#")
        if d.ndim != 2 or d.shape[1] < 6:
            return None, None
        x, y, ux, uy = d[:, 0], d[:, 1], d[:, 3], d[:, 4]
        r = _np.hypot(x, y)
        m = r > 1e-9
        if not m.any():
            return None, None
        ut = _np.abs((-ux[m] * y[m] + uy[m] * x[m]) / r[m])
        return float(ut.mean()), float(ut.max())
    except Exception as e:  # noqa: BLE001
        print("swirlU read failed (non-fatal): %s" % e)
        return None, None


def extract_objectives(case_dir: Path, Do_mm: float, *,
                       z_exit_m: float | None = None) -> dict:
    """Zhang objectives from the finished VOF case:
      film_mm        annular film at the orifice from filmPlane areaIntegrate(alpha.water)
      air_core_mm    swirl-section (air-core) diameter on the same plane = 2*R_core
      u_theta_*      circumferential (swirl) velocity mean/max from the swirlU surface
      vordruck_MPa   supply pressure from inletPressure areaAverage(p)
    Returns last values + the inletPressure/filmPlane time series for the chart."""
    p_series = _fo_series(case_dir, "inletPressure")
    f_series = _fo_series(case_dir, "filmPlane")
    vordruck_mpa = p_series[-1][1] / 1e6 if p_series else None

    film_mm = air_core_mm = None
    if f_series:
        a_wet = f_series[-1][1]                      # liquid area [m^2]
        r_o = 0.5 * Do_mm * 1e-3                      # orifice radius [m]
        disc = r_o * r_o - max(a_wet, 0.0) / math.pi
        if disc >= 0.0:
            r_core = math.sqrt(disc)
            film_mm = (r_o - r_core) * 1e3            # annular film [mm]
            air_core_mm = 2.0 * r_core * 1e3          # air-core / swirl-section diameter [mm]
    u_mean, u_max = _read_swirl_u(case_dir)
    cone = _read_cone_alpha(case_dir, z_exit_m) if z_exit_m is not None else {}
    return {"film_mm": film_mm, "air_core_mm": air_core_mm,
            "u_theta_mean": u_mean, "u_theta_max": u_max,
            "vordruck_MPa": vordruck_mpa,
            "cone_angle_deg": cone.get("mean_deg") if cone else None,
            "cone_planes": cone.get("planes") if cone else None,
            "p_series": p_series, "f_series": f_series}


def parse_cells(case_dir: Path) -> int | None:
    log = case_dir / "logs" / "checkMesh.log"
    if not log.is_file():
        return None
    # First match is the "Mesh stats" block ("    cells: 835362"). Later
    # matches come from the "Failed checks" block ("Concave cells ...
    # number of cells: 1684") and are sub-counts, not the total.
    m = re.findall(r"cells:\s+(\d+)", log.read_text(encoding="utf-8", errors="ignore"))
    return int(m[0]) if m else None


def render_chart(png_path: str, obj: dict, project: str) -> str | None:
    p, f = obj["p_series"], obj["f_series"]
    if not p and not f:
        return None
    fig, ax1 = plt.subplots(figsize=(8, 5))
    if p:
        ax1.plot([t * 1e3 for t, _ in p], [v / 1e6 for _, v in p],
                 color="#1f4e79", lw=2, label="Vordruck (inlet p)")
    ax1.set_xlabel("time [ms]")
    ax1.set_ylabel("supply pressure [MPa]", color="#1f4e79")
    ax1.grid(alpha=0.3)
    if f:
        r_o = None
        ax2 = ax1.twinx()
        # convert wetted area -> annular film thickness for the second axis
        ys = []
        for _, a in f:
            # Do is encoded via the last extract; approximate with first if needed
            ys.append(a)
        ax2.plot([t * 1e3 for t, _ in f], ys, color="#b5651d", lw=1.5, ls="--",
                 label="wetted area [m^2]")
        ax2.set_ylabel("wetted area at orifice [m^2]", color="#b5651d")
    ax1.set_title("%s -- A3 pressure-swirl nozzle (interFoam VOF)" % project)
    fig.tight_layout()
    try:
        fig.savefig(png_path, dpi=140)
    finally:
        plt.close(fig)
    return png_path


# ---------------------------------------------------------------------------
# results zip
# ---------------------------------------------------------------------------
# Goal: download is a self-contained, ParaView-openable OpenFOAM case.
# Includes 0/, constant/, system/, logs/, postProcessing/ and all written
# time directories (e.g. the latest reconstructed snapshot). Skips
# processor*/ (parallel decomposition, redundant after reconstructPar and
# usually deleted by clean_processor_directories anyway) and a few state
# files used by the streamlit harness.
_ZIP_SKIP_TOP_PREFIXES = ("processor",)
_ZIP_SKIP_TOP_NAMES = {"command_kill", "command_finished", ".pv_script_failed"}


def zip_results(series_path: str, results_path: str, case_dir: Path) -> None:
    import zipfile
    series = Path(series_path)
    try:
        with zipfile.ZipFile(Path(results_path) / "results.zip", "w",
                             zipfile.ZIP_DEFLATED) as zf:
            # 1) Charts + per-step progress notes from the streamlit run.
            prog = series / "progress"
            if prog.is_dir():
                for p in prog.rglob("*"):
                    if p.is_file():
                        zf.write(p, p.relative_to(series))
            # 2) The case itself -- ALL top-level entries except processor*/
            #    and a few harness state files. That covers 0/, constant/,
            #    system/, logs/, postProcessing/, every written time
            #    directory, the .foam anchor, overview.png/html, pv_*.py.
            for child in sorted(case_dir.iterdir()):
                if child.name.startswith(_ZIP_SKIP_TOP_PREFIXES):
                    continue
                if child.name in _ZIP_SKIP_TOP_NAMES:
                    continue
                if child.is_dir():
                    for p in child.rglob("*"):
                        if p.is_file():
                            zf.write(p, p.relative_to(series))
                elif child.is_file():
                    zf.write(child, child.relative_to(series))
    except Exception as e:  # noqa: BLE001
        print("zip_results failed (non-fatal): %s" % e)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def run_worker(series_path: str, *, do_source_of: bool = True,
               do_run: bool = True) -> None:
    series_path = os.path.abspath(series_path).rstrip("/\\")
    progress_path = os.path.join(series_path, "progress")
    results_path = os.path.join(series_path, "results")
    for p in (progress_path, results_path):
        os.makedirs(p, exist_ok=True)

    info = ["py_A3_nozzle.py is running now."]
    update_progress(progress_path, info)

    j = unroll_interface(series_path)
    project = str(j.get("project name", os.path.basename(series_path)))
    geom = {k: float(j[GEOM_LABELS[k]]) for k in GEOM_LABELS}
    geom["n_inlet"] = int(round(geom["n_inlet"]))
    # The form takes theta (Zhang's expansion half-angle); R_exit is derived
    # so the downstream params.py contract stays unchanged.
    theta_deg = float(j[THETA_LABEL])
    geom["R_exit"] = 0.5 * geom["Do"] + geom["Lk"] * math.tan(math.radians(theta_deg) / 2.0)
    # optional DOE shortcut A_iB_jC_kD_l overrides Do(B)/alpha(C)/Ds(D) + the outlet
    # expansion beta(A) -> R_exit. Empty -> the individual A/B/C/D + theta fields stay.
    doe_note = ""
    doe_code = str(j.get("DOE code (optional, e.g. A4B1C2D4)", "")).strip()
    if doe_code:
        doe = parse_doe_code(doe_code)
        if doe is None:
            raise ValueError("DOE code %r not understood (expected A_iB_jC_kD_l)" % doe_code)
        geom["Do"], geom["alpha"], geom["Ds"] = doe["Do"], doe["alpha"], doe["Ds"]
        geom["R_exit"] = 0.5 * doe["Do"] + geom["Lk"] * math.tan(math.radians(doe["beta"]) / 2.0)
        doe_note = ("%s -> beta=%g Do=%g alpha=%g Ds=%g => R_exit=%.4g mm (beta=0 = straight)"
                    % (doe["label"], doe["beta"], doe["Do"], doe["alpha"], doe["Ds"], geom["R_exit"]))
    mesh = {"N_Do": int(float(j["base cells over Do (N_Do)"])),
            "bl_cells": int(float(j["BL/film band cells (bl_cells)"]))}
    u_inlet = float(j["inlet velocity [m/s]"])
    end_time = float(j["end time [s]"])
    n_proc = max(1, min(128, int(float(j.get("number of processors", 128)))))
    contact = str(j.get("contact email", "")).strip()

    # Optional ambient-Erweiterung (Task 3) -> 2 zusaetzliche z-Stationen
    # stromabwaerts des Nozzle-Exits + Sprueh-Kegelwinkel-Auswertung.
    # spray_halfangle steuert die Frustum-Aufweitung (R_amb auto wenn 0).
    enable_ambient = bool(int(float(j.get("enable ambient (0/1)", 0) or 0)))
    spray_halfangle = float(j.get("spray_halfangle [deg]", 47.0))
    topology = str(j.get("TOPOLOGY", "auto")).strip().lower()
    cone_half_angle = float(j.get("cone_half_angle [deg]", 40.0))
    inner_amb_shrink = float(j.get("inner_amb_shrink", 0.9))
    ambient = None
    if enable_ambient:
        L_amb_mm = float(j.get("L_amb [mm]", 30.0))
        R_amb_user = float(j.get("R_amb [mm]", 0.0))
        if R_amb_user <= 0.0:
            margin_deg = 5.0
            R_amb_mm = geom["R_exit"] + L_amb_mm * math.tan(
                math.radians(spray_halfangle + margin_deg))
        else:
            R_amb_mm = R_amb_user
        ambient = {"L_amb": L_amb_mm, "R_amb": R_amb_mm,
                   "spray_halfangle": spray_halfangle,
                   "topology": topology,
                   "cone_half_angle_deg": cone_half_angle,
                   "inner_amb_shrink": inner_amb_shrink}

    # Mesh/BC-Toggles (Lessons 028/029/038-v2). Defaults: BL refine 2x + Zhang BCs + no-prefill.
    bl_passes = max(0, min(2, int(float(j.get("BL refinement passes (0..2)", 2) or 0))))
    zhang_bcs = bool(int(float(j.get("Zhang outlet BCs (0/1)", 1) or 0)))
    disable_prefill = bool(int(float(j.get("Disable prefill (0/1)", 1) or 0)))

    info += [
        "=" * 60,
        "A3_nozzle -- Zhang 2023 pressure-swirl nozzle (interFoam VOF)",
        "=" * 60,
        "Project   : %s" % project,
        "DOE-Code  : %s" % (doe_note or "(none -- individual A3 params)"),
        "Geometry  : Di=%(Di)g Ds=%(Ds)g Ls=%(Ls)g alpha=%(alpha)g Do=%(Do)g "
        "Lo=%(Lo)g Lk=%(Lk)g n_inlet=%(n_inlet)g L_exit=%(L_exit)g R_exit=%(R_exit)g [mm,deg]"
        % geom,
        "Mesh      : N_Do=%(N_Do)d bl_cells=%(bl_cells)d (variant #6: O-grid + snappy inlets)"
        % mesh,
        "Flow      : U_inlet=%g m/s, endTime=%g s, mpirun -np %d" % (u_inlet, end_time, n_proc),
        "Ambient   : %s" % (
            "L_amb=%(L_amb)g mm  R_amb=%(R_amb)g mm  (spray cone-angle planes active)" % ambient
            if ambient else "off (no spray cone-angle; internal-flow only)"),
        "Objectives: film [mm] + air-core [mm] + U_theta [m/s] + Vordruck [MPa]"
        + (" + cone-angle [deg]" if ambient else ""),
        "=" * 60,
    ]
    update_progress(progress_path, info)

    if do_source_of:
        py_OF_utils.source_OF(source_version=OF_BASHRC)

    case_dir = Path(series_path) / ("case_%s" % project.replace(" ", "_"))
    png_path = os.path.join(progress_path, "A3_nozzle_convergence.png")

    if killed(series_path):
        info.append("\n>>> command_kill before build -- stopping.")
        update_progress(progress_path, info)
        touch(os.path.join(series_path, "command_finished"))
        return

    info.append("\n>>> Building case + variant-#6 mesh ...")
    update_progress(progress_path, info)
    build_case_dir(case_dir, geom=geom, mesh=mesh, ambient=ambient,
                   u_inlet=u_inlet, end_time=end_time, n_proc=n_proc)
    info.append("    case + mesh inputs written: %s" % case_dir.name)
    update_progress(progress_path, info)

    # Configuration overview snapshot (PNG + HTML in the case dir,
    # PNG also copied into progress/ for the streamlit live view).
    _write_overview(case_dir, progress_path, project=project,
                    doe_code=doe_code, geom=geom, mesh=mesh,
                    flow={"u_inlet": u_inlet, "end_time": end_time, "n_proc": n_proc},
                    toggles={"bl_passes": bl_passes, "zhang_bcs": int(zhang_bcs),
                              "disable_prefill": int(disable_prefill),
                              "enable_ambient": int(bool(ambient))},
                    ambient=ambient)
    info.append("    overview.png + overview.html written")
    update_progress(progress_path, info)

    summary = {"sweep": {"workflow": "A3_nozzle", "project": project,
                         "geometry": geom, "mesh": mesh,
                         "ambient": ambient, "u_inlet": u_inlet,
                         "end_time": end_time, "n_proc": n_proc},
               "cases": []}

    if do_run:
        if killed(series_path):
            info.append("\n>>> command_kill before solve -- stopping.")
            update_progress(progress_path, info)
            touch(os.path.join(series_path, "command_finished"))
            return
        info.append("\n>>> blockMesh -> snappyHexMesh -> checkMesh -> interFoam ...")
        if bl_passes > 0:
            info.append("    BL refinement: %dx refineWallLayer on (wall_straight wall_expansion) "
                        "+ outerWall split (wall_top buffer at z>-0.7mm to avoid Lesson-038 crash)"
                        % bl_passes)
        if zhang_bcs:
            info.append("    Zhang outlet BCs: p_rgh fixedValue 0 + alpha inletOutlet (Lesson 029)")
        if disable_prefill:
            info.append("    No prefill: alpha=0 everywhere, aircore builds naturally (Lesson 028)")
        update_progress(progress_path, info)
        mesh_and_solve(series_path, case_dir, n_proc, geom=geom,
                       bl_passes=bl_passes, zhang_bcs=zhang_bcs,
                       disable_prefill=disable_prefill)
        # Auto-Postprocessing: solver-Figuren (alpha/U/p) rendern. Die Mesh-Views
        # sind schon vor decomposePar entstanden -- hier kommen nur die Solver-
        # ergebnis-PNGs hinzu. Failures sind non-fatal (Ergebnisse bleiben).
        try:
            _run_pv_figure_renders(case_dir)
        except Exception as e:  # noqa: BLE001
            info.append("    pvbatch-figure-renders failed (non-fatal): %s" % e)
            update_progress(progress_path, info)

        cells = parse_cells(case_dir)
        z_exit_m = -(float(geom["Lo"]) + float(geom["Lk"]) + float(geom["L_exit"])) * 1e-3
        obj = extract_objectives(case_dir, Do_mm=geom["Do"], z_exit_m=z_exit_m)
        try:
            render_chart(png_path, obj, project)
        except Exception as e:  # noqa: BLE001
            print("chart render failed (non-fatal): %s" % e)
        _f = lambda v, p="%.3f": (p % v) if v is not None else "n/a"  # noqa: E731
        info += [
            "\n>>> Results:",
            "    cells        : %s" % (cells if cells is not None else "n/a"),
            "    film         : %s mm  (Zhang baseline ~0.65)" % _f(obj["film_mm"]),
            "    air-core dia : %s mm  (swirl-section diameter at orifice)" % _f(obj["air_core_mm"]),
            "    U_theta      : mean %s / max %s m/s  (circumferential velocity)"
            % (_f(obj["u_theta_mean"], "%.2f"), _f(obj["u_theta_max"], "%.2f")),
            "    Vordruck     : %s MPa" % _f(obj["vordruck_MPa"]),
        ]
        if ambient:
            info.append("    cone-angle   : %s deg  (mean ueber 5 Planes z=5..25mm)"
                        % _f(obj["cone_angle_deg"], "%.1f"))
        update_progress(progress_path, info)
        summary["cases"].append({
            "name": case_dir.name, "cells": cells,
            "film_mm": obj["film_mm"], "air_core_mm": obj["air_core_mm"],
            "u_theta_mean": obj["u_theta_mean"], "u_theta_max": obj["u_theta_max"],
            "vordruck_MPa": obj["vordruck_MPa"],
            "cone_angle_deg": obj["cone_angle_deg"],
            "cone_planes": obj["cone_planes"],
            **{k: geom[k] for k in GEOM_PARAMS}, **mesh,
            **({"L_amb": ambient["L_amb"], "R_amb": ambient["R_amb"]} if ambient else {}),
        })
    else:
        info.append("    (solver run skipped -- offline build only)")
        update_progress(progress_path, info)
        summary["cases"].append({"name": case_dir.name, "built": True})

    with open(os.path.join(progress_path, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    zip_results(series_path, results_path, case_dir)
    touch(os.path.join(series_path, "command_finished"))
    print("Finished! (contact=%s)" % (contact or "-"))


if __name__ == "__main__":
    series = check_args()
    # honour the run_workflow.py --dry-build convention (A2_DO_SOURCE_OF /
    # A2_DO_RUN = "0" -> build case + mesh inputs, skip OF source + solver).
    _do_source = os.environ.get("A2_DO_SOURCE_OF", "1") != "0"
    _do_run = os.environ.get("A2_DO_RUN", "1") != "0"
    try:
        run_worker(series, do_source_of=_do_source, do_run=_do_run)
    except Exception:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        try:
            prog = os.path.join(os.path.abspath(series), "progress")
            os.makedirs(prog, exist_ok=True)
            with open(os.path.join(prog, "9_error.txt"), "w", encoding="utf-8") as f:
                f.write("py_A3_nozzle.py failed:\n%s\n" % traceback.format_exc())
        except Exception:
            pass
        touch(os.path.join(os.path.abspath(series), "command_finished"))
        sys.exit(1)
