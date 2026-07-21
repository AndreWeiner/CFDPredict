#!/usr/bin/env python3
"""Build a complete OpenFOAM case for A2 from params.json + case_template/.

Workflow:
    1. Load params.json (or alternative file via --params).
    2. Apply optional CLI overrides (each scalar param has a matching --name flag).
    3. Build ElbowGeometry + FlowParams from the spec.
    4. Copy case_template/ to target/.
    5. Generate blockMeshDict, optionally topoSetDict + createBafflesDict
       (if vane_r* parameters are present).
    6. Substitute @U_INLET@, @K_INIT@, @OMEGA_INIT@, @NU@ in template files.
    7. Expand "vane.*" regex BC block to explicit vaneN_master / vaneN_slave entries.

Quasi-2D, k-omega-SST, simpleFoam. Optional guide vanes via createBaffles.
"""
from __future__ import annotations

import argparse
import re
import shutil
import stat
from pathlib import Path

from gen_blockmesh import ElbowGeometry, render as render_blockmesh
from gen_createbaffles import render as render_createbaffles, render_toposet
from gen_refine_layer import render_patch_arg as render_layer_patches
import params
from params import FlowParams


DEFAULT_PARAMS = Path(__file__).parent / "params.json"


def _substitute(text: str, mapping: dict[str, str]) -> str:
    for key, value in mapping.items():
        text = text.replace(f"@{key}@", value)
    return text


_VANE_REGEX_BLOCK = re.compile(
    r'( *)"vane\.\*"\s*\n\1\{\n(.*?)\1\}\n',
    re.DOTALL,
)
_FB_BLOCK = re.compile(
    r'( *)frontAndBack\s*\n\1\{\n(.*?)\1\}\n',
    re.DOTALL,
)


def _strip_frontandback(text: str) -> str:
    """For 3D cases (nz>1), the frontAndBack patch is merged into 'walls', so
    the BC entry must be removed."""
    return _FB_BLOCK.sub("", text)


def _expand_vane_bcs(text: str, n_vanes: int) -> str:
    """Expand a ``"vane.*"`` regex BC block into explicit per-baffle entries.

    createBaffles writes ``type calculated`` entries for new patches into existing
    field files, and explicit patch names take precedence over regex matches.
    By pre-writing explicit entries here (before createBaffles runs), we keep
    the proper wall BCs instead of getting overwritten.
    """
    m = _VANE_REGEX_BLOCK.search(text)
    if not m:
        return text
    indent, body = m.group(1), m.group(2)
    if n_vanes == 0:
        return _VANE_REGEX_BLOCK.sub("", text)
    parts = []
    for i in range(1, n_vanes + 1):
        for side in ("master", "slave"):
            parts.append(f"{indent}vane{i}_{side}\n{indent}{{\n{body}{indent}}}\n")
    return _VANE_REGEX_BLOCK.sub("".join(parts), text)


def build(target: Path, g: ElbowGeometry, f: FlowParams,
          n_layer_splits: int = 0, layer_thickness: float = 0.5) -> dict[str, float]:
    template = Path(__file__).parent / "case_template"
    if not template.is_dir():
        raise FileNotFoundError(f"case_template not found at {template}")

    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(template, target)

    U, k, omega = f.inlet_state(g.H)
    mapping = {
        "U_INLET":         f"{U:.6g}",
        "K_INIT":          f"{k:.6g}",
        "OMEGA_INIT":      f"{omega:.6g}",
        "NU":              f"{f.nu:.6g}",
        "N_LAYER_SPLITS":  str(n_layer_splits),
        "LAYER_PATCHES":   render_layer_patches(g),
        "LAYER_THICKNESS": f"{layer_thickness:g}",
    }

    for path in target.rglob("*"):
        if not path.is_file():
            continue
        if path.name in ("Allrun", "Allclean"):
            path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if "@" in text:
            path.write_text(_substitute(text, mapping), encoding="utf-8", newline="\n")

    (target / "system" / "blockMeshDict").write_text(
        render_blockmesh(g), encoding="utf-8", newline="\n"
    )
    if g.vane_radii:
        (target / "system" / "topoSetDict").write_text(
            render_toposet(g), encoding="utf-8", newline="\n"
        )
        (target / "system" / "createBafflesDict").write_text(
            render_createbaffles(g), encoding="utf-8", newline="\n"
        )


    n_vanes = len(g.vane_radii)
    is_3d = g.nz > 1
    for field_name in ("U", "p", "k", "omega", "nut"):
        path = target / "0" / field_name
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        new_text = _expand_vane_bcs(text, n_vanes)
        if is_3d:
            new_text = _strip_frontandback(new_text)
        if new_text != text:
            path.write_text(new_text, encoding="utf-8", newline="\n")

    return {"U_inlet": U, "k_init": k, "omega_init": omega, "Re": f.Re, "nu": f.nu}


_SCALAR_OVERRIDES = [
    ("H", float), ("W", float), ("R", float),
    ("L_in", float), ("L_out", float),
    ("nx_in", int), ("nx_bend", int), ("nx_out", int),
    ("ny", int), ("nz", int),
    ("Re", float), ("nu", float),
    ("turbulent_intensity", float),
    ("n_layer_splits", int), ("layer_thickness", float),
]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target", type=Path, required=True)
    p.add_argument("--params", type=Path, default=DEFAULT_PARAMS,
                   help=f"JSON params file (default: {DEFAULT_PARAMS.name})")

    ov = p.add_argument_group("overrides", "replace values from params.json")
    for name, t in _SCALAR_OVERRIDES:
        ov.add_argument(f"--{name.replace('_', '-')}", type=t, default=None, dest=name)
    ov.add_argument("--vane-radii", type=float, nargs="+", default=None,
                    dest="vane_radii",
                    help="override vane_r1, vane_r2, ... positionally")
    p.add_argument("--auto-layers", action="store_true",
                   help="n_layer_splits Re-adaptiv aus y+_target berechnen "
                        "(wall_resolution.py) statt aus params.json")
    p.add_argument("--y-plus-target", type=float, default=1.0,
                   help="Ziel-y+ für --auto-layers (default 1.0)")

    args = p.parse_args()

    spec = params.load(args.params)
    params.apply_overrides(spec, {name: getattr(args, name)
                                  for name, _ in _SCALAR_OVERRIDES})
    if args.vane_radii is not None:
        params.set_vane_radii(spec, args.vane_radii)

    g = params.to_geometry(spec)
    f = params.to_flow(spec)
    layer_thickness = float(spec.value("layer_thickness"))
    from wall_resolution import layer_splits_for_Re, predicted_yplus
    # Default: params.json honorieren (auto_layers/y_plus_target, sonst fixes n_layer_splits).
    n_layer_splits = params.resolve_layers(spec, g, f)
    if args.n_layer_splits is not None:
        n_layer_splits = args.n_layer_splits                 # expliziter CLI-Fix-Override
    elif args.auto_layers:
        n_layer_splits = layer_splits_for_Re(f.Re, f.nu, g.H, # CLI-Recompute @ Ziel-y+
                                             args.y_plus_target, g.ny)
    if n_layer_splits > 0:
        yp = predicted_yplus(f.Re, f.nu, g.H, g.ny, n_layer_splits)
        print(f"  layers (Re-adaptiv): Re={f.Re:.3g} -> {n_layer_splits} "
              f"refineWallLayer-Pass(es), y+_pred~{yp:.1f} "
              f"(nut=nutUSpalding -> yPlus-FO meldet echtes y+).")
    info = build(args.target, g, f,
                 n_layer_splits=n_layer_splits, layer_thickness=layer_thickness)

    print(f"built case at {args.target}")
    print(f"  params:   {args.params}")
    print(f"  geometry: H={g.H}, R={g.R}, L_in={g.L_in}, L_out={g.L_out}")
    if g.vane_radii:
        print(f"  vanes:    {len(g.vane_radii)} at radii "
              f"{sorted(g.vane_radii, reverse=True)}")
    n_bands = len(g.vane_radii) + 1
    ny_per_band = max(1, g.ny // n_bands)
    total = n_bands * (g.nx_in + g.nx_bend + g.nx_out) * ny_per_band * g.nz
    print(f"  mesh:     {n_bands} band(s), {ny_per_band} cells across each, {total} cells total")
    print(f"  flow:     Re={info['Re']:.3g}, nu={info['nu']:.3g}")
    print(f"            U_inlet={info['U_inlet']:.4f} m/s")
    print(f"            k_init={info['k_init']:.4e}, omega_init={info['omega_init']:.4e}")
    if n_layer_splits > 0:
        print(f"  layers:   {n_layer_splits} refineWallLayer pass(es) "
              f"@ thickness={layer_thickness} -> wall cell x {layer_thickness**n_layer_splits:.3f}")
    print(f"  run with: cd {args.target} && ./Allrun")


if __name__ == "__main__":
    main()
