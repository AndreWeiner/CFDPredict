#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a complete OpenFOAM case for A2 brand_topology from params.json.

Minimal orchestrator (analogous to A2-mother's build_case.py, but slimmer):
    1. Load params.json.
    2. Build BrandTopology.
    3. Copy case_template/ to target/.
    4. Generate system/blockMeshDict via gen_blockmesh.
    5. If inner_arc_fractions non-empty: system/topoSetDict + createBafflesDict
       via gen_createbaffles.
    6. Substitute @U_INLET@, @K_INIT@, @OMEGA_INIT@, @NU@,
       @N_LAYER_SPLITS@, @LAYER_PATCHES@, @LAYER_THICKNESS@ in template files.
    7. Expand "vane.*" regex BC block to explicit vane{i}_master/vane{i}_slave.
    8. Strip frontAndBack BC if nz > 1 (3D -> walls).
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import stat
from pathlib import Path

from gen_blockmesh import BrandTopology, write_block_mesh_dict, _from_params_json
from gen_createbaffles import (render_toposet, render_createbaffles,
                                n_vanes as count_vanes)


_VANE_REGEX_BLOCK = re.compile(
    r'( *)"vane\.\*"\s*\n\1\{\n(.*?)\1\}\n',
    re.DOTALL,
)
_FB_BLOCK = re.compile(
    r'( *)frontAndBack\s*\n\1\{\n(.*?)\1\}\n',
    re.DOTALL,
)


def _substitute(text: str, mapping: dict[str, str]) -> str:
    for key, value in mapping.items():
        text = text.replace(f"@{key}@", value)
    return text


def _expand_vane_bcs(text: str, n_v: int) -> str:
    m = _VANE_REGEX_BLOCK.search(text)
    if not m:
        return text
    indent, body = m.group(1), m.group(2)
    if n_v == 0:
        return _VANE_REGEX_BLOCK.sub("", text)
    parts = []
    for i in range(1, n_v + 1):
        for side in ("master", "slave"):
            parts.append(f"{indent}vane{i}_{side}\n{indent}{{\n{body}{indent}}}\n")
    return _VANE_REGEX_BLOCK.sub("".join(parts), text)


def _strip_frontandback(text: str) -> str:
    return _FB_BLOCK.sub("", text)


def _layer_patches_str(n_v: int) -> str:
    """Patch list for refineWallLayer: walls + all vane patches."""
    parts = ["walls"]
    for i in range(1, n_v + 1):
        parts.append(f"vane{i}_master")
        parts.append(f"vane{i}_slave")
    return "(" + " ".join(parts) + ")"


def _inlet_state(geo: BrandTopology, Re: float, nu: float, ti: float) -> tuple[float, float, float]:
    """Return (U_inlet, k_init, omega_init) for the k-omega-SST BCs."""
    U = Re * nu / geo.H
    k = 1.5 * (ti * U) ** 2                  # standard inlet k from intensity
    # omega from length scale = 0.07*H (typical):
    Lt = 0.07 * geo.H
    Cmu = 0.09
    omega = (k ** 0.5) / (Cmu ** 0.25 * Lt)
    return U, k, omega


def build(target: Path, geo: BrandTopology,
          Re: float, nu: float, turbulent_intensity: float,
          n_layer_splits: int = 0, layer_thickness: float = 0.5) -> dict:
    template = Path(__file__).parent / "case_template"
    if not template.is_dir():
        raise FileNotFoundError(f"case_template not found at {template}")

    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(template, target)

    U, k, omega = _inlet_state(geo, Re, nu, turbulent_intensity)
    n_v = count_vanes(geo)

    mapping = {
        "U_INLET":         f"{U:.6g}",
        "K_INIT":          f"{k:.6g}",
        "OMEGA_INIT":      f"{omega:.6g}",
        "NU":              f"{nu:.6g}",
        "N_LAYER_SPLITS":  str(n_layer_splits),
        "LAYER_PATCHES":   _layer_patches_str(n_v),
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

    # Geometry dicts
    write_block_mesh_dict(geo, target / "system" / "blockMeshDict")
    if n_v > 0:
        (target / "system" / "topoSetDict").write_text(
            render_toposet(geo), encoding="utf-8", newline="\n")
        (target / "system" / "createBafflesDict").write_text(
            render_createbaffles(geo), encoding="utf-8", newline="\n")

    # BC expansion
    is_3d = geo.nz > 1
    for field_name in ("U", "p", "k", "omega", "nut"):
        path = target / "0" / field_name
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        new_text = _expand_vane_bcs(text, n_v)
        if is_3d:
            new_text = _strip_frontandback(new_text)
        if new_text != text:
            path.write_text(new_text, encoding="utf-8", newline="\n")

    return {
        "U_inlet": U, "k_init": k, "omega_init": omega,
        "Re": Re, "nu": nu, "n_vanes": n_v,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target", type=Path, required=True)
    ap.add_argument("--params", default="params.json")
    args = ap.parse_args()

    with open(args.params, "r", encoding="utf-8") as f:
        spec = json.load(f)
    p = spec["parameters"]
    def v(name): return p[name]["value"]

    geo = _from_params_json(args.params)
    info = build(args.target, geo,
                 Re=v("Re"), nu=v("nu"),
                 turbulent_intensity=v("turbulent_intensity"),
                 n_layer_splits=v("n_layer_splits"),
                 layer_thickness=v("layer_thickness"))

    print(f"built case at {args.target}")
    print(f"  params:   {args.params}")
    print(f"  geometry: H={geo.H}, W={geo.W}, R={geo.R}, "
          f"L_in={geo.L_in}, L_trail={geo.L_trail}, L_out={geo.L_out}")
    print(f"  vanes:    {info['n_vanes']} "
          f"(from inner_arc_fractions = {geo.inner_arc_fractions})")
    print(f"  flow:     Re={info['Re']:.3g}, nu={info['nu']:.3g}")
    print(f"            U_inlet={info['U_inlet']:.4f} m/s")
    print(f"            k_init={info['k_init']:.4e}, "
          f"omega_init={info['omega_init']:.4e}")
    print(f"  run with: cd {args.target} && ./Allrun")


if __name__ == "__main__":
    main()
