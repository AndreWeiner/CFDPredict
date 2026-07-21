#!/usr/bin/env python3
"""Parameter loading + introspection for A2 cases.

Reads params.json (schema-rich values + bounds + roles), produces
ElbowGeometry and FlowParams for build_case.py, and exposes the
design-variable subset for later Bayes-optimization integration
(Ax / BoTorch in Iteration 3).
"""
from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from gen_blockmesh import ElbowGeometry


@dataclass
class FlowParams:
    Re: float = 1.0e5
    nu: float = 1.5e-5
    turbulent_intensity: float = 0.05
    mixing_length_factor: float = 0.07

    def inlet_state(self, H: float) -> tuple[float, float, float]:
        U = self.Re * self.nu / H
        k = 1.5 * (self.turbulent_intensity * U) ** 2
        L = self.mixing_length_factor * H
        omega = math.sqrt(k) / (0.09 ** 0.25 * L)
        return U, k, omega


@dataclass
class ParamSpec:
    case_type: str
    version: int
    parameters: dict[str, dict[str, Any]]
    constraints: list[str] = field(default_factory=list)
    objectives: list[dict[str, Any]] = field(default_factory=list)

    def value(self, name: str) -> Any:
        return self.parameters[name]["value"]

    def set_value(self, name: str, value: Any) -> None:
        if name not in self.parameters:
            raise KeyError(f"unknown parameter {name!r}")
        self.parameters[name]["value"] = value

    def design_params(self) -> dict[str, dict[str, Any]]:
        return {n: p for n, p in self.parameters.items() if p.get("role") == "design"}

    def fixed_params(self) -> dict[str, dict[str, Any]]:
        return {n: p for n, p in self.parameters.items() if p.get("role") == "fixed"}


def load(path: Path) -> ParamSpec:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return ParamSpec(
        case_type=raw.get("case_type", ""),
        version=raw.get("version", 1),
        parameters=raw["parameters"],
        constraints=raw.get("constraints", []),
        objectives=raw.get("objectives", []),
    )


def dump(spec: ParamSpec, path: Path) -> None:
    Path(path).write_text(
        json.dumps(
            {
                "version": spec.version,
                "case_type": spec.case_type,
                "parameters": spec.parameters,
                "constraints": spec.constraints,
                "objectives": spec.objectives,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
        newline="\n",
    )


_VANE_KEY = re.compile(r"^vane_r(\d+)$")


def _vane_keys(spec: ParamSpec) -> list[str]:
    return sorted(
        (k for k in spec.parameters if _VANE_KEY.match(k)),
        key=lambda k: int(_VANE_KEY.match(k).group(1)),
    )


def _vane_ext(spec: ParamSpec) -> list[float]:
    """Downstream-Verlängerung je Blech aus vane_ext{N}-Keys (default 0),
    ausgerichtet auf die vane_r{N}-Reihenfolge. Rückwärtskompatibel: fehlen die
    Keys, sind alle Verlängerungen 0 (unmodifizierte Full-Arc-Bleche)."""
    exts = []
    for k in _vane_keys(spec):
        idx = _VANE_KEY.match(k).group(1)
        ek = f"vane_ext{idx}"
        exts.append(float(spec.value(ek)) if ek in spec.parameters else 0.0)
    return exts


def to_geometry(spec: ParamSpec) -> ElbowGeometry:
    vane_radii = [spec.value(k) for k in _vane_keys(spec)]
    return ElbowGeometry(
        H=spec.value("H"),
        W=spec.value("W"),
        R=spec.value("R"),
        L_in=spec.value("L_in"),
        L_out=spec.value("L_out"),
        nx_in=spec.value("nx_in"),
        nx_bend=spec.value("nx_bend"),
        nx_out=spec.value("nx_out"),
        ny=spec.value("ny"),
        nz=spec.value("nz"),
        vane_radii=vane_radii,
        vane_ext=_vane_ext(spec),
    )


def to_flow(spec: ParamSpec) -> FlowParams:
    return FlowParams(
        Re=spec.value("Re"),
        nu=spec.value("nu"),
        turbulent_intensity=spec.value("turbulent_intensity"),
    )


def resolve_layers(spec: ParamSpec, g: ElbowGeometry, f: FlowParams) -> int:
    """Anzahl refineWallLayer-Paesse, die build() bekommt.

    Default (auto_layers=true): Re-adaptiv aus (Re, y_plus_target) berechnet, sodass
    y+ ueber den Re-Bereich konstant bleibt (Dresden-Anforderung). Andernfalls der
    feste n_layer_splits-Wert aus params.json. Beides rueckwaertskompatibel: fehlen
    die neuen Keys, faellt es auf n_layer_splits zurueck."""
    auto = bool(spec.value("auto_layers")) if "auto_layers" in spec.parameters else False
    if not auto:
        return int(spec.value("n_layer_splits"))
    from wall_resolution import layer_splits_for_Re
    y_plus = float(spec.value("y_plus_target")) if "y_plus_target" in spec.parameters else 1.0
    return layer_splits_for_Re(f.Re, f.nu, g.H, y_plus, g.ny)


def apply_overrides(spec: ParamSpec, overrides: dict[str, Any]) -> None:
    """Apply scalar overrides in-place. None values are ignored."""
    for name, value in overrides.items():
        if value is not None and name in spec.parameters:
            spec.set_value(name, value)


def set_vane_radii(spec: ParamSpec, radii: list[float]) -> None:
    """Replace vane_r1, vane_r2, ... with the given list (in order)."""
    keys = _vane_keys(spec)
    if len(keys) != len(radii):
        raise ValueError(
            f"params.json has {len(keys)} vane_r* slot(s), got {len(radii)} radii"
        )
    for k, r in zip(keys, radii):
        spec.set_value(k, r)
