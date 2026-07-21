#!/usr/bin/env python3
"""Re-abhängige Wandauflösung — hält y⁺ über den Re-Bereich konstant.

Hintergrund (Dresden/Andre Weiner): Damit Diskretisierungsfehler Re-unabhängig
bleiben (Vergleichbarkeit/Universalität → ROM-Transfer), wird die wandnahe
Auflösung an Re gekoppelt, sodass y⁺ ≈ konstant bleibt — statt das Gitter fix zu
lassen (dann schwankt y⁺ über den Re-Bereich um Faktor ~60).

Etabliertes Verfahren ("y⁺-Rechner"):
    C_f(Re)  →  u_τ = U·√(C_f/2)  →  Δy_zentrum = y⁺·ν/u_τ
    erste ZellHÖHE = 2·Δy_zentrum  (y⁺ wird am ERSTEN ZELLZENTRUM ausgewertet)

Umsetzung in unserer Pipeline: über `refineWallLayer` (jeder Pass halbiert die
wandadjazente Zelle bei thickness 0.5) — deterministisch, Bayes-/ROM-safe.
`n_layer_splits` wird aus (Re, y⁺_target, ny) berechnet.

⚠ y⁺ vs. y*: Der OpenFOAM-`yPlus`-FO meldet bei **nutkWallFunction** y* (k-basiert,
C_μ^¼√k), NICHT das τ_w-basierte y⁺, das diese C_f-Korrelation schätzt. Beide sind
im Gleichgewichts-Log-Bereich gleich, divergieren aber an Ablösung/Staupunkt.
→ Für konsistente Verifikation im Re-Sweep auf **nutUSpaldingWallFunction**
umstellen (FO meldet dann echtes y⁺) oder y⁺≈1 aufgelöst rechnen.

Korrelation: C_f ≈ 0.058·Re⁻⁰·² (flat-plate, Engineering-Schätzung). Der
`yPlus`-FO ist die Ground-Truth zur Nachkorrektur.
"""
from __future__ import annotations

import argparse
import math

CMU = 0.09


def skin_friction(Re: float) -> float:
    """Reibungsbeiwert C_f (flat-plate-Korrelation, Engineering-Schätzung)."""
    return 0.058 * Re ** -0.2


def u_tau(Re: float, nu: float, H: float) -> float:
    U = Re * nu / H
    return U * math.sqrt(skin_friction(Re) / 2.0)


def first_cell_height_for_yplus(Re: float, nu: float, H: float,
                                y_plus: float) -> float:
    """Wandadjazente Zell**höhe**, sodass y⁺ am Zellzentrum = y_plus.

    y⁺ wird am ersten Zellzentrum ausgewertet → Zentrumsabstand = höhe/2,
    daher höhe = 2 · (y⁺·ν/u_τ)."""
    y_centre = y_plus * nu / u_tau(Re, nu, H)
    return 2.0 * y_centre


def predicted_yplus(Re: float, nu: float, H: float, ny: int,
                    n_splits: int) -> float:
    """y⁺ (am Zellzentrum), das ny + n_splits refineWallLayer-Pässe liefern."""
    cell_h = (H / ny) / (2 ** n_splits)
    y_centre = cell_h / 2.0
    return y_centre * u_tau(Re, nu, H) / nu


def layer_splits_for_Re(Re: float, nu: float, H: float, y_plus: float,
                        ny: int, max_splits: int = 12) -> int:
    """Anzahl refineWallLayer-Pässe, damit y⁺ ≲ y_plus (Zweierpotenz-Granularität).

    base_cell = H/ny (uniforme blockMesh-Zellhöhe quer zum Kanal)."""
    base_cell = H / ny
    target_h = first_cell_height_for_yplus(Re, nu, H, y_plus)
    if target_h >= base_cell:
        return 0
    n = math.ceil(math.log2(base_cell / target_h))
    return max(0, min(n, max_splits))


def _table(re_list, nu, H, y_plus, ny):
    print(f"# Re-adaptive Wandauflösung — Ziel y+={y_plus}, ny={ny}, H={H}, nu={nu}")
    print(f"# C_f ~ 0.058*Re^-0.2 ; base cell H/ny = {H/ny:.4g}")
    print(f"{'Re':>9} {'U[m/s]':>8} {'u_tau':>9} {'dy_cell':>10} "
          f"{'splits':>7} {'y+_pred':>8}")
    for Re in re_list:
        U = Re * nu / H
        ut = u_tau(Re, nu, H)
        n = layer_splits_for_Re(Re, nu, H, y_plus, ny)
        yp = predicted_yplus(Re, nu, H, ny, n)
        dy = first_cell_height_for_yplus(Re, nu, H, y_plus)
        print(f"{Re:9.0e} {U:8.3f} {ut:9.4f} {dy:10.2e} {n:7d} {yp:8.1f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--y-plus", type=float, default=1.0, help="Ziel-y+ (default 1.0)")
    ap.add_argument("--ny", type=int, default=20, help="Querzellen über H (default 20)")
    ap.add_argument("--H", type=float, default=1.0)
    ap.add_argument("--nu", type=float, default=1.5e-5)
    ap.add_argument("--re", type=float, nargs="+",
                    default=[1e4, 3e4, 1e5, 3e5, 1e6])
    args = ap.parse_args()
    _table(args.re, args.nu, args.H, args.y_plus, args.ny)


if __name__ == "__main__":
    main()
