# -*- coding: utf-8 -*-
"""
A3 SwirlNozzleInjector -- zentrale Parametrisierung (Single Source of Truth).

Wird von sketch.py (Querschnitt, xy) und sketch_side.py (z-Aufbau, x-z) importiert.

Offizielle Geometrie nach Zhang et al. 2023, Fig. 1(b) "Scale of pressure swirl
nozzle". Alle inneren O-Grid-Masze werden als feste VERHAELTNISSE zu diesen
offiziellen Parametern abgeleitet. Primaere Laengenskala: Do (Orifice-Durchmesser).

Zhang Fig.1(b) Bezeichner
-------------------------
    Di      Tangentialeinlass-Durchmesser (n_inlet Stueck, rund)
    Ds      Drallkammer-Durchmesser
    Ls      Drallkammer-Laenge
    alpha   Kontraktions-Halbwinkel (von der Horizontalen),  Ds -> Do
    Do      Orifice-Durchmesser           [PRIMAERE LAENGENSKALA]
    Lo      Orifice-Laenge
    Lk      Expansions-Laenge
    theta   Expansions-VOLLwinkel  (Fig.1b: theta; 0 = Baseline gerade, >0 divergent)

Einheiten: mm, Grad.
"""
import math

# ======================================================================
# 1) OFFIZIELLE Parameter (Zhang 2023, Table 1 Baseline)
# ======================================================================
Di      = 2.0
Ds      = 10.0
Ls      = 6.0
alpha   = 45.0       # Kontraktions-Halbwinkel
Do      = 5.0        # <<< PRIMAERE LAENGENSKALA
Lo      = 20.0
Lk      = 6.0        # Expansions-Laenge (user L3 = 6 mm; Zhang Baseline 5.0)
n_inlet = 4

# ======================================================================
# 2) PHYSIK fuer BL-/Film-Aufloesung
# ======================================================================
film_thickness = 0.65    # mm, erwartete Fluessigfilmdicke am Austritt (Zhang Baseline)
film_cover     = 1.30    # BL-Band deckt das film_cover-fache der Filmdicke ab (>1)

# fuer y+-Schaetzung (Rohrstroemung, Blasius)
U_ref   = 18.46          # m/s  (Zhang velocity-inlet)
nu_ref  = 1.0e-6         # m^2/s (Wasser, 20 C)

# ======================================================================
# 3) MESH-Aufloesung
# ======================================================================
N_Do         = 40        # Basis-Zellen ueber Do (Grundgitter-Aufloesung)
bl_cells     = 12        # GLEICHMAESSIGE radiale Zellen im wandnahen BL-/Film-Band
                         #   (= Ring2, p4->p8); deckt den Film uniform ab.
                         #   12 -> 24 ausprobiert 2026-06-03: erzeugte 832 sevNonOrtho
                         #   Faces nach snappy (vs 0 bei bl_cells=12). Verworfen.
                         #   Film im Ambient ist mit bl_cells=12 nur ~1.5 Cells dick;
                         #   bessere Optionen: snappy addLayers mit firstLayerThickness
                         #   oder gezieltes refineWallLayer auf Frustum-Mantel-cellSet.
                         # TODO 2026-06-04 (Martin): Auch am aufweitenden Nozzle-Exit
                         #   (Lampenschirm, Expansion-Segment Lk) ist der Film nur knapp
                         #   mit Cells unterlegt - hier z-Cells am Auslass splitten oder
                         #   z-grading in Expansion-Segment einfuehren, sodass die Wand-
                         #   region feinere Cells bekommt waehrend der Film duenner wird.
trans_ratio  = 0.60      # Ring1 (Uebergang Quadrat->Kreis1) als Anteil von t2(BL)
buffer_ratio = 0.40      # t3 (Kreis2->zentr. Quadrat, Aussen-O-Grid) als Anteil von t2

# ======================================================================
# 4) z-AUFBAU (Zhang-Profil, Stroemung aufwaerts; z=0 am Orifice-Eintritt)
#    Drallkammer -> Kontraktion(alpha) -> Orifice(Lo) -> Expansion(Lk,theta) -> Exit
#    Tangentialeinlaesse (Di, n_inlet) werden NICHT in blockMesh definiert,
#    sondern von snappyHexMesh aus der Drallkammerwand gestanzt.
# ======================================================================
L_exit = 0.0             # Gerade nach der Expansion -- NICHT im Zhang-Paper (Fig.1b
                         # kennt nur Di/Ds/Ls/alpha/Do/Lo/Lk/theta). Default 0 = paper-
                         # konform; >0 = Martin-User-Variante mit zusaetzlichem L4-Stueck.
R_exit = 5.0             # Bohrungsradius am Expansionsende        (user 0.005 m); treibt theta

# ======================================================================
# 4b) AMBIENT-Erweiterung fuer Sprueh-Kegelwinkel (Task 3, 2026-06-01)
# ======================================================================
# Wenn AMBIENT=True: stromabwaerts des Nozzle-Exits 2 zusaetzliche z-Stationen
# (amb_mid, amb_far) mit aufgeweitetem Scale -> die O-Grid-Frusta oeffnen sich in
# den Spray-Konus; die alpha=0.25-Kontur kann auf Plane-Plots ausgewertet werden.
# Aussen-Raster der Ambient-Segmente + Top der amb_far landen im 'outlet'-Patch
# (atmosphaerische pressureInletOutletVelocity-BC, identisch zur Exit-BC).
AMBIENT      = False         # Default off (alte Cases unveraendert)
L_amb        = 30.0          # axiale Ambient-Laenge nach exit_top [mm]
R_amb        = 17.5          # radialer Ambient-Radius [mm] -- AUTO bei AMBIENT
                              # (siehe SPRAY_HALFANGLE_DEG); manuell nur wenn
                              # A3_R_amb env explicit gesetzt
amb_mid_frac = 0.5           # z-Position der Zwischenstation als Anteil von L_amb
amb_mid_frac_s = 0.5         # Scale-Anteil zw. R_exit-Scale und R_amb-Scale

# Spray-Cone-Geometrie fuer Auto-Mesh-Auffaecherung im Ambient:
# Der Frustum-Mantel (Bohrungs-Exit -> Ambient-Aussen) muss den vollen Spray-
# Cone abdecken + Marge -- sonst laufen Lamelle + Primaertropfen aus dem fein
# aufgeloesten Mesh-Bereich raus (visuell sichtbar in der 102_PEGASUS-Sim:
# Tropfen treffen die Aussenwand-Cells bei groberem Mesh). Default 47 deg =
# Zhang Tab.4 A1B2C2D2 Cone-Half-Angle (94.4°/2). Pro DOE-Variante anders --
# manuell ueberschreibbar via env A3_SPRAY_HALFANGLE_DEG.
SPRAY_HALFANGLE_DEG = 47.0   # Halbwinkel des Sprueh-Cones [deg]
SPRAY_MARGIN_DEG    = 5.0    # zusaetzliche Marge gegen Wobble/Spread

# ======================================================================
# 4c) TOPOLOGIE-VARIANTE
# ======================================================================
# v1 (Legacy ohne AMBIENT): uniform 21-Block-Topologie ueber alle z-Stationen,
#     Outer-Raster wird im Orifice/Expansion-Bereich von snappy weggeschnitten
#     (-> Cell-Verschwendung im blockMesh, kein Solver-Impact).
# v2 (NEUER DEFAULT bei AMBIENT=1, Lesson 103_PEGASUS): stations-abhaengige
#     Topologie. Inner-O-Grid + Ring3 faechern im Ambient frustum-artig auf,
#     Outer-Raster bleibt konstant. Loest die Frustum-Knick-Geometrie deutlich
#     besser auf (102_ max non-orth 78 deg -> 103_ 65 deg, avg 22 deg -> 16 deg)
#     und haelt Lamelle + Primaertropfen im fein aufgeloesten Bereich.
#       Drallkammer + Kontraktion: Core+Ring1+Ring2+Trans+Out (21)
#       Orifice + Expansion:       Core+Ring1+Ring2            (9)
#       Exit + Ambient:            +Ring3 (R3, konstant)       (13)
TOPOLOGY     = "v1"          # OHNE AMBIENT bleibt v1; mit AMBIENT siehe Auto-Switch unten
R3           = 10.0          # Ring3 Radius [mm] -- nur TOPOLOGY=v2 (= 2*R_exit)
cone_half_angle_deg = 40.0   # Ring3-Frustum half-angle [deg] -- 40 deg deckt
                              # Spray-Cone 47 deg + 5 deg Marge (103_PEGASUS); 60 deg
                              # war Ueberweite, Aussenrand zu weit weg.
inner_amb_shrink    = 0.9    # Inner-O-Grid-Aussenrand (Ring2) wird im Ambient
                              # leicht nach innen gezogen (~10%) damit Ring3 Pufferzone
                              # bekommt + spitzer Block-Vertex am Frustum-Knick
                              # entschaerft wird.

# --- DOE / Sweep-Overrides: Umgebungsvariablen A3_<Param> ueberschreiben die
#     Baseline oben (additiv; ohne gesetzte Var bleibt das Verhalten exakt
#     gleich). Erlaubt Parameter-Sweeps ohne diese Datei zu editieren. ---
import os as _os
import math as _math
_INT_OVERRIDES = ("n_inlet", "N_Do", "bl_cells")
_R_AMB_USER_OVERRIDE = _os.environ.get("A3_R_amb") is not None
_TOPOLOGY_USER_OVERRIDE = _os.environ.get("A3_TOPOLOGY") is not None
for _k in ("Di", "Ds", "Ls", "alpha", "Do", "Lo", "Lk", "n_inlet", "L_exit", "R_exit",
           "N_Do", "bl_cells", "L_amb", "R_amb", "R3", "cone_half_angle_deg",
           "SPRAY_HALFANGLE_DEG", "SPRAY_MARGIN_DEG", "inner_amb_shrink"):
    _v = _os.environ.get("A3_" + _k)
    if _v is not None:
        globals()[_k] = int(_v) if _k in _INT_OVERRIDES else float(_v)
_v = _os.environ.get("A3_TOPOLOGY")
if _v is not None:
    TOPOLOGY = _v
_v = _os.environ.get("A3_AMBIENT")
if _v is not None:
    AMBIENT = _v not in ("0", "false", "False", "")

# Auto-TOPOLOGY: bei AMBIENT ist v2 der validierte Default (Lesson 103_PEGASUS);
# expliziter A3_TOPOLOGY-Env-Override bleibt erhalten.
if AMBIENT and not _TOPOLOGY_USER_OVERRIDE:
    TOPOLOGY = "v2"

# Auto-R_amb: bei AMBIENT nehmen wir den Spray-Cone als Massgabe fuer den
# Frustum-Aussenrand (vs. festen Default 17.5 mm). Manueller Override bleibt
# erhalten wenn A3_R_amb env explicit gesetzt war.
if AMBIENT and not _R_AMB_USER_OVERRIDE:
    R_amb = R_exit + L_amb * _math.tan(_math.radians(
        SPRAY_HALFANGLE_DEG + SPRAY_MARGIN_DEG))


# ======================================================================
# ABGELEITETE Querschnitts-Geometrie (O-Grid, radial von Wand nach innen)
# ======================================================================
R_wall    = Do / 2.0                    # Bohrungswand = Kreis 2 = p8..p11
h_base    = Do / N_Do                   # Basis-Zellgroesse
t2        = film_cover * film_thickness # BL-/Film-Band (Ring2) radiale Dicke
t1        = trans_ratio * t2            # Uebergangsband (Ring1)
R1        = R_wall - t2                 # Kreis 1
R0_corner = R1 - t1                     # Quadrat-Eckradius (Diagonale)
a         = R0_corner / math.sqrt(2.0)  # Quadrat-Halbseite (Schale 0)
t3        = buffer_ratio * t2           # Kreis 2 -> zentrales Aussen-Quadrat
c         = R_wall + t3                 # zentrales Aussen-Quadrat Halbseite

h_wall    = t2 / bl_cells               # gleichmaessige wandnormale Zellhoehe im BL-Band
N_core    = max(2, round(2.0 * a / h_base))
n_r2      = bl_cells                    # BL-Band radiale Zellen (uniform)
# Uebergangsband radial 2.5x feiner als h_base nahelegen wuerde. ParaView-
# Diagnose chamberBottom-Slice (2026-06-03): Trans-Cells aktuell viel
# groesser als BL-Band-Cells, erzeugen visible Aspect-Sprung; deutet auf
# Auflösungs-Diskontinuitaet hin, an der die VOF-Front blutet (1.5x reichte
# nicht).  Faktor 2.5 -> ~10 cells, h_r1 ~0.05mm = feiner als h_wall.
n_r1      = max(1, round(2.5 * t1 / h_base))

# z-Profil der Bohrungswand (Radius je Abschnitt)
Rs = Ds / 2.0                                          # Drallkammer-Radius
Lc = (Rs - R_wall) * math.tan(math.radians(alpha))     # Kontraktionslaenge aus alpha
                                                       #   (alpha von Horizontale: Lz = dr*tan(alpha))
R_outer = 3.0 * c                                       # Hintergrund-/Aussenrand (3x3-Raster), konstant in z

# Expansions-Vollwinkel aus Zielradius R_exit (theta = Apex/Vollwinkel von der Achse):
theta = (math.degrees(2.0 * math.atan((R_exit - R_wall) / Lk))
         if R_exit > R_wall else 0.0)

# z-Grenzen (z=0 = Orifice-Eintritt)
z_swirl_bot = -(Lc + Ls)     # Boden Drallkammer (Wand; Zulauf via snappy)
z_contr_bot = -Lc            # Drallkammer-Oberkante = Kontraktions-Unterkante
z_orif_bot  = 0.0
z_orif_top  = Lo
z_exp_top   = Lo + Lk
z_exit_top  = Lo + Lk + L_exit

# Ambient-Stationen liegen oberhalb von z_exit_top (im original Frame); nach FLIP
# (180 deg um Y) landen sie stromabwaerts bei z = -(z_exit_top + L_amb).
z_amb_mid   = z_exit_top + amb_mid_frac * L_amb
z_amb_far   = z_exit_top + L_amb

# y+-Schaetzung (Blasius-Cf, glattes Rohr Do)
Re_orifice = U_ref * (Do * 1e-3) / nu_ref
Cf         = 0.046 * Re_orifice ** -0.2
u_tau      = U_ref * math.sqrt(Cf / 2.0)
yplus_wall = (h_wall * 1e-3) * u_tau / nu_ref


def feasible():
    """True, wenn die radiale Schalen-Abfolge konsistent ist."""
    return (a > 0.0) and (R0_corner > 0.0) and (R1 > R0_corner) and (R_wall > R1)


def report():
    print("== OFFIZIELL (Zhang Fig.1b, mm/deg) ==")
    print(f"  Di={Di} Ds={Ds} Ls={Ls} alpha={alpha}  Do={Do} Lo={Lo} Lk={Lk} n={n_inlet}")
    print(f"  theta (abgeleitet aus R_exit={R_exit}) = {theta:.1f} deg")
    print("== QUERSCHNITT (abgeleitet) ==")
    print(f"  R_wall = Do/2 = {R_wall:.4f}")
    print(f"  BL-Band  t2 = {t2:.4f}  (>= film {film_thickness}? {t2 >= film_thickness})"
          f"  n_r2 = {n_r2} (uniform)")
    print(f"     h_wall = {h_wall:.4f} mm  ->  y+ ~ {yplus_wall:.0f}  (Re={Re_orifice:.0f})")
    print(f"  Uebergang t1 = {t1:.4f}  n_r1 = {n_r1}")
    print(f"  Kreis1 R1 = {R1:.4f}   Quadrat-Ecke R0 = {R0_corner:.4f}   a = {a:.4f}")
    print(f"  zentr. Quadrat c = {c:.4f}  (t3 = {t3:.4f})")
    print(f"  N_core = {N_core}   h_base = {h_base:.4f} mm  (N_Do = {N_Do})")
    print("== VERHAELTNISSE (zu R_wall = Do/2) ==")
    print(f"  a/R_wall = {a/R_wall:.3f}   t1/R_wall = {t1/R_wall:.3f}"
          f"   t2/R_wall = {t2/R_wall:.3f}   t3/R_wall = {t3/R_wall:.3f}")
    print(f"  feasible = {feasible()}")
    print("== z-AUFBAU (mm; z=0 = Orifice-Eintritt) ==")
    print(f"  Drallkammer Rs={Rs:.2f}  z[{z_swirl_bot:.2f} .. {z_contr_bot:.2f}]  (Ls={Ls})")
    print(f"  Kontraktion Rs->R_wall  z[{z_contr_bot:.2f} .. 0]  Lc={Lc:.2f} (alpha={alpha})")
    print(f"  Orifice     R_wall      z[0 .. {z_orif_top:.2f}]  (Lo={Lo})")
    print(f"  Expansion   R_wall->{R_exit:.1f}  z[{z_orif_top:.2f} .. {z_exp_top:.2f}]"
          f"  Lk={Lk} theta={theta:.1f}")
    print(f"  Exit        R_exit={R_exit:.1f}  z[{z_exp_top:.2f} .. {z_exit_top:.2f}]"
          f"  (L_exit={L_exit})")
    if AMBIENT:
        print(f"  Ambient     R_amb ={R_amb:.1f}  z[{z_exit_top:.2f} .. {z_amb_far:.2f}]"
              f"  (L_amb={L_amb}, mid@{z_amb_mid:.2f})")
    print(f"  Aussenrand R_outer=3c={R_outer:.2f} (konstant, nur Nozzle-Segmente)")


if __name__ == "__main__":
    report()
