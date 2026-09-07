# A3 Pressure-Swirl Nozzle — Zhang 2023

Parametrische Druck-Drall-Düse (Pressure-Swirl Atomizer) nach
**Zhang et al. 2023, *Atomization and Sprays***. Tangentiale Wasser-
einläufe erzeugen einen Drall in der Drallkammer, der bei der Kontraktion
zur Bohrung beschleunigt und am Bohrungsaustritt einen rotierenden
Wasserfilm bildet. Im Ambient bildet sich daraus der typische
Hohlkegel-Spray.

## Geometrie (Zhang Fig. 1b)

Achssymmetrischer Halbschnitt von oben nach unten (Strömungsrichtung):

- **Drallkammer** (Durchmesser **Ds**, Länge **Ls**) mit
  **n_inlet** Tangentialeinläufen vom Durchmesser **Di**
- **Kontraktion** Ds → Do mit Halbwinkel **α**
- **Bohrung** Durchmesser **Do**, Länge **Lo**
- **Expansion** (Lampenschirm) Länge **Lk**, Auslassradius **R_exit**;
  Öffnungs-Halbwinkel **θ** = atan((R_exit − Do/2) / Lk).

## Zhang DOE (Taguchi L16)

Vier Faktoren mit je 4 Levels — der DOE-Code im Interface (z.B. *A4B1C2D4*)
überschreibt die Einzelwerte unten:

| Letter | Variable | Level 1 | Level 2 | Level 3 | Level 4 |
|--------|----------|---------|---------|---------|---------|
| A      | θ [°]    | 0       | 10      | 20      | 30      |
| B      | Do [mm]  | 4       | 5       | 6       | 7       |
| C      | α [°]    | 30      | 45      | 60      | 75      |
| D      | Ds [mm]  | 9       | 10      | 11      | 12      |

Zhang's Baseline-Konfiguration (Section 2.2, Fig. 3) ist **A1B2C2D2**
(θ=0, Do=5, α=45, Ds=10). Sein Optimum (Table 6) ist **A4B1C2D4**.

## Numerik (Default-Setup, Lessons 029-048)

- **Solver**: interFoam (VOF, MULES), endTime ~25 ms
- **Turbulenz**: RNG k-ε + multiphaseStabilizedTurbulence
- **Schemes**: linearUpwind grad(U), upwind k/ε, vanLeer alpha — der
  div(rhoPhi,U)-Scheme ist DER dominante Hebel auf p_inlet
- **PIMPLE**: 1 outer / 2 inner / 1 nonOrth — verifiziert stabil
- **MULES**: nAlphaCorr 1 / nAlphaSubCycles 3
- **BC**: Zhang-paper-konform (p_rgh outlet fixedValue 0, alpha
  inletOutlet 0); Toggle "Zhang outlet BCs" = 1
- **Inlet-Turbulenz**: Intensity 1% (statt 5%) — Lesson 041

## Validierung (1:1 Zhang Baseline A1B2C2D2, Lk=5, 12 Cases)

- **Film** (Sample-Linie z=−24.9 mm): 047 RNG + Zhang outlet → **−1%
  vs Zhang 0.655 mm** ⭐
- **U_θ Peak** (z=−24.9): 048 kOmegaSST → **−11% vs Zhang ±30 m/s** ⭐
- **p_inlet**: strukturelle **+24-34%** Abweichung von Zhang's 1.0 MPa,
  robust gegen Wall-Function / Inlet-Turb / k-ε-Discretization
  (vermutlich OpenFOAM-vs-Fluent RNG-Implementation). 049
  (compressibleInterFoam) testet ob Air-Core-Kompression das schließt.

## Output

- `case_<name>/overview.png` + `overview.html` — Konfigurations-Snapshot
- `progress/A3_nozzle_convergence.png` — Live-Konvergenz
- `progress/<step>_info.txt` — Live-Log pro Pipeline-Schritt
- `results/` — finale Tabelle (Film, U_θ, p_inlet, Cone-Angle)

## Quellen

- Zhang, R. et al. (2023). *Effect of Structural Parameters on the
  Atomization Performance of a Pressure-Swirl Nozzle*. Atomization and
  Sprays.
- Lefebvre & McDonell (2017). *Atomization and Sprays*, 2nd Ed.
