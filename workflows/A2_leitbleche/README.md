# A2 Leitbleche -- Mother-Topologie (konzentrische Vanes)

Streamlit-Demo-Workflow für CFDPredict A2 in der **Mutter-Topologie**:
konzentrische Bogenwände um einen gemeinsamen Bend-Mittelpunkt, jedes
Blech mit eigenem Radius `vane_r_i` und optionaler Downstream-Verlängerung
`vane_ext_i` parallel zur Auslaufwand.

Schwester-Workflow: [`../A2_brand_topology/`](../A2_brand_topology/) --
Brand-Paper-konforme parallel-translatierte Topologie für
Lehrbuch-Reproduktionen (Brand Kap. 7.3). Beide Workflows produzieren
`postProcessing/`-Outputs + `summary.json` + `summary.csv` mit
identischen Spalten, so dass sich Sweeps aus beiden in `pandas` mergen
lassen.

## Wann diese Topologie wählen?

Siehe `applications/A2_leitbleche/BAYES_TOPOLOGY_CHOICE.md` im CFDPredict-
Repo für die ausführliche Argumentation. Kurz:

| Anwendungsfall | Topologie |
|---|---|
| **Bayes-Optimierung / Aktives Lernen** (AP2/AP3) | **Mother** -- fixed-dim GP-Input via `vane_r_i` + `vane_ext_i` |
| **ROM-Trainingsdaten erzeugen** | **Mother** -- gleiche Topologie ⇒ konsistente Snapshots |
| **Brand-Paper-Reproduktion** (Showcase, Lehrbuch) | brand_topology |
| **Parameterstudie bei fester Vane-Anzahl** | egal (Mother flexibler) |

## Parameter (Form-Felder, MVP mit N_max=3)

Drei Vane-Slots, in radialer Reihenfolge von aussen (`vane_r1`, an der
Aussenwand) nach innen (`vane_r3`, an der Innenwand). Validierung:

- Radius muss in `(R-H/2, R+H/2)` liegen (strikt zwischen den Wänden).
- Reihenfolge muss strikt absteigend sein: `vane_r1 > vane_r2 > vane_r3`.
- Downstream-Verlängerung `vane_ext_i` in `[0, L_out]`.
- `number of vanes` (0..3) wählt, wie viele der drei aktiv sind. 0 = leerer
  Bogen (Baseline für Vergleichsläufe).

> N_max>3? Aktuell in der Streamlit-Form fest. Für N_max=6 oder 8 (Bayes
> mit grösserem Parameterraum) `interface.json` direkt um `vane_r4/...` +
> `vane_ext4/...` erweitern. Der Generator (`gen_blockmesh.py`) hat keine
> N-Obergrenze.

## Auto-Layers (Re-adaptive Wandauflösung)

Der Default `auto layers = true` berechnet aus `(Re, y+_target, ny)` die
Anzahl `refineWallLayer`-Passes so, dass y+ am ersten Zellzentrum
Re-unabhängig ungefähr bei `y+_target` landet (Dresden-Empfehlung: Gitter
muss Re-adaptiv sein, sonst schwankt y+ über den Re-Bereich um Faktor
~60 und Diskretisierungsfehler werden Re-abhängig).

- `y+_target = 40` (Default): Wandfunktions-Regime, robust + billig.
  Geeignet für Bayes/ROM/Standard-Bewertung.
- `y+_target = 1`: wand-aufgelöst (5-10x teurer pro Case, braucht
  potentialFoam-Initialisierung), für Validierung gegen Paper.

`auto layers = false` deaktiviert das, dann zählt nur `wall-layer splits`.

## Solver-Pipeline

Aus dem `case_template/Allrun`:

1. `blockMesh` -- konzentrische Bogen-Bands inkl. Vane-Interfaces als
   internal Faces; cellZones `band{j}` pro Radialband.
2. `topoSet` -- per Blech `boxToFace`: Bogen + Downstream-Verlängerung
   `L_ext_i` (`x>=0` filtert Inlet raus, `y<=R+ext` begrenzt Auslauf).
3. `createBaffles` -- Vane-Interfaces zu `vane{i}_master/_slave`-Wänden.
4. `refineWallLayer` -- N-mal an `walls` + `vaneN_master/_slave`, jeder
   Pass halbiert die wandnächste Zelle.
5. `decomposePar` + parallel `potentialFoam -initialiseUBCs` (Init) +
   `mpirun simpleFoam` (Solver) + `reconstructPar -latestTime`.

functionObjects in `controlDict`:

- `inletPressure`, `outletPressure` -- `areaAverage(p)` an inlet/outlet
- `outletVorticity` -- `areaAverage(magVorticity)` (sekundärströmung)
- `outletUniformity` -- `uniformity(U)` (kundenrelevante Metrik gamma)
- `yPlus` -- Diagnose pro Wand-Patch

## Outputs

`results.zip` enthält:

```
progress/
  1_info.txt             # live-Log (Setup + Run-Status)
  summary.json           # sweep_meta + per-case stats (machine-readable)
  summary.csv            # flach, sweep_meta pro Zeile dupliziert (pandas)
  <case>_pv_meshView.png # Gitterstruktur (Top-Down, Surface With Edges)
  <case>_pv_streamlines.png  # mid-slice Streamlines (diagonal seed)
  <case>_pv_3slice_iso.png   # 3 orthogonale Slices mit Vane-Tubes
<case>/                   # ParaView-öffenbar via paraFoam / *.foam
  system/, constant/, 0/, <latestTime>/, postProcessing/, logs/, *.foam
```

## Bayes-Loop-Verwendung (Python-API)

```python
from tools.run_workflow import run_workflow_headless
import json

def evaluate(theta):
    """theta = (vane_r1, vane_r2, vane_r3, vane_ext1, vane_ext2, vane_ext3)"""
    settings = json.loads((HERE / "examples/default.json").read_text())
    for i, r in enumerate(theta[:3], 1):
        settings[f"vane_r{i} [m]"] = r
    for i, e in enumerate(theta[3:], 1):
        settings[f"vane_ext{i} [m]"] = e
    result = run_workflow_headless(
        "A2_leitbleche", settings,
        workdir=Path(f"./runs/iter_{iter_id:04d}"))
    if result["exit_code"] != 0:
        return PENALTY                # solver crash / out-of-range / kill
    case = result["summary"]["cases"][0]
    return -case["uniformity"], case["dP"]   # maximize gamma, minimize dP
```

Siehe auch das geplante Dresden-Bundle (Phase 1) für ein vollständiges
ParEGO-Skelett.

## Smoke-Test (ohne OpenFOAM)

```bash
cd workflows/A2_leitbleche
python test_build_offline.py
```

Mit dem Repo-CLI:

```bash
python tools/run_workflow.py A2_leitbleche <settings.json> --dry-build
```

## Bekannte Einschränkungen / Gotchas

- **`R > H/2` ist Pflicht — wird seit 2026-07 erzwungen** (Worker-Guard +
  Defense-in-depth in `gen_blockmesh._radii`). Bei `R ≤ H/2` bricht der Job
  jetzt vor dem Build mit klarer Meldung ab (`9_error.txt`), statt dass
  `blockMesh` mit `inward-pointing faces` stirbt. Physikalisch validiertes
  Regime bleibt R/H ≳ 0.6; für H=4 z. B. R=2.4 (R/H=0.6), Vanes 3.4/2.4/1.4.
  Historie: Prod-Crash 2026-06-09 (`20260609_095558_A2-Brand-Leitbleche`,
  H=4 mit Default R=1.5 → R/H=0.375).
- **Mesh-/Solver-Fehler sind seit 2026-07 fail-loud.** `run_case` wertet den
  Allrun-Exit-Code aus und prüft `constant/polyMesh` + mindestens einen
  `Time =`-Schritt in `logs/solver.log`; bei Verstoß endet der Job als FAILED
  (`9_error.txt` + Log-Tail des gescheiterten Schritts, nonzero exit) statt
  „Finished!" mit `null`-summary und leerer results.zip. Gleiches Muster im
  Schwester-Workflow `A2_brand_topology` (dort mark-and-continue pro
  Sweep-Case, Gesamtjob endet FAILED, gesunde Cases bleiben in results.zip).

## Referenzen

- `applications/A2_leitbleche/README.md` -- Iterationen 1-8 der A2-Entwicklung
- `applications/A2_leitbleche/VALIDATION.md` -- Brand-2020-Trend-Validierung
- `applications/A2_leitbleche/BAYES_TOPOLOGY_CHOICE.md` -- Topologie-Empfehlung
- `SYSTEM_REQUIREMENTS.md` (Repo-root) -- OF / Python / ParaView
