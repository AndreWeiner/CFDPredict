# A2_brand_topology -- Brand 2020 Kap. 7 Reproduktion

Streamlit-Workflow für CFDPredict A2 (Leitbleche in einem 90°-Kanal-Bogen).
Reproduziert die Brand-2020-Serien aus Kapitel 7.3.

## Geometrie-Konvention (Brand 2020 S. 114)

- Windkanal-Seitenlänge: **4 m** (`H_phys = 4 m`)
- Innen- und Außen-Verrundungs-Radius beide: **1 m** (`r_bend = R - H/2 = 0.25·H`)
- Default `R/H = 0.75` (sodass `R - H/2 = 0.25 = r_bend`)
- Leitbleche: mittlerer Radius 1 m, Länge Viertelkreis, Dicke 90 mm
  *(in diesem Generator als dimensionslose Baffles realisiert, Dicke nicht modelliert)*

Brand-Pos-Konvention: Pos N = `(N + 1) · 100 mm` Abstand zur Innenwand.
Normalisiert auf `H = 1`: `fraction = (N + 1) / 40`.

| Brand-Pos | Abstand (mm) | fraction (H=1) |
|-----------|--------------|----------------|
| 0         |   100        | 0.025          |
| 1         |   200        | 0.050          |
| 5         |   600        | 0.150          |
| 11        | 1 200        | 0.300          |
| 19        | 2 000        | 0.500          |
| 29        | 3 000        | 0.750          |

## Drei Brand-Serien als Workflow-Inputs

Die Brand-Aufsätze 7.3.2 / 7.3.3 / 7.3.4 lassen sich über das `number of vanes`-Feld
abbilden:

### Serie A: Single-Vane-Position-Sweep (Brand Fig. 7.14 / 7.15)

Mehrere Cases mit jeweils **genau 1 Vane** an verschiedenen Pos.

```
number of vanes  =  P[1, 3, 5, 11, 19, 29]
```

Gibt 6 Cases, jeder mit 1 Vane an Brand-Pos 1/3/5/11/19/29. Chart-x-Achse
schaltet automatisch auf „Brand-Pos N".

### Serie B: Multi-Vane-Anzahl-Sweep (Brand Fig. 7.17)

Mehrere Cases mit unterschiedlicher **Anzahl equidistanter Vanes**.

```
number of vanes  =  [0, 3, 5, 7, 11]
```

Gibt 5 Cases mit 0/3/5/7/11 equidistanten Vanes (Fraktionen i/(N+1)).

### Serie C: Multi-Vane mit expliziten Brand-Positionen (Brand Fig. 7.18)

**Ein** Case mit allen 7 Brand-Bleche gleichzeitig (z. B. Pos 1/5/11/19/29
plus weitere).

```
number of vanes  =  P[1, 5, 11, 19, 25, 29, 33]
```

…oder mit expliziten Pathline-Bruchteilen:
```
number of vanes  =  [0.025, 0.125, 0.300, 0.500, 0.750]
```

### Modifizierte Leitbleche (Brand Fig. 7.20 / 7.21)

Brand verlängert das Vane-Ende um 500 mm downstream. Im Workflow:

```
trail length [H]  =  0.125     # = 500 mm / 4000 mm
number of vanes   =  P[1, 5, 11, 19, 29]
```

Bei `trail length > 0` werden die Vanes als gerade Strecken nach dem
Bogen-Ende fortgesetzt (parallel zur Auslaufwand).

## Eingabemodi (Auto-Detection per Prefix / Typ)

| Eingabe                              | Modus                                         |
|--------------------------------------|-----------------------------------------------|
| `"3"`                                | 1 Case mit 3 equidistanten Vanes              |
| `"[0, 3, 5]"`                        | 3 Cases, equidistant (Anzahl-Sweep)           |
| `"[0.025, 0.125, 0.3]"`              | 1 Case mit 3 expliziten Pathline-Bruchteilen  |
| `"P11"`                              | 1 Case, 1 Vane bei Brand-Pos 11               |
| `"P[1, 5, 11]"`                      | 3 Cases, je 1 Vane (Brand-Pos-Sweep)          |

Cartesian-Product mit `cells across inlet N_inlet` (Liste) für
Gitterkonvergenz-Studien.

## Offline-Tests

```bash
cd workflows/A2_brand_topology
python test_build_offline.py
```

8 Smoke-Tests inkl. Mode-1/2/3-Parser, `brand_pos_to_fraction`,
Pos-Single-Vane-Sweep, und „modifizierte Leitbleche" mit `Trail_Lane`-Blöcken.
OpenFOAM ist **nicht** erforderlich (`do_run=False`).

## Solver-Lauf

Auf der Streamlit-VM (`a demo VM`) oder a production HPC server mit OpenFOAM v2512.
`STREAMLIT_OPENFOAM_BASHRC` Env-Var honoriert.

## Verwandte Dokumente

- `CFDPredict:applications/A2_leitbleche/brand_topology/README.md` -- Generator-Doku
- `CFDPredict:applications/A2_leitbleche/brand_topology/handoff.md` -- Hintergrund + Bugs
- Brand-Showcase reproduction (internal DHCAE reference)
