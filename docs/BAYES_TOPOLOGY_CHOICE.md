# A2 — Topologie-Wahl für Bayes-Optimierung / Aktives Lernen

**Adressat:** Projektpartner TU Dresden (André Weiner)
**Stand:** 2026-05-26
**Autor:** Martin Becker (DHCAE)

## Kurzempfehlung

Für die Bayes-Optimierung und das aktive Lernen in AP2/AP3 empfehlen wir die
**konzentrische A2-Mutter-Topologie** (`applications/A2_leitbleche/`) gegenüber
der parallel-translatierten **A2-Brand-Topologie**
(`applications/A2_leitbleche/brand_topology/`).

Die Brand-Topologie ist und bleibt nützlich — für Brand-konforme
Reproduktionen (Kap. 7.3 Showcase im Web-Portal) und für rein visuelle Parameterstudien
mit fester Vane-Anzahl. Für die GP-/EI-Schleife ist sie strukturell
ungeeignet.

## Warum die Brand-Topologie für Bayes nicht passt

In der Brand-Topologie liegen die Leitbleche als parallel-translatierte
Bogenwände im Bogenbereich, alle mit identischem Krümmungsradius `r_bend`.
Die Auslegungs­variable des Generators ist die Liste
`inner_arc_fractions: list[float]` in (0,1) — sie kodiert die Lane-Trennlinien
entlang einer Pathline­koordinate.

**Konsequenz für die Optimierung:**

1. **Variable Länge.** Eine N-Vane-Konfiguration und eine (N±1)-Vane-Konfiguration
   leben in unterschiedlich-dimensionalen Räumen. Der GP-Kern (RBF, Matérn, …)
   ist auf einem Fixed-Dimensional Inputraum definiert; eine Liste variabler
   Länge ist kein gangbarer GP-Input ohne komplexe Tricks (Set-Kernels, RKHS
   für Maße, etc.).

2. **Topologie-Sprung pro Blech.** Jeder zusätzliche Bruchteil in der Liste
   fügt einen weiteren Block zur blockMesh-Topologie hinzu (Lane → zwei Lanes
   mit Trennlinie + neuem Vane-Patch). Es gibt keine kontinuierliche „weiche"
   Deaktivierung eines Blechs — beim Entfernen verschwindet der zugehörige
   Block samt seinen Nachbarschafts­beziehungen.

3. **Konstantes `r_bend`.** Alle Bleche teilen dasselbe `r_bend`; das ist
   die Brand-Treue, schränkt aber das Design­vokabular ein (keine
   Radius-Variation pro Blech).

Theoretisch kann man die Lücke umgehen — z. B. mit einem GP über die Anzahl
N (kategorisch) gekoppelt mit einem GP über die Lane-Bruchteile bei festem N,
oder mit kontextuellen Bandit-Ansätzen. Wir halten das für den
ZIM-Projektkontext (10–16 Designparameter, ROM-Trainingsdaten erzeugen)
für unnötige Komplexität.

## Warum die A2-Mutter passt

In der Mutter-Topologie liegen die Bleche als **konzentrische Bogen­wände um
einen gemeinsamen Bogen­mittelpunkt** — geometrisch das, was Brand 2020 Fig.
7.14 als Untersuchungs­variante zeichnet (verschiedene Krümmungsradien für
ein einzelnes Blech). Der Generator parametrisiert pro Blech:

```
vane_r_i     ∈ [R - H/2 + ε,  R + H/2 - ε]    # Radius des Blechs (kontinuierlich)
vane_ext_i   ∈ [0,            L_ext_max]      # downstream-Verlängerung (kontinuierlich)
```

Bei festem `N_max` (z. B. 6 oder 8) ist der Inputraum
`R^(2·N_max)` — ein klassisches GP-Setup. Ein Blech wird durch
`vane_ext_i = 0` und/oder durch Aus­schieben des Radius an die
Wand kontinuierlich „deaktiviert" — keine Topologie-Sprünge.

**Brand-Trend-Validierung** wurde auf der Mutter durchgeführt und liegt
dokumentiert vor (`applications/A2_leitbleche/VALIDATION.md`, Iter 8):
Vorzeichen-, Positions- und Anzahl-Trends aus Brand Kap. 7 sind in 2D+3D
reproduziert. Die Single-Vane-Positionsstudie (`_obsolet_test_brand_classic/`
in der Showcase-Sammlung; mit korrigiertem H_phys=4 m nachzuziehen) zeigt
eine **nicht-monotone ω-Antwortlandschaft** über die Vane-Position
(globales Minimum bei Pos 29, lokales Minimum bei Pos 5) — genau die
Art von nicht-trivialer Antwortfunktion, bei der GP/EI ihren Mehrwert
zeigen.

## Was das praktisch heißt

- **AP2/AP3 (DHCAE-Anteil, Bayes-Optimierung + ROM-Trainingsdaten):**
  Mutter-Topologie mit `N_max = 6..8` `vane_r{i}` + `vane_ext{i}`
  ⇒ 12–16 Bayes-Variablen. Im Parameterbudget (Dresden-Empfehlung: ≤ 10,
  max. 16 Designparameter pro Anwendung) gut darstellbar.

- **Brand-Showcase (DHCAE-Anteil, AP5 Cloud-Integration):**
  Brand-Topologie über das Streamlit-Portal — drei deployed Workflow-Konfigurationen
  reproduzieren Kap. 7.3.2 /
  7.3.3 / 7.3.4 (Single-Vane-Sweep / Anzahl-Sweep / modifizierte Vanes).
  Quantitativer Brand-Vergleich qualitativ (Brand-Original ist Discovery-
  Live, low-fidelity).

- **Schnittstelle:** Beide Topologien teilen `params.json`-Schema-Konvention,
  Patch-Naming (`inlet/outlet/walls/vane{i}_master|_slave`) und das
  `case_template/`. Postprocessing-Pipeline (functionObjects, surfaceFieldValue
  für ω/Uniformity/dP, summary.json+csv) ist identisch — Bayes- und Showcase-
  Daten sind direkt vergleichbar bzw. kombinierbar.

## Offene Punkte / zur Abstimmung mit TUD

1. **`N_max` für Bayes-Mutter** — 6 oder 8? Brand zeigt die Asymptote in
   Fig. 7.17 ab N ≈ 6, unsere eigene 7.3.3-Reproduktion (13 Cases
   N=0..12) bestätigt das mit Sättigung von dP_kin ab N=5 auf ~0.44.
   `N_max=8` gibt etwas Sicherheit, kostet aber 4 zusätzliche
   Bayes-Variablen.

2. **Sentinel für „Vane deaktiviert"** — sauberer ist
   `vane_ext_i = 0 ∧ vane_r_i an die Wand geschoben`, alternativ ein
   diskreter `active_i ∈ {0,1}` zusätzlich (gemischt-kontinuierlich GP via
   BoTorch). Wir tendieren zu Variante 1 (rein kontinuierlich, kein
   gemischtes Modell).

3. **Ziel­größe(n).** Aktuell sind drei FOs produktiv (dP_kin, ω_outlet,
   Uniformity). Für die ZIM-Kundenmetrik (Propeller­anströmung) ist die
   Uniformity γ die abrechenbare Größe. Für die Bayes-Schleife würden wir
   ParEGO mit Hauptziel `γ↑` + Nebenziel `dP_kin↓` vorschlagen — die ROM
   muss dann beides erlauben.

## Referenzen

- `applications/A2_leitbleche/README.md` — A2-Mutter-Topologie
- `applications/A2_leitbleche/brand_topology/README.md` — Brand-Topologie
- `applications/A2_leitbleche/brand_topology/handoff.md` — Streamlit-Integrationsstand
- `applications/A2_leitbleche/VALIDATION.md` — Brand-Trend-Validierung Iter 8
- Brand-Reproduktions-Showcase im Streamlit-Portal (interne DHCAE-Referenz)
