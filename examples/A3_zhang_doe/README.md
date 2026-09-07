# Zhang 2023 L16-DOE Case-Konfigurationen

Vorgenerierte `interface.json`-Dateien für die 16 Cases des Taguchi-L16(4^4)-DOE
aus Zhang et al. 2023 (Sektion 4.3, Tabelle 3) — plus eine 17. Verifikation für
das vom Paper vorhergesagte Optimum `A4B1C2D4` (SMD 40.52 µm).

## Faktoren / Levels

| Code | Bedeutung              | Level 1 | Level 2 | Level 3 | Level 4 |
|------|------------------------|---------|---------|---------|---------|
| A    | Expansionswinkel β [°] | 0       | 10      | 20      | 30      |
| B    | Orifice-Durchmesser Do [mm] | 4  | 5       | 6       | 7       |
| C    | Kontraktionswinkel α [°] | 30    | 45      | 60      | 75      |
| D    | Drallkammer-Ø Ds [mm]   | 9       | 10      | 11      | 12      |

A=β treibt `R_exit` (= Do/2 + Lk·tan(β/2)). Die Worker-Funktion `parse_doe_code()`
decodet `A_iB_jC_kD_l` → überschreibt `Do`, `alpha`, `Ds`, `R_exit` in der
geladenen Schema-Datei.

## Inhalt

16 JSONs nach dem Standard-L16(4^4)-Taguchi-Array (`A1B1C1D1.json` ...
`A4B4C1D3.json`) plus die Verifikation `A4B1C2D4.json`. Alle haben
`enable ambient = 1`, `end time = 0.150 s`, `nProc = 128` — die schweren
Solver-Settings für Sprühkegel + Vordruck-Konvergenz.

`generate.py` regeneriert sie aus `../interface.json` (passt sich bei zukünftigen
Schema-Aenderungen automatisch an).

## Verwendung

### Einzeln (manuell via Portal)

```text
1. Streamlit-Portal öffnen → Workflow A3_nozzle
2. interface.json hochladen, z.B. A4B1C2D4.json
3. Submit
```

### Batch (headless, via tools/run_workflow.py)

Sequenziell alle 17:
```bash
cd <repo-root>
for cfg in workflows/A3_nozzle/zhang_doe_configs/*.json; do
  name=$(basename "$cfg" .json)
  python tools/run_workflow.py A3_nozzle \
    --interface "$cfg" \
    --series "runs/zhang_doe_$name"
done
```

Oder einzeln:
```bash
python tools/run_workflow.py A3_nozzle \
  --interface workflows/A3_nozzle/zhang_doe_configs/A4B1C2D4.json \
  --series runs/zhang_optimum
```

### Anmerkungen

- **Laufzeit pro Case** auf unserem produktiven HPC-Server (128 Kerne): ~3 h für `end time = 0.150 s`
  (Mesh ~460k Zellen, VOF + RNG-k-ε). 17 Cases sequenziell = ~50 h.
- **Bayes-Loop:** für eine Optimierung über die Levels (nicht-Taguchi) ist
  `tools/run_workflow.py:run_workflow_headless()` die Python-API; die
  Level-Tabellen in `py_A3_nozzle._DOE_A/B/C/D` machen die Discretisierung.
- **A4B1C2D4 ist NICHT im Standard-L16-Array** — es ist die vom Taguchi-Auswertung
  vorhergesagte Optimum-Kombination und wurde von Zhang als separater
  Verifikationslauf gerechnet. SMD 40.52 µm dort = unser Referenzziel für die
  PEGASUS-VOF-Lagrange-Validierung (Task 4, separates Setup mit größerer Domain).
