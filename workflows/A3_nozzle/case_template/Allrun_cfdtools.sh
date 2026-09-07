#!/bin/bash
# A3_nozzle -- Allrun-Skript fuer manuelle Builds, DHCAE-Default 128 Kerne.
# Wird von py_A3_nozzle.py NICHT verwendet (das macht alles ueber den Worker,
# empfohlen: tools/run_workflow.py). Nuetzlich fuer manuelle Sweeps/Tests ohne
# Streamlit-UI. Portabel: setze STREAMLIT_OPENFOAM_BASHRC auf euren
# etc/bashrc-Pfad, sonst greift der cfdtools-Default (v2506/v2512 unter /opt).
#
# Pipeline:
#   1. blockMesh + surfaceFeatureExtract + snappyHexMesh -overwrite
#   2. createPatch (case-template Disk falls vorhanden) -overwrite
#   3. checkMesh
#   4. BL-Split (Lesson 038-v2, nur wenn topoSetDict_bl_split + createPatchDict_bl_split
#      im system/ liegen, OPT-OUT durch Loeschen): topoSet + createPatch + 2x
#      refineWallLayer + 0/-Patch fuer wall_top/wall_straight/wall_expansion
#   5. Zhang-konforme BCs (Lesson 047, KRITISCH): patcht 0/p_rgh outlet zu
#      fixedValue 0 + 0/alpha.water outlet zu inletOutlet 0. ueberspringt mit
#      OPT-OUT-Flag SKIP_ZHANG_BCS=1.
#   6. setFields (wenn dict da)
#   7. decomposePar -force + interFoam parallel
#
# Use OPT-OUTs by environment vars:
#   SKIP_BL_SPLIT=1     skip BL refinement
#   SKIP_ZHANG_BCS=1    keep DHCAE advective outlet BCs (NICHT Zhang-konform)
#   N_PROCS=N           override decomposePar processor count (default 128)

set -eo pipefail   # NOT -u: OF bashrc references WM_PROJECT_DIR before setting it

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$CASE_DIR"
mkdir -p logs

if [[ -n "${STREAMLIT_OPENFOAM_BASHRC:-}" ]]; then
    source "$STREAMLIT_OPENFOAM_BASHRC"
else
    source /opt/OpenFOAM/OpenFOAM-v2506/etc/bashrc 2>/dev/null \
      || source /opt/OpenFOAM/OpenFOAM-v2512/etc/bashrc 2>/dev/null \
      || true
fi

N_PROCS="${N_PROCS:-128}"

echo "==[1/7] blockMesh =="
blockMesh > logs/blockMesh.log 2>&1
echo "==[2/7] surfaceFeatureExtract =="
surfaceFeatureExtract > logs/surfaceFeatureExtract.log 2>&1
echo "==[3/7] snappyHexMesh =="
snappyHexMesh -overwrite > logs/snappy.log 2>&1
if [[ -f system/createPatchDict ]]; then
    echo "==[3b]   createPatch (case-template) =="
    createPatch -overwrite > logs/createPatch.log 2>&1
fi
echo "==[4/7] checkMesh (pre-BL) =="
checkMesh -allGeometry -allTopology > logs/checkMesh_preBL.log 2>&1
tail -5 logs/checkMesh_preBL.log

if [[ -z "${SKIP_BL_SPLIT:-}" && -f system/topoSetDict_bl_split && -f system/createPatchDict_bl_split ]]; then
    echo "==[5/7] BL-Split (Lesson 038-v2) =="
    topoSet -dict system/topoSetDict_bl_split > logs/topoSet.log 2>&1
    createPatch -dict system/createPatchDict_bl_split -overwrite > logs/createPatch_bl.log 2>&1
    refineWallLayer -overwrite '(wall_straight wall_expansion)' 0.5 > logs/refineWallLayer_1.log 2>&1
    refineWallLayer -overwrite '(wall_straight wall_expansion)' 0.5 > logs/refineWallLayer_2.log 2>&1
    # 0/* Felder ergaenzen um neue Patches (kloning wall-BC -> wall_top/_straight/_expansion)
    python3 << 'PYEOF'
import re
from pathlib import Path
for fn in ("U", "k", "epsilon", "nut", "p_rgh", "alpha.water"):
    fp = Path("0") / fn
    if not fp.exists():
        continue
    s = fp.read_text(encoding="utf-8")
    if "wall_top" in s:
        continue  # already patched
    # single-line wall BC (e.g. "    wall { type noSlip; }")
    m = re.search(r"^(\s*)wall(\s*\{[^}]*\}\s*)$", s, flags=re.M)
    if m:
        indent, bc = m.group(1), m.group(2).rstrip("\n")
        block = (indent + "wall_top" + bc + "\n" +
                 indent + "wall_straight" + bc + "\n" +
                 indent + "wall_expansion" + bc)
        s2 = re.sub(r"^(\s*wall\s*\{[^}]*\}\s*)$",
                    r"\1\n" + block, s, count=1, flags=re.M)
        fp.write_text(s2, encoding="utf-8", newline="\n")
        continue
    # multi-line wall BC
    m = re.search(r"^(\s*)wall\s*\{([^}]*)\}", s, flags=re.M | re.S)
    if m:
        indent, body = m.group(1), m.group(2)
        block = (indent + "wall_top\n" + indent + "{" + body + "}\n" +
                 indent + "wall_straight\n" + indent + "{" + body + "}\n" +
                 indent + "wall_expansion\n" + indent + "{" + body + "}")
        s2 = re.sub(r"^(\s*wall\s*\{[^}]*\})",
                    r"\1\n" + block, s, count=1, flags=re.M | re.S)
        fp.write_text(s2, encoding="utf-8", newline="\n")
PYEOF
    checkMesh -allGeometry -allTopology > logs/checkMesh_postBL.log 2>&1
    tail -5 logs/checkMesh_postBL.log
else
    echo "==[5/7] BL-Split SKIPPED (SKIP_BL_SPLIT=1 oder Dicts fehlen) =="
fi

if [[ -z "${SKIP_ZHANG_BCS:-}" ]]; then
    echo "==[6/7] Zhang-konforme outlet-BCs (Lesson 047, KRITISCH!) =="
    python3 << 'PYEOF'
import re
from pathlib import Path
# p_rgh outlet -> fixedValue 0
p = Path("0/p_rgh")
if p.exists():
    s = p.read_text(encoding="utf-8")
    s2 = re.sub(r"^(\s*)outlet\s*\{[^}]*\}",
                r"\1outlet         { type fixedValue; value uniform 0; }",
                s, count=1, flags=re.M | re.S)
    p.write_text(s2, encoding="utf-8", newline="\n")
    print("  0/p_rgh outlet -> fixedValue 0")
# alpha.water outlet -> inletOutlet 0
p = Path("0/alpha.water")
if p.exists():
    s = p.read_text(encoding="utf-8")
    s2 = re.sub(r"^(\s*)outlet\s*\{[^}]*\}",
                r"\1outlet\n\1{\n\1    type            inletOutlet;\n"
                r"\1    inletValue      uniform 0;\n\1    value           uniform 0;\n\1}",
                s, count=1, flags=re.M | re.S)
    p.write_text(s2, encoding="utf-8", newline="\n")
    print("  0/alpha.water outlet -> inletOutlet 0")
PYEOF
else
    echo "==[6/7] Zhang-BCs SKIPPED (SKIP_ZHANG_BCS=1) -- DHCAE advective outlet behalten =="
fi

if [[ -f system/setFieldsDict ]]; then
    echo "==     setFields =="
    setFields > logs/setFields.log 2>&1
fi
rm -f 0/cellLevel 0/pointLevel   # snappy leftover, mismatch nach refineWallLayer

echo "==[7/7] decomposePar -force + interFoam (np=$N_PROCS) =="
decomposePar -force > logs/decomposePar.log 2>&1
nohup mpirun --bind-to none -np $N_PROCS interFoam -parallel > logs/solver.log 2>&1 &
sleep 5
echo "  active mpirun in this case:"
for p in $(pgrep -f "mpirun.*interFoam"); do
    cwd=$(readlink /proc/$p/cwd 2>/dev/null || true)
    [[ "$cwd" == "$CASE_DIR" ]] && echo "    PID $p"
done
echo "==[done] tail logs/solver.log for live progress."
