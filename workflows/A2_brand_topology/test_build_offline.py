#!/usr/bin/env python3
"""Offline smoke for A2_brand_topology: build the workflow's 6-tuple
schema into a series-dir (as the streamlit app does with
inject_schema:true), run the worker with do_source_of=False, do_run=False,
and assert that the case dirs are built.

OpenFOAM is NOT required. Run from inside the workflow dir:
    cd workflows/A2_brand_topology && python test_build_offline.py
"""
from __future__ import annotations

import copy
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import py_A2_brand_topology as worker

SCHEMA = json.loads((HERE / "interface.json").read_text(encoding="utf-8"))


def make_series(tmp: Path, name: str, **overrides) -> Path:
    """Materialize a copy of the workflow's 6-tuple schema, patching the
    default-value slot (index 1) for the given field names."""
    series = tmp / name
    series.mkdir()
    schema = copy.deepcopy(SCHEMA)
    for field, value in overrides.items():
        assert field in schema, f"unknown field {field!r}"
        schema[field][1] = value
    (series / "interface.json").write_text(
        json.dumps(schema), encoding="utf-8")
    return series


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)

        # 1) empty bend, single mesh -> 1 case, no createBaffles
        series = make_series(
            tmp, "empty_bend",
            **{"number of vanes": "0", "cells across inlet N_inlet": "8",
               "cells in z (nz)": "4"})
        worker.run_worker(str(series), do_source_of=False, do_run=False)
        case_dirs = sorted([p for p in series.iterdir()
                            if p.is_dir() and p.name.startswith("0_N")])
        assert len(case_dirs) == 1, \
            f"expected 1 case dir, got {[c.name for c in case_dirs]}"
        c0 = case_dirs[0]
        assert (c0 / "system" / "blockMeshDict").is_file()
        assert not (c0 / "system" / "createBafflesDict").is_file(), \
            "empty bend should not produce createBafflesDict"
        assert (series / "command_finished").is_file()
        print(f"[ok] empty-bend: {c0.name}")

        # 2) 2-vane case + N_inlet series "[8, 12]" -> 1 x 2 = 2 cases
        series2 = make_series(
            tmp, "two_vanes_mesh_sweep",
            **{"number of vanes": "2",
               "cells across inlet N_inlet": "[8, 12]",
               "cells in z (nz)": "4"})
        worker.run_worker(str(series2), do_source_of=False, do_run=False)
        case_dirs2 = sorted([p for p in series2.iterdir()
                             if p.is_dir() and p.name.startswith(("0_N", "1_N"))])
        assert len(case_dirs2) == 2, \
            f"expected 2 cases, got {[c.name for c in case_dirs2]}"
        for c in case_dirs2:
            assert (c / "system" / "blockMeshDict").is_file()
            assert (c / "system" / "createBafflesDict").is_file(), \
                f"{c.name} (N=2) should have createBafflesDict"
        print(f"[ok] N=2, N_inlet=[8,12]: {[c.name for c in case_dirs2]}")

        # 3) vane-sweep "[0, 3]" + single mesh -> 2 cases
        series3 = make_series(
            tmp, "vane_sweep",
            **{"number of vanes": "[0, 3]",
               "cells across inlet N_inlet": "8",
               "cells in z (nz)": "4"})
        worker.run_worker(str(series3), do_source_of=False, do_run=False)
        case_dirs3 = sorted([p for p in series3.iterdir()
                             if p.is_dir() and p.name.startswith(("0_N", "1_N"))])
        assert len(case_dirs3) == 2, \
            f"expected 2 cases, got {[c.name for c in case_dirs3]}"
        zero = next(c for c in case_dirs3 if "_N0_" in c.name)
        three = next(c for c in case_dirs3 if "_N3_" in c.name)
        assert not (zero / "system" / "createBafflesDict").is_file(), \
            "N=0 case must not have createBafflesDict"
        assert (three / "system" / "createBafflesDict").is_file(), \
            "N=3 case must have createBafflesDict"
        print(f"[ok] vane-sweep [0, 3]: empty + 3-vane")

        # 4) explicit non-equidistant fractions (Mode 2): Brand Fig 7.14
        brand = [0.02703, 0.08108, 0.16216, 0.32432, 0.59459]
        series4 = make_series(
            tmp, "brand_fig714",
            **{"number of vanes": str(brand),
               "cells across inlet N_inlet": "20",
               "cells in z (nz)": "4",
               "bend radius ratio R/H": "2.0"})
        worker.run_worker(str(series4), do_source_of=False, do_run=False)
        case_dirs4 = [p for p in series4.iterdir()
                      if p.is_dir() and p.name.startswith("0_N")]
        assert len(case_dirs4) == 1, \
            f"expected exactly 1 Mode-2 case, got {[c.name for c in case_dirs4]}"
        c4 = case_dirs4[0]
        # case dir name marks custom fractions with the 'x' suffix
        assert "x_" in c4.name, \
            f"Mode-2 dir should carry the 'x' tag, got {c4.name}"
        assert (c4 / "system" / "createBafflesDict").is_file()
        bmd = (c4 / "system" / "blockMeshDict").read_text(encoding="utf-8")
        # all five Brand pathline fractions must end up resolved -- spot-
        # check by counting baffle patches in createBafflesDict
        baf = (c4 / "system" / "createBafflesDict").read_text(encoding="utf-8")
        n_baffle_groups = baf.count("type            faceZone;")
        assert n_baffle_groups == len(brand), \
            (f"expected {len(brand)} baffle groups, got "
             f"{n_baffle_groups} in {c4.name}/createBafflesDict")
        print(f"[ok] Brand Fig 7.14 explicit-fractions: {c4.name}, "
              f"{n_baffle_groups} baffle groups")

        # 4b) verify results.zip exists, has progress + system + constant +
        # initial fields, and is bounded. (Mode 2 case from step 4 is
        # offline-built only -> no logs/, no postProcessing/, no
        # reconstructed solver time; 0/ + constant/ come from
        # case_template.)
        import zipfile
        rzip = series4 / "results" / "results.zip"
        assert rzip.is_file(), f"results.zip missing for {series4.name}"
        sz_mb = rzip.stat().st_size / 1024 / 1024
        assert sz_mb < 50, f"results.zip is {sz_mb:.2f} MB, expected < 50"
        with zipfile.ZipFile(rzip) as zf:
            names = zf.namelist()
        assert any(n.startswith("progress/") and n.endswith(".txt")
                   for n in names), \
            f"progress/*.txt missing in results.zip: {names[:5]}"
        assert any("/system/blockMeshDict" in n for n in names), \
            f"system/blockMeshDict missing in results.zip: {names[:5]}"
        assert any("/constant/" in n for n in names), \
            f"constant/ missing in results.zip (case should be openable " \
            f"in ParaView): {names[:5]}"
        assert any("/0/U" in n or "/0/p" in n for n in names), \
            f"0/ initial fields missing in results.zip: {names[:5]}"
        # processor*/ must still be excluded (those are the parallel split
        # mesh + fields, redundant with the reconstructed time)
        assert not any("/processor" in n for n in names), \
            f"processor*/ unexpectedly inside results.zip: " \
            f"{[n for n in names if '/processor' in n][:5]}"
        print(f"[ok] results.zip {sz_mb*1024:.1f} kB, "
              f"{len(names)} entries, ParaView-openable")

        # 5) parser unit-checks (no build): backwards compat + Mode 2 + edge
        from py_A2_brand_topology import parse_vane_specs
        assert parse_vane_specs("3") == [3]
        assert parse_vane_specs("[0, 3, 5]") == [0, 3, 5]
        assert parse_vane_specs("[3.0, 5.0]") == [3, 5]   # int-valued floats
        assert parse_vane_specs("0") == [0]
        assert parse_vane_specs("") == [0]
        s = parse_vane_specs("[0.5, 0.7]")
        assert len(s) == 1 and isinstance(s[0], list) and s[0] == [0.5, 0.7], s
        s = parse_vane_specs("[0.7, 0.5, 0.7]")   # dedupe + sort
        assert s == [[0.5, 0.7]], s
        s = parse_vane_specs("[1.5, 0.5]")   # out-of-range filtered
        assert s == [[0.5]], s
        print("[ok] parse_vane_specs covers Mode 1 + Mode 2 + edges")

        # 6) Mode 3 (Brand-Pos) parser
        from py_A2_brand_topology import (brand_pos_to_fraction,
                                          resolve_vane_spec)
        assert parse_vane_specs("P11") == [("P", 11)]
        assert parse_vane_specs("p11") == [("P", 11)]   # case-insensitive prefix
        assert parse_vane_specs("P[1, 5, 11]") == [
            ("P", 1), ("P", 5), ("P", 11)]
        assert parse_vane_specs("P[11, 5, 1]") == [
            ("P", 11), ("P", 5), ("P", 1)]    # order preserved
        # Brand-Pos -> fraction mapping (H_phys = 4 m, 100 mm step)
        assert brand_pos_to_fraction(0) == 0.025          # Pos 0 = 100 mm
        assert brand_pos_to_fraction(1) == 0.05           # Pos 1 = 200 mm
        assert brand_pos_to_fraction(11) == 0.3           # Pos 11 = 1200 mm
        assert brand_pos_to_fraction(29) == 0.75          # Pos 29 = 3000 mm
        # resolve_vane_spec for Mode 3
        assert resolve_vane_spec(("P", 11)) == [0.3]
        assert resolve_vane_spec(("P", [1, 5, 11])) == [0.05, 0.15, 0.3]
        print("[ok] Mode 3 (Brand-Pos) parser + brand_pos_to_fraction")

        # 7) Brand single-vane Pos-sweep "P[1, 11, 29]" -> 3 cases each
        # with EXACTLY ONE vane (Brand Fig. 7.14 / 7.15 reproduction)
        series_brand_sweep = make_series(
            tmp, "brand_pos_sweep",
            **{"number of vanes": "P[1, 11, 29]",
               "cells across inlet N_inlet": "8",
               "cells in z (nz)": "4",
               "bend radius ratio R/H": "0.75"})
        worker.run_worker(str(series_brand_sweep),
                          do_source_of=False, do_run=False)
        case_dirs_brand = sorted([p for p in series_brand_sweep.iterdir()
                                  if p.is_dir() and p.name.startswith(
                                      ("0_P", "1_P", "2_P"))])
        assert len(case_dirs_brand) == 3, \
            (f"expected 3 Brand-Pos cases, got "
             f"{[c.name for c in case_dirs_brand]}")
        for c in case_dirs_brand:
            # each case should have createBafflesDict for its 1 vane
            assert (c / "system" / "createBafflesDict").is_file(), \
                f"{c.name} (1 vane) should have createBafflesDict"
            baf = (c / "system" / "createBafflesDict").read_text(
                encoding="utf-8")
            n_baffle_groups = baf.count("type            faceZone;")
            assert n_baffle_groups == 1, \
                (f"{c.name} (Brand-Pos single vane) expected exactly 1 "
                 f"baffle group, got {n_baffle_groups}")
        # case dir names encode the Brand-Pos: 0_P1_N1_..., 1_P11_N1_..., etc.
        pos_in_names = sorted([c.name.split("_")[1][1:]
                               for c in case_dirs_brand])
        assert pos_in_names == ["1", "11", "29"], pos_in_names
        print(f"[ok] Brand-Pos single-vane sweep P[1, 11, 29]: "
              f"3 cases, {[c.name for c in case_dirs_brand]}")

        # 8) Brand 'modifizierte Leitbleche' = vane extension via trail length
        # Fig. 7.20/7.21: 500 mm extension at H_phys = 4 m -> trail = 0.125 H
        series_mod = make_series(
            tmp, "brand_modifizierte_leitbleche",
            **{"number of vanes": "P[1, 11, 29]",
               "trail length [H]": "0.125",
               "cells across inlet N_inlet": "8",
               "cells in z (nz)": "4",
               "bend radius ratio R/H": "0.75"})
        worker.run_worker(str(series_mod), do_source_of=False, do_run=False)
        case_dirs_mod = sorted([p for p in series_mod.iterdir()
                                if p.is_dir() and p.name.startswith(
                                    ("0_P", "1_P", "2_P"))])
        assert len(case_dirs_mod) == 3, \
            (f"expected 3 modified-vane cases, got "
             f"{[c.name for c in case_dirs_mod]}")
        # spot-check that blockMeshDict for the Pos11 case has the trail
        # block (Trail_Lane{n}) for the 2 lanes around the vane (n=1, n=2)
        c_pos11 = next(c for c in case_dirs_mod if "_P11_" in c.name)
        bmd = (c_pos11 / "system" / "blockMeshDict").read_text(
            encoding="utf-8")
        # the brand_topology generator names trail blocks Trail_Lane{n};
        # with 1 vane there are 2 lanes -> Trail_Lane1 + Trail_Lane2 expected
        assert "Trail_Lane1" in bmd, \
            f"Trail_Lane1 missing in {c_pos11.name}/system/blockMeshDict"
        assert "Trail_Lane2" in bmd, \
            f"Trail_Lane2 missing in {c_pos11.name}/system/blockMeshDict"
        print(f"[ok] Brand 'modifizierte Leitbleche' (L_trail=0.125 H): "
              f"{len(case_dirs_mod)} cases incl. Trail_Lane blocks")

    print("all offline smoke checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
