#!/usr/bin/env python3
"""Offline smoke for A2_leitbleche (Mother): build the 6-tuple schema
into a series-dir (as the streamlit app does with inject_schema:true),
run the worker with do_source_of=False, do_run=False, and assert that
the case dirs are built with the right ingredients.

OpenFOAM is NOT required. Run from inside the workflow dir:
    cd workflows/A2_leitbleche && python test_build_offline.py
"""
from __future__ import annotations

import copy
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import py_A2_leitbleche as worker

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

        # 1) empty bend (n_vanes=0) -- baseline, no createBaffles needed
        series = make_series(
            tmp, "empty_bend",
            **{"number of vanes": "0",
               "cells across channel ny": "20",
               "cells in z (nz)": "1"})
        worker.run_worker(str(series), do_source_of=False, do_run=False)
        case_dirs = [p for p in series.iterdir()
                     if p.is_dir() and p.name.startswith("case_")]
        assert len(case_dirs) == 1, \
            f"expected 1 case dir, got {[c.name for c in case_dirs]}"
        c0 = case_dirs[0]
        assert (c0 / "system" / "blockMeshDict").is_file()
        assert not (c0 / "system" / "createBafflesDict").is_file(), \
            "n_vanes=0 must not produce createBafflesDict"
        assert (series / "command_finished").is_file()
        print(f"[ok] empty bend: {c0.name}")

        # 2) 3-vane default config -- createBaffles + topoSet present
        series2 = make_series(
            tmp, "three_vanes",
            **{"number of vanes": "3",
               "cells across channel ny": "20",
               "cells in z (nz)": "1"})
        worker.run_worker(str(series2), do_source_of=False, do_run=False)
        case_dirs2 = [p for p in series2.iterdir()
                      if p.is_dir() and p.name.startswith("case_")]
        assert len(case_dirs2) == 1, \
            f"expected 1 case dir, got {[c.name for c in case_dirs2]}"
        c2 = case_dirs2[0]
        assert (c2 / "system" / "createBafflesDict").is_file()
        assert (c2 / "system" / "topoSetDict").is_file()
        bd = (c2 / "system" / "createBafflesDict").read_text(encoding="utf-8")
        # 3 vane interfaces -> 3 faceZone blocks
        n_groups = bd.count("type            faceZone;")
        assert n_groups == 3, f"expected 3 baffle groups, got {n_groups}"
        # 0/U should have vane{1..3}_master and _slave entries (expanded from "vane.*")
        u0 = (c2 / "0" / "U").read_text(encoding="utf-8")
        for i in (1, 2, 3):
            for side in ("master", "slave"):
                tag = f"vane{i}_{side}"
                assert tag in u0, f"0/U missing entry for {tag}"
        print(f"[ok] 3-vane default: {c2.name}, {n_groups} baffle groups, "
              f"0/U has vane{{1..3}}_{{master,slave}}")

        # 3) vane_ext > 0 -- downstream extension via boxToFace
        series3 = make_series(
            tmp, "with_ext",
            **{"number of vanes": "2",
               "vane_ext1 [m]": "1.0",
               "vane_ext2 [m]": "0.5",
               "cells across channel ny": "20",
               "cells in z (nz)": "1"})
        worker.run_worker(str(series3), do_source_of=False, do_run=False)
        case_dirs3 = [p for p in series3.iterdir()
                      if p.is_dir() and p.name.startswith("case_")]
        assert len(case_dirs3) == 1, \
            f"expected 1 case dir, got {[c.name for c in case_dirs3]}"
        ts = (case_dirs3[0] / "system" / "topoSetDict").read_text(
            encoding="utf-8")
        # boxToFace bounds should reflect R + ext = 1.5 + 1.0 = 2.5 for vane 1
        assert "2.5" in ts or "2.500" in ts, \
            f"topoSetDict should reference y_hi=R+ext=2.5; not found"
        print(f"[ok] vane_ext=1.0+0.5: {case_dirs3[0].name}, boxToFace bounds present")

        # 4) 3D mode (nz>1) strips frontAndBack BC
        series4 = make_series(
            tmp, "three_d",
            **{"number of vanes": "1",
               "cells across channel ny": "12",
               "cells in z (nz)": "8"})
        worker.run_worker(str(series4), do_source_of=False, do_run=False)
        case_dirs4 = [p for p in series4.iterdir()
                      if p.is_dir() and p.name.startswith("case_")]
        u3d = (case_dirs4[0] / "0" / "U").read_text(encoding="utf-8")
        assert "frontAndBack" not in u3d, \
            "0/U for 3D case must not contain frontAndBack block"
        print(f"[ok] 3D mode: {case_dirs4[0].name}, frontAndBack stripped from 0/U")

        # 5) validation: out-of-range radius raises
        series5 = make_series(
            tmp, "bad_radius",
            **{"number of vanes": "1",
               "vane_r1 [m]": "2.5",   # outside (R-H/2, R+H/2) = (1.0, 2.0)
               "cells across channel ny": "12",
               "cells in z (nz)": "1"})
        try:
            worker.run_worker(str(series5), do_source_of=False, do_run=False)
            # Worker catches exceptions and writes 9_error.txt -- check it
            err = series5 / "progress" / "9_error.txt"
            # In the no-except path, run_worker reraises -- but the
            # __main__ in py_A2_leitbleche catches them. Direct call from
            # test code reraises, so this branch shouldn't be hit. If it
            # is, fail loudly:
            assert err.is_file(), \
                "out-of-range radius should have errored"
        except ValueError as e:
            assert "must lie strictly inside" in str(e), \
                f"unexpected ValueError: {e}"
            print(f"[ok] out-of-range radius rejected: {e}")

        # 6) validation: wrong ordering raises
        series6 = make_series(
            tmp, "bad_order",
            **{"number of vanes": "3",
               "vane_r1 [m]": "1.3",
               "vane_r2 [m]": "1.5",
               "vane_r3 [m]": "1.7",
               "cells across channel ny": "12",
               "cells in z (nz)": "1"})
        try:
            worker.run_worker(str(series6), do_source_of=False, do_run=False)
            assert False, "wrong ordering should have errored"
        except ValueError as e:
            assert "strictly decreasing" in str(e), \
                f"unexpected ValueError: {e}"
            print(f"[ok] wrong ordering rejected: {e}")

    print("all offline smoke checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
