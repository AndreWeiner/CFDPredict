# Authoring a Streamlit-Portal workflow

A workflow is a self-contained directory under `workflows/` in the Streamlit
app tree. The portal discovers it automatically via its `workflow.yaml`
manifest and renders a UI page from `interface.json`.

This document is the reference for the manifest fields and the directory
layout. When you add a new workflow, follow the structure below.

---

## 1. Directory layout

```
workflows/
└── myWorkflow/
    ├── workflow.yaml                 # MANIFEST — discovered by the portal
    ├── interface.json                # 6-tuple form schema (required)
    ├── py_myWorkflow.py              # entry script (subprocess target)
    ├── 0_template_<solver>/          # OpenFOAM template case (optional)
    │   ├── 0/, constant/, system/
    │   └── *.foam                    # ParaView marker file
    ├── pv_*.py                       # pvbatch render scripts (optional)
    ├── foam_dictionary.py            # shared utilities (per-workflow copy
    ├── myUtils.py / pyPipe.py        # is fine until we factor out a libs/)
    ├── pipeGeo.stl                   # default geometry (optional)
    └── interfaceUnrolled2.json       # GENERATED at runtime — do not commit
```

Generated/transient files that should NOT be checked in:
`__pycache__/`, `finished/`, `interfaceUnrolled2.json`, any `*.foam`-marker
files in run dirs (the one in `0_template_*/` is fine).

---

## 2. `workflow.yaml` field reference

| Field | Type | Default | Purpose |
|---|---|---|---|
| `title` | str | dirname | UI title shown on the workflow page |
| `description` | str | "" | Caption below the title (one or two sentences) |
| `entry_script` | str | required | Path to the worker script, relative to the workflow dir. Spawned via `subprocess.Popen([python_exec, entry_script, series_dir])` |
| `python_exec` | str | `sys.executable` of the streamlit process | Full path to the Python interpreter for the worker. Use `/usr/bin/python3` if your worker needs system packages (reportlab, paraview, foamCore wrappers). Empty string = inherit |
| `input_filename` | str | `input.json` | Filename written into the job's series dir with the form values. Most of our workflows use `interface.json` |
| `schema_file` | str | `interface.json` | Which file holds the 6-tuple form schema |
| `inject_schema` | bool | false | If `true`, values are written back into the 6-tuple schema (for `unroll_dict`); if `false`, a flat values dict is saved |
| `version` | str | `default` | Tag stored in the job log alongside the workflow name (e.g. `v2406`, `demo`) |
| `uploads` | bool | false | Show the file-uploader widget |
| `upload_extensions` | list[str] | `[stl, zip, tar, gz, tgz, obj, step, stp]` | Allowed file extensions when `uploads: true` |
| `max_concurrent_per_user` | int | `1` | Soft-lock limit. Warns + requires override checkbox when this user already has ≥ this many running jobs of this workflow. `0` = disabled. Admins always bypass |

Future fields we'll likely need (not yet implemented):
- `openfoam_version`: which OF tree to source for the worker (`v2406` vs `v2506`)
- `default_cores`: cores hint, for resource-aware scheduling
- `est_wallclock_min`: rough estimate for queue UI / billing

---

## 3. The 6-tuple schema (`interface.json`)

Each form field is a list of exactly six elements:

```
"field name": [position, default, type, options, selector, tooltip]
```

| Slot | Meaning |
|---|---|
| `position` | int — sort order on the form (sparse numbers OK: 10, 20, 30…) |
| `default` | any — initial value shown in the widget |
| `type` | str or list[str] — widget kind (see below) |
| `options` | depends on type — extra widget config |
| `selector` | dict or null — only meaningful for dropdown types; contains nested sub-schemas |
| `tooltip` | str — hover help text, may be empty |

### Type values

| `type` | Widget | `options` use | `selector` use |
|---|---|---|---|
| `"string"` | `st.text_input` | list[str] of example values shown in tooltip | — |
| `"text_area"` | `st.text_area` | int height in pixels (default 120) | — |
| `"integer"` | `st.number_input` (step=1) | — | — |
| `"float"` | `st.number_input` (`format="%g"`) | — | — |
| `"boolean"` | `st.checkbox` | — | — |
| `"valid_path"` | `st.text_input` | — | — |
| `"valid_file"` | informative banner | — | — |
| `"separator"` | `st.markdown("#### name")` heading | — | — |
| `"runPythonScript"` | caption only, value passes through | — | — |
| `list[str]` | `st.selectbox` with the list as choices | — | dict `{choice: [sub_name, sub_schema]}` for nested forms per choice |

### Dropdown + nested sub-schema example

```json
"project type": [
  1, "single",
  ["single", "serie", "optimization"],
  "single",
  {
    "single":       ["singleCoeffs",       { …sub-schema… }],
    "serie":        ["serieCoeffs",        { …sub-schema… }],
    "optimization": ["optimizationCoeffs", { …sub-schema… }]
  },
  "Choose how many runs to make"
]
```

The portal renders the chosen choice's sub-schema in a bordered container,
nested as `values["singleCoeffs"]["geometry definition"]…`.

---

## 4. Worker contract

`entry_script` is invoked as:

```
<python_exec> <entry_script> <series_dir>
```

`<series_dir>` is `JOBS_DIR / "<timestamp>_<project_name>"`. The series dir
has been created and contains:

- `<input_filename>` — JSON of the form values (or 6-tuple-injected schema
  if `inject_schema: true`)
- `progress/` — empty subdir, the worker writes `*_info.txt` and `*.png`
  here for the live log and the VTK viewer
- `uploads/` — present and populated only if `uploads: true` and the user
  attached files
- `owner` — username (one line)
- `workflow` — `<wf_name>\t<version>` (one line)

The worker must:

1. Read `<series_dir>/<input_filename>`
2. Do its work
3. Write progress via `progress/<step>_info.txt` so the live log shows it
4. On success, **touch `<series_dir>/command_finished`** as the very last step
5. Honour kill: if `<series_dir>/command_kill` appears, exit cleanly and
   touch `command_finished` (this marks the job as "killed" in the UI)

PDF reports, plots, and packaged output go into `<series_dir>/results/` so
the download-link helper picks them up.

---

## 5. End-to-end checklist for a new workflow

- [ ] Pick a short, lowercase name (`myWorkflow`)
- [ ] Create `workflows/myWorkflow/` with `workflow.yaml` and `interface.json`
- [ ] Write `py_myWorkflow.py` that follows the worker contract above
- [ ] Test locally: `python py_myWorkflow.py /tmp/test_series/` against a
      manually-prepared series dir with a sample `interface.json`
- [ ] Restart streamlit (`sudo systemctl restart streamlit.service`); the
      portal auto-discovers workflows on each request
- [ ] Grant ACL: an admin must add a row to `user_workflows` for non-admin
      users to see the workflow
- [ ] Set `max_concurrent_per_user: 1` unless the workflow is genuinely
      lightweight and parallel-safe
- [ ] Smoke-test from the UI with the test user
