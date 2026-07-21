"""A2 Leitbleche (Mother-Topologie) -- mesh + patch preview (top-down).

Runs after blockMesh + createBaffles, BEFORE the solver. Shows:
- Whole mesh (internal cells + every boundary patch) as Surface With Edges,
  light off-white surface with strong black cell edges -> the cell grid
  (resolution, refineWallLayer split, bend topology) is the focus.

Output: pv_meshView.png next to the .foam file.

Vane position is intentionally not over-emphasised here -- the 3-slice
iso view and the streamlines image cover that. The mesh-preview's job
is to catch broken meshes (wrong R/H, cells degenerated at the corner,
missing refinement) before the solver runs.
"""
import os
from paraview.simple import (
    OpenFOAMReader, Show, GetActiveViewOrCreate, SaveScreenshot,
    Render, UpdatePipeline, _DisableFirstRenderCameraReset,
)

_DisableFirstRenderCameraReset()

foam_file = next(f for f in os.listdir() if f.endswith(".foam"))
print(f"[pv_meshView] reading {foam_file}")

reader = OpenFOAMReader(registrationName=foam_file, FileName=foam_file)
reader.CaseType = "Reconstructed Case"
reader.SkipZeroTime = 0
reader.CellArrays = []
# OpenFOAMReader 'DecomposePolyhedra' property: present in PV 5.x, removed in
# PV 6.0+ where even hasattr raises NotSupportedException. Swallow both.
try:
    for attr in ("Decomposepolyhedra", "DecomposePolyhedra"):
        if hasattr(reader, attr):
            try: setattr(reader, attr, 0)
            except Exception: pass
            break
except Exception:
    pass

available = list(reader.MeshRegions.Available)
patch_regions = [r for r in available if r.startswith("patch/")]
reader.MeshRegions = ["internalMesh"] + patch_regions
UpdatePipeline(time=None, proxy=reader)
print(f"[pv_meshView] {len(patch_regions)} patches loaded")

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1280, 800]
view.Background = [1.0, 1.0, 1.0]
view.UseColorPaletteForBackground = 0
view.OrientationAxesVisibility = 1
view.CameraParallelProjection = 1

disp = Show(reader, view, "UnstructuredGridRepresentation")
disp.Representation = "Surface With Edges"
disp.AmbientColor = [0.95, 0.95, 0.95]
disp.DiffuseColor = [0.95, 0.95, 0.95]
disp.EdgeColor = [0.0, 0.0, 0.0]
disp.LineWidth = 1.0

view.ResetCamera()
focal = list(view.CameraFocalPoint)
view.CameraPosition = [focal[0], focal[1], focal[2] + 20.0]
view.CameraFocalPoint = focal
view.CameraViewUp = [0.0, 1.0, 0.0]
Render()

SaveScreenshot("pv_meshView.png", view, ImageResolution=[1280, 800],
               TransparentBackground=0)
print("[pv_meshView] wrote pv_meshView.png")
