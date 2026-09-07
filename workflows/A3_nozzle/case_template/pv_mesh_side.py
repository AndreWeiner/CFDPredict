"""A3 nozzle -- mesh side view (full domain).

Surface-with-edges rendering of the whole nozzle from the side (camera looking
along +y, slice normal-y at y=0). Lets the user check overall topology +
cell distribution from inlet at top to outlet at bottom.

Output: pv_mesh_side.png next to the .foam file.
"""
import os
from paraview.simple import (
    OpenFOAMReader, Slice, Show, GetActiveViewOrCreate, SaveScreenshot,
    Render, UpdatePipeline, ResetCamera, _DisableFirstRenderCameraReset,
)

_DisableFirstRenderCameraReset()

foam_file = next(f for f in os.listdir() if f.endswith(".foam"))
print("[pv_mesh_side] reading %s" % foam_file)

reader = OpenFOAMReader(registrationName=foam_file, FileName=foam_file)
reader.CaseType = "Reconstructed Case"
reader.SkipZeroTime = 0
reader.CellArrays = []
try:
    for attr in ("Decomposepolyhedra", "DecomposePolyhedra"):
        if hasattr(reader, attr):
            try: setattr(reader, attr, 0)
            except Exception: pass
            break
except Exception:
    pass
reader.MeshRegions = ["internalMesh"]
UpdatePipeline(time=None, proxy=reader)

slc = Slice(registrationName="slice_y0", Input=reader)
slc.SliceType = "Plane"
slc.SliceType.Origin = [0.0, 0.0, 0.0]
slc.SliceType.Normal = [0.0, 1.0, 0.0]
slc.Triangulatetheslice = 0  # keep polygons -> Polyederzellen sichtbar
UpdatePipeline(time=None, proxy=slc)

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [800, 1280]
view.Background = [1.0, 1.0, 1.0]
view.UseColorPaletteForBackground = 0
view.OrientationAxesVisibility = 1
view.CameraParallelProjection = 1

disp = Show(slc, view, "UnstructuredGridRepresentation")
disp.Representation = "Surface With Edges"
disp.AmbientColor = [0.95, 0.95, 0.95]
disp.DiffuseColor = [0.95, 0.95, 0.95]
disp.EdgeColor = [0.0, 0.0, 0.0]
disp.LineWidth = 0.5

ResetCamera(view)
# Camera looks along +y -> view xz-plane. Zoom to fit.
view.CameraPosition = [0.0, 0.05, -0.008]
view.CameraFocalPoint = [0.0, 0.0, -0.008]
view.CameraViewUp = [0.0, 0.0, 1.0]
view.CameraParallelScale = 0.020

Render(view)
SaveScreenshot("pv_mesh_side.png", view, ImageResolution=[800, 1280])
print("[pv_mesh_side] saved pv_mesh_side.png")
