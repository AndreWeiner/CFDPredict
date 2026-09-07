"""A3 nozzle -- mesh zoom on Drallkammer<->bore transition (normal-x slice).

Surface-with-edges rendering on slice x=0 (camera along +x). Zoomed on the
transition region z ~ -2..+2 mm where the orifice connects to the swirl
chamber bottom. This is where the BL refinement split (wall_top buffer +
wall_straight refined) is most visible -- and where Lesson-038 placed the
48 negative-volume cells when refineWallLayer was applied to the full
outerWall without the buffer.

Output: pv_mesh_zoom_x.png
"""
import os
from paraview.simple import (
    OpenFOAMReader, Slice, Show, GetActiveViewOrCreate, SaveScreenshot,
    Render, UpdatePipeline, ResetCamera, _DisableFirstRenderCameraReset,
)

_DisableFirstRenderCameraReset()

foam_file = next(f for f in os.listdir() if f.endswith(".foam"))
print("[pv_mesh_zoom_x] reading %s" % foam_file)

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

slc = Slice(registrationName="slice_x0", Input=reader)
slc.SliceType = "Plane"
slc.SliceType.Origin = [0.0, 0.0, 0.0]
slc.SliceType.Normal = [1.0, 0.0, 0.0]
slc.Triangulatetheslice = 0
UpdatePipeline(time=None, proxy=slc)

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1280, 1280]
view.Background = [1.0, 1.0, 1.0]
view.UseColorPaletteForBackground = 0
view.OrientationAxesVisibility = 1
view.CameraParallelProjection = 1

disp = Show(slc, view, "UnstructuredGridRepresentation")
disp.Representation = "Surface With Edges"
disp.AmbientColor = [0.95, 0.95, 0.95]
disp.DiffuseColor = [0.95, 0.95, 0.95]
disp.EdgeColor = [0.0, 0.0, 0.0]
disp.LineWidth = 1.0

# Zoom on Drallkammer<->bore transition (z ~ -2..+3 mm, y ~ -4..+4 mm).
view.CameraPosition = [0.05, 0.0, 0.0]
view.CameraFocalPoint = [0.0, 0.0, 0.0]
view.CameraViewUp = [0.0, 0.0, 1.0]
view.CameraParallelScale = 0.0035   # ~7 mm vertical extent

Render(view)
SaveScreenshot("pv_mesh_zoom_x.png", view, ImageResolution=[1280, 1280])
print("[pv_mesh_zoom_x] saved pv_mesh_zoom_x.png")
