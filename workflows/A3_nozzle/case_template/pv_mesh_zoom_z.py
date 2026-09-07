"""A3 nozzle -- normal-z mesh slice in crinkle mode (Polyederzellen sichtbar).

Crinkle-slice = no clipping at the plane, instead show whole cells whose
centroid is within slice-cell-thickness of the plane. Reveals the true
polyhedral cell structure (vs the triangulated polygons of a regular slice).

Default position: outlet plane (z = -(Lo + Lk + L_exit) ~ -26 mm). For the
default geometry with L_exit=0 and Lk=6, this matches the cone-of-cells at
the orifice exit. For the transition view (z ~ 0), edit the SLICE_Z below.

Output: pv_mesh_zoom_z.png
"""
import os
from paraview.simple import (
    OpenFOAMReader, Slice, Show, GetActiveViewOrCreate, SaveScreenshot,
    Render, UpdatePipeline, ResetCamera, _DisableFirstRenderCameraReset,
)

_DisableFirstRenderCameraReset()

SLICE_Z = -0.0258   # m -- just inside the outlet for the default geometry
                    # (Outlet bei z=-0.026; -25.8mm faengt die letzte Cell-Reihe ein).

foam_file = next(f for f in os.listdir() if f.endswith(".foam"))
print("[pv_mesh_zoom_z] reading %s" % foam_file)

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

slc = Slice(registrationName="slice_z_crinkle", Input=reader)
slc.SliceType = "Plane"
slc.SliceType.Origin = [0.0, 0.0, SLICE_Z]
slc.SliceType.Normal = [0.0, 0.0, 1.0]
slc.Crinkleslice = 1   # KEY: show whole cells, not triangulated cross-sections
UpdatePipeline(time=None, proxy=slc)

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1024, 1024]
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

# Camera straight along -z onto the orifice cross-section.
view.CameraPosition = [0.0, 0.0, 0.05]
view.CameraFocalPoint = [0.0, 0.0, SLICE_Z]
view.CameraViewUp = [0.0, 1.0, 0.0]
view.CameraParallelScale = 0.004    # ~ 4 mm half-width

Render(view)
SaveScreenshot("pv_mesh_zoom_z.png", view, ImageResolution=[1024, 1024])
print("[pv_mesh_zoom_z] saved pv_mesh_zoom_z.png at z=%g" % SLICE_Z)
