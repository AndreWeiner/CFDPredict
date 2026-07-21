"""A2 Leitbleche (Mother) -- mid-slice streamlines, diagonal seed line."""
import math
import os
from paraview.simple import (
    OpenFOAMReader, Slice, Line, StreamTracerWithCustomSource, ExtractBlock,
    Transform, Show, ColorBy, GetColorTransferFunction, GetScalarBar,
    GetActiveViewOrCreate, SaveScreenshot, Render, UpdatePipeline,
    _DisableFirstRenderCameraReset,
)

_DisableFirstRenderCameraReset()

foam_file = next(f for f in os.listdir() if f.endswith(".foam"))
print(f"[pv_streamlines] reading {foam_file}")

reader = OpenFOAMReader(registrationName=foam_file, FileName=foam_file)
reader.CaseType = "Reconstructed Case"
reader.SkipZeroTime = 1
reader.CellArrays = ["U", "p"]
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
vane_patches = [r for r in patch_regions if "/vane" in r]
reader.MeshRegions = ["internalMesh"] + patch_regions
UpdatePipeline(time=None, proxy=reader)

bounds = reader.GetDataInformation().GetBounds()
xmin, xmax, ymin, ymax, zmin, zmax = bounds
z_mid = 0.5 * (zmin + zmax)
# A2 geometry: pre-bend channel at y in [-H/2, +H/2]
H = -2.0 * ymin
R = xmax + ymin
r_inner = R - 0.5 * H
r_outer = R + 0.5 * H
ang = math.radians(-45.0)
inner_mid = (r_inner * math.cos(ang), R + r_inner * math.sin(ang), z_mid)
outer_mid = (r_outer * math.cos(ang), R + r_outer * math.sin(ang), z_mid)
print(f"[pv_streamlines] H={H:.3f} R={R:.3f} seed {inner_mid} -> {outer_mid}")

internal_only = ExtractBlock(registrationName="internal", Input=reader)
internal_only.Selectors = ["/Root/internalMesh"]
UpdatePipeline(time=None, proxy=internal_only)

slice_ = Slice(registrationName="z_mid_backdrop", Input=internal_only)
slice_.SliceType = "Plane"
slice_.SliceType.Origin = [0.0, 0.0, z_mid]
slice_.SliceType.Normal = [0.0, 0.0, 1.0]

# In-plane streamlines: integrate ON the z_mid slice with Surface Streamlines
# (the velocity is projected onto the cut plane), so the lines stay in the
# plane instead of weaving out of it through the Dean secondary flow -- the
# clean Brand Fig. 7.16 look, and no z-lift/occlusion hack needed. Seed from a
# diagonal line across the channel, on the plane.
seed = Line(registrationName="seed", Point1=list(inner_mid), Point2=list(outer_mid))
seed.Resolution = 40

streamlines = StreamTracerWithCustomSource(registrationName="streamlines",
                                           Input=slice_, SeedSource=seed)
streamlines.Vectors = ["POINTS", "U"]
streamlines.SurfaceStreamlines = 1
streamlines.IntegrationDirection = "BOTH"
streamlines.IntegratorType = "Runge-Kutta 4-5"
streamlines.MaximumStreamlineLength = 12.0 * (xmax - xmin + ymax - ymin)
streamlines.InitialStepLength = 0.01
streamlines.MaximumSteps = 8000
UpdatePipeline(time=None, proxy=streamlines)

# Backdrop: a grey copy of the slice pushed slightly behind (-z, away from the
# +z camera) so the coplanar in-plane streamlines render cleanly on top.
backdrop = Transform(registrationName="backdrop", Input=slice_)
backdrop.Transform = "Transform"
backdrop.Transform.Translate = [0.0, 0.0, -0.02 * (zmax - zmin)]
UpdatePipeline(time=None, proxy=backdrop)

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1280, 1000]
view.Background = [1.0, 1.0, 1.0]
view.UseColorPaletteForBackground = 0
view.OrientationAxesVisibility = 1
view.CameraParallelProjection = 1

slice_disp = Show(backdrop, view, "GeometryRepresentation")
slice_disp.Representation = "Surface"
slice_disp.AmbientColor = [0.94, 0.94, 0.94]
slice_disp.DiffuseColor = [0.94, 0.94, 0.94]

sl_disp = Show(streamlines, view, "GeometryRepresentation")
sl_disp.Representation = "Surface"
sl_disp.LineWidth = 2.0
sl_disp.RenderLinesAsTubes = 1
ColorBy(sl_disp, ("POINTS", "U", "Magnitude"))
u_lut = GetColorTransferFunction("U")
u_lut.ApplyPreset("Blue to Red Rainbow", True)
sl_disp.RescaleTransferFunctionToDataRange(True, False)
sl_disp.SetScalarBarVisibility(view, True)

sbar = GetScalarBar(u_lut, view)
sbar.TitleColor = [0.0, 0.0, 0.0]
sbar.LabelColor = [0.0, 0.0, 0.0]
sbar.Title = "|U| [m/s]"

if vane_patches:
    vane_selectors = [f"//{vp.split("/")[-1]}" for vp in vane_patches]
    vane_block = ExtractBlock(registrationName="vanes", Input=reader)
    vane_block.Selectors = vane_selectors
    vane_disp = Show(vane_block, view, "GeometryRepresentation")
    vane_disp.Representation = "Surface"
    vane_disp.AmbientColor = [0.15, 0.15, 0.15]
    vane_disp.DiffuseColor = [0.15, 0.15, 0.15]

# Camera: ResetCamera first to fit, then SET ViewUp last (ResetCamera can
# clobber the up vector if called after).
view.ResetCamera()
focal = list(view.CameraFocalPoint)
view.CameraPosition = [focal[0], focal[1], focal[2] + 20.0]
view.CameraFocalPoint = focal
# Brand-style: pre-bend points -X visually, post-bend rotated 90° to lay
# along the down direction. World +X axis appears as image LEFT (inlet at
# right side of image, flow going left), world +Y appears as image DOWN
# (post-bend going downward).
view.CameraViewUp = [-1.0, 0.0, 0.0]
# parallel-scale chosen to fit the full L-shape
view.CameraParallelScale = 0.5 * max(xmax - xmin, ymax - ymin) * 1.05
Render()

SaveScreenshot("pv_streamlines.png", view, ImageResolution=[1280, 1000],
               TransparentBackground=0)
print("[pv_streamlines] wrote pv_streamlines.png")
