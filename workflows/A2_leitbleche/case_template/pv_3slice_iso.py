"""A2 Leitbleche (Mother) -- isometric view with three orthogonal mid-axis slices."""
import os
from paraview.simple import (
    OpenFOAMReader, Slice, ExtractBlock, FeatureEdges, Tube, Show, ColorBy,
    GetColorTransferFunction, GetScalarBar, GetActiveViewOrCreate,
    SaveScreenshot, Render, UpdatePipeline, _DisableFirstRenderCameraReset,
)

_DisableFirstRenderCameraReset()

foam_file = next(f for f in os.listdir() if f.endswith(".foam"))
print(f"[pv_3slice_iso] reading {foam_file}")

reader = OpenFOAMReader(registrationName=foam_file, FileName=foam_file)
reader.CaseType = "Reconstructed Case"
reader.SkipZeroTime = 1
reader.CellArrays = ["U", "p"]
# Decomposepolyhedra was removed in ParaView 5.13. setattr in try/except --
# in PV 6.0.1 hasattr() itself raises NotSupportedException for retired props.
for attr in ("Decomposepolyhedra", "DecomposePolyhedra"):
    try:
        setattr(reader, attr, 0)
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
H = -2.0 * ymin if ymin < 0 else 1.0
x_mid = 0.5 * (xmin + xmax)
y_mid = 0.5 * (ymin + ymax)
z_mid = 0.5 * (zmin + zmax)
print(f"[pv_3slice_iso] mids: x={x_mid:.2f} y={y_mid:.2f} z={z_mid:.2f} H={H:.2f}")

internal_only = ExtractBlock(registrationName="internal", Input=reader)
internal_only.Selectors = ["/Root/internalMesh"]
UpdatePipeline(time=None, proxy=internal_only)

slices_spec = [
    ("slice_x", [x_mid, 0.0, 0.0], [1.0, 0.0, 0.0]),
    ("slice_y", [0.0, y_mid, 0.0], [0.0, 1.0, 0.0]),
    ("slice_z", [0.0, 0.0, z_mid], [0.0, 0.0, 1.0]),
]

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1280, 1000]
view.Background = [1.0, 1.0, 1.0]
view.UseColorPaletteForBackground = 0
view.OrientationAxesVisibility = 1

slice_disps = []
for name, origin, normal in slices_spec:
    s = Slice(registrationName=name, Input=internal_only)
    s.SliceType = "Plane"
    s.SliceType.Origin = origin
    s.SliceType.Normal = normal
    UpdatePipeline(time=None, proxy=s)
    d = Show(s, view, "GeometryRepresentation")
    d.Representation = "Surface"
    ColorBy(d, ("POINTS", "U", "Magnitude"))
    slice_disps.append(d)

u_lut = GetColorTransferFunction("U")
u_lut.ApplyPreset("Blue to Red Rainbow", True)
slice_disps[2].RescaleTransferFunctionToDataRange(True, False)
slice_disps[2].SetScalarBarVisibility(view, True)

sbar = GetScalarBar(u_lut, view)
sbar.TitleColor = [0.0, 0.0, 0.0]
sbar.LabelColor = [0.0, 0.0, 0.0]
sbar.Title = "|U| [m/s]"

# Vanes: zero-thickness baffles seen near-edge-on at iso angle barely render.
# Extract vane surfaces -> FeatureEdges -> Tube to make them visible from
# any viewing angle.
if vane_patches:
    vane_selectors = [f"//{vp.split("/")[-1]}" for vp in vane_patches]
    vane_block = ExtractBlock(registrationName="vanes", Input=reader)
    vane_block.Selectors = vane_selectors
    UpdatePipeline(time=None, proxy=vane_block)
    # Surface itself, opaque dark grey
    vane_surf = Show(vane_block, view, "GeometryRepresentation")
    vane_surf.Representation = "Surface"
    vane_surf.AmbientColor = [0.10, 0.10, 0.10]
    vane_surf.DiffuseColor = [0.10, 0.10, 0.10]
    # Plus tubed feature edges to give the vanes a clear silhouette
    edges = FeatureEdges(registrationName="vane_edges", Input=vane_block)
    edges.BoundaryEdges = 1
    edges.FeatureEdges = 1
    edges.NonManifoldEdges = 0
    edges.ManifoldEdges = 0
    UpdatePipeline(time=None, proxy=edges)
    tubes = Tube(registrationName="vane_tubes", Input=edges)
    tubes.Radius = 0.012 * H
    UpdatePipeline(time=None, proxy=tubes)
    tubes_disp = Show(tubes, view, "GeometryRepresentation")
    tubes_disp.Representation = "Surface"
    tubes_disp.AmbientColor = [0.0, 0.0, 0.0]
    tubes_disp.DiffuseColor = [0.0, 0.0, 0.0]

view.ResetCamera()
focal = list(view.CameraFocalPoint)
extent = max(xmax - xmin, ymax - ymin, zmax - zmin)
d = 2.5 * extent
view.CameraPosition = [focal[0] + d, focal[1] - 0.6*d, focal[2] + 0.7*d]
view.CameraFocalPoint = focal
view.CameraViewUp = [0.0, 0.0, 1.0]
view.CameraParallelProjection = 1
Render()

SaveScreenshot("pv_3slice_iso.png", view, ImageResolution=[1280, 1000],
               TransparentBackground=0)
print("[pv_3slice_iso] wrote pv_3slice_iso.png")
