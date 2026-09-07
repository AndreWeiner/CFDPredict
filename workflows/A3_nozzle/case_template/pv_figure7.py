"""A3 nozzle -- Zhang Fig.7 style: z-normal slice at outlet, alpha + U.

Two side-by-side panels stacked vertically. Outlet plane (z ~ -25.8 mm for
the default geometry). Colormap Blue to Red Rainbow.

  panel 1: alpha.water     [0..1]     -- film cross-section (annulus around air-core)
  panel 2: |U|             [auto]     -- velocity magnitude in the outlet section

Output: pv_figure7.png
"""
import os
import sys
from paraview.simple import (
    OpenFOAMReader, Slice, Calculator, Show, Hide, GetActiveViewOrCreate,
    SaveScreenshot, Render, UpdatePipeline, GetColorTransferFunction,
    GetScalarBar, _DisableFirstRenderCameraReset,
)

_DisableFirstRenderCameraReset()

SLICE_Z = -0.0258   # m -- just inside the outlet for default geometry

foam_file = next(f for f in os.listdir() if f.endswith(".foam"))
print("[pv_figure7] reading %s" % foam_file)

reader = OpenFOAMReader(registrationName=foam_file, FileName=foam_file)
reader.CaseType = "Reconstructed Case"
reader.SkipZeroTime = 0
reader.CellArrays = ["alpha.water", "U"]
try:
    for attr in ("Decomposepolyhedra", "DecomposePolyhedra"):
        if hasattr(reader, attr):
            try: setattr(reader, attr, 0)
            except Exception: pass
            break
except Exception:
    pass
reader.MeshRegions = ["internalMesh"]
# Skip cleanly if the solver never ran (e.g. mpirun aborted at init due to
# concurrent 128-rank load on cfdtools). U is not in the case -> the
# Calculator(mag(U)) below would error with "Undefined symbol: 'U'".
tsteps = reader.TimestepValues if hasattr(reader, "TimestepValues") else []
solver_times = [t for t in tsteps if t > 0.0]
if not solver_times:
    print("[pv_figure7] no solver output (only t=0 or empty case) -- skipping")
    sys.exit(0)
UpdatePipeline(time=solver_times[-1], proxy=reader)

slc = Slice(registrationName="slice_z_outlet", Input=reader)
slc.SliceType = "Plane"
slc.SliceType.Origin = [0.0, 0.0, SLICE_Z]
slc.SliceType.Normal = [0.0, 0.0, 1.0]
slc.Triangulatetheslice = 0
UpdatePipeline(time=None, proxy=slc)

calc_U = Calculator(registrationName="Umag", Input=slc)
calc_U.ResultArrayName = "Umag"
calc_U.Function = "mag(U)"
UpdatePipeline(time=None, proxy=calc_U)


def _render_panel(input_proxy, field, png, lo=None, hi=None):
    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [800, 800]
    view.Background = [1.0, 1.0, 1.0]
    view.UseColorPaletteForBackground = 0
    view.OrientationAxesVisibility = 1
    view.CameraParallelProjection = 1
    disp = Show(input_proxy, view, "GeometryRepresentation")
    disp.Representation = "Surface"
    disp.ColorArrayName = ["CELLS", field]
    lut = GetColorTransferFunction(field)
    lut.ApplyPreset("Blue to Red Rainbow", True)
    if lo is not None and hi is not None:
        lut.RescaleTransferFunction(lo, hi)
    else:
        disp.RescaleTransferFunctionToDataRange(True, False)
    disp.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    bar.Title = field
    bar.WindowLocation = "Lower Right Corner"
    view.CameraPosition = [0.0, 0.0, 0.05]
    view.CameraFocalPoint = [0.0, 0.0, SLICE_Z]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraParallelScale = 0.004
    Render(view)
    SaveScreenshot(png, view, ImageResolution=[800, 800])
    Hide(input_proxy, view)


_render_panel(slc,    "alpha.water", "_fig7_alpha.png", lo=0.0, hi=1.0)
_render_panel(calc_U, "Umag",        "_fig7_U.png")

try:
    from PIL import Image
    imgs = [Image.open(p) for p in ("_fig7_alpha.png", "_fig7_U.png")]
    w = sum(im.size[0] for im in imgs)
    h = max(im.size[1] for im in imgs)
    out = Image.new("RGB", (w, h), "white")
    x = 0
    for im in imgs:
        out.paste(im, (x, 0))
        x += im.size[0]
    out.save("pv_figure7.png")
    print("[pv_figure7] saved pv_figure7.png (alpha + Umag stitched)")
except Exception as e:
    try:
        os.replace("_fig7_alpha.png", "pv_figure7.png")
    except OSError:
        pass
    print("[pv_figure7] PIL not available (%s) -- using alpha-only panel" % e)
