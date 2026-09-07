"""A3 nozzle -- Zhang Fig.4 style: normal-x slice with alpha/U/p triptych.

Three vertical sub-images stacked into one PNG. Same slice plane (x=0), latest
time. Colormap "Blue to Red Rainbow" -- the classic ParaView rainbow LUT.

  panel 1: alpha.water     [0..1]     -- water film vs air-core
  panel 2: |U|             [0..U_max] -- velocity magnitude
  panel 3: p_rgh           [auto]     -- gauge pressure (rho*g*h removed)

Output: pv_figure4.png  (single PNG, 1200x1500, 3 rows)
"""
import os
import sys
import math
from paraview.simple import (
    OpenFOAMReader, Slice, Calculator, Show, Hide, GetActiveViewOrCreate,
    SaveScreenshot, Render, UpdatePipeline, GetColorTransferFunction,
    GetScalarBar, AssignViewToLayout, _DisableFirstRenderCameraReset,
)

_DisableFirstRenderCameraReset()

foam_file = next(f for f in os.listdir() if f.endswith(".foam"))
print("[pv_figure4] reading %s" % foam_file)

reader = OpenFOAMReader(registrationName=foam_file, FileName=foam_file)
reader.CaseType = "Reconstructed Case"
reader.SkipZeroTime = 0
reader.CellArrays = ["alpha.water", "U", "p_rgh"]
try:
    for attr in ("Decomposepolyhedra", "DecomposePolyhedra"):
        if hasattr(reader, attr):
            try: setattr(reader, attr, 0)
            except Exception: pass
            break
except Exception:
    pass
reader.MeshRegions = ["internalMesh"]
# Move to latest non-zero time. If the solver never ran (e.g. mpirun aborted at
# init, two concurrent 128-rank jobs on cfdtools), the case only has t=0 (or
# nothing), so U/p_rgh aren't populated -- the Calculator(mag(U)) below would
# error out with "Undefined symbol: 'U'". Skip cleanly in that case.
tsteps = reader.TimestepValues if hasattr(reader, "TimestepValues") else []
solver_times = [t for t in tsteps if t > 0.0]
if not solver_times:
    print("[pv_figure4] no solver output (only t=0 or empty case) -- skipping")
    sys.exit(0)
UpdatePipeline(time=solver_times[-1], proxy=reader)

slc = Slice(registrationName="slice_x0", Input=reader)
slc.SliceType = "Plane"
slc.SliceType.Origin = [0.0, 0.0, 0.0]
slc.SliceType.Normal = [1.0, 0.0, 0.0]
slc.Triangulatetheslice = 0
UpdatePipeline(time=None, proxy=slc)

# |U| via Calculator (raw U is vector).
calc_U = Calculator(registrationName="Umag", Input=slc)
calc_U.ResultArrayName = "Umag"
calc_U.Function = "mag(U)"
UpdatePipeline(time=None, proxy=calc_U)


def _render_panel(input_proxy, field, png, lo=None, hi=None):
    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [800, 500]
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
    bar.ComponentTitle = ""
    bar.WindowLocation = "Lower Right Corner"
    bar.Orientation = "Vertical"
    # camera along +x; show xz mirror (the slice is yz-plane).
    view.CameraPosition = [0.05, 0.0, -0.008]
    view.CameraFocalPoint = [0.0, 0.0, -0.008]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraParallelScale = 0.020
    Render(view)
    SaveScreenshot(png, view, ImageResolution=[800, 500])
    Hide(input_proxy, view)
    return view


# alpha.water 0..1
_render_panel(slc,    "alpha.water", "_fig4_alpha.png", lo=0.0, hi=1.0)
# |U|
_render_panel(calc_U, "Umag",        "_fig4_U.png")
# p_rgh
_render_panel(slc,    "p_rgh",       "_fig4_p.png")

# Stitch the 3 PNGs into one figure with PIL if available -- else leave the
# 3 separate files and a marker PNG for downstream zip-pickup.
try:
    from PIL import Image
    imgs = [Image.open(p) for p in ("_fig4_alpha.png", "_fig4_U.png", "_fig4_p.png")]
    w = max(im.size[0] for im in imgs)
    h = sum(im.size[1] for im in imgs)
    out = Image.new("RGB", (w, h), "white")
    y = 0
    for im in imgs:
        out.paste(im, (0, y))
        y += im.size[1]
    out.save("pv_figure4.png")
    print("[pv_figure4] saved pv_figure4.png (3 panels stitched)")
except Exception as e:
    # Fallback: rename alpha panel as the main output.
    try:
        os.replace("_fig4_alpha.png", "pv_figure4.png")
    except OSError:
        pass
    print("[pv_figure4] PIL not available (%s) -- using alpha-only panel" % e)
