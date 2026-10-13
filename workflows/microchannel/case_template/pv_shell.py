"""thermalShell micro-channel -- shell (foil) temperature, one PNG.

Run by images.render_images() with pvbatch from the case dir:

    pvbatch --force-offscreen-rendering pv_shell.py \
        --T4 <time> --3T4 <time>

Reads postProcessing/shell/<time>/*.vtp (controlDict FO "shell": areaWrite
of Ts_outerShell on the whole wetted wall -- crest fillet, slanted flank,
trough) and writes shell_temperature.png: one window under one title, one
row per phase (top T/4 forward pass, bottom 3T/4 reverse pass), each the
wall from above,
tilted to the side and slightly from the inlet so it reads as 3D. The wall
is ~40x longer than wide, so the length (z) is compressed for display --
the tick labels give the true distance from the inlet and a note states the
factor. Shown mirrored about the symmetry plane x = 0 (the other
half-channel's wall, physically identical): that puts the flank's outer
face towards the camera with the inlet on the left. The trough (bottom)
lies under the flank and is mostly hidden from above. Colour range per row.
"""
import argparse
import glob
import os

from paraview.simple import (
    AssignViewToLayout, CellDatatoPointData, ColorBy, CreateLayout,
    CreateView, GetColorTransferFunction, GetScalarBar, Render,
    SaveScreenshot, Show, Text, Transform, XMLPolyDataReader,
    _DisableFirstRenderCameraReset,
)
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.numpy_interface import dataset_adapter as dsa
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--T4", default=None, help="time dir of the T/4 shell data")
ap.add_argument("--3T4", dest="T34", default=None,
                help="time dir of the 3T/4 shell data")
args = ap.parse_args()

PHASES = [("T4", args.T4, "T/4 (top)"), ("3T4", args.T34, "3T/4 (bottom)")]
FIELD = "Ts_outerShell"
OUT = "shell_temperature.png"
IMG_W = 1400          # image width [px]; the height follows the content
BAR_PX = 210          # left band for the colour bar [px]
TITLE_PX = 50         # title band above the top row [px]
ROW_TOP_PX = 55       # band above every wall: room for the colour-bar title
NOTE_PX = 30          # band for the compression note under the last row [px]
LABEL_FONT = 22       # tick labels [pt]
LENGTH_RATIO = 4.0    # displayed length / cross-section size (z compression)
# camera offset from the focal point (display coordinates: mirrored in x,
# z compressed): from above (+y), tilted towards the flank's outer side and
# from the inlet (-z), so the slanted wall faces the viewer and reads as 3D
TO_CAM = np.array([-1.0, 1.0, -0.7])
N_TICKS = 6           # distance labels at 0, l/5, .., l


def load(path):
    """Points (m) and the cell field of the shell."""
    r = vtkXMLPolyDataReader()
    r.SetFileName(path)
    r.Update()
    d = dsa.WrapDataObject(r.GetOutput())
    return np.asarray(d.Points), np.asarray(d.CellData[FIELD])


def shell_file(tname):
    files = sorted(glob.glob(os.path.join("postProcessing", "shell", tname,
                                          "*.vtp")))
    return files[0] if files else None


rows = [(key, tname, caption, shell_file(tname))
        for key, tname, caption in PHASES if tname is not None]
rows = [r for r in rows if r[3] is not None]
if not rows:
    raise SystemExit("[pv_shell] no shell data -- nothing to render")

# geometry (identical for every phase) and the shared camera / scale
pts, _ = load(rows[0][3])
ex, ey = np.ptp(pts[:, 0]), np.ptp(pts[:, 1])
z0, z1 = float(pts[:, 2].min()), float(pts[:, 2].max())
length = z1 - z0
sz = min(1.0, LENGTH_RATIO * max(ex, ey) / length)      # z compression
disp_pts = pts * [-1.0, 1.0, sz]                       # mirrored, z compressed

# z exactly horizontal on screen (inlet left, outlet right); up follows
to_cam = TO_CAM / np.linalg.norm(TO_CAM)
zdir = np.array([0.0, 0.0, 1.0])
right = zdir - zdir.dot(to_cam) * to_cam
right /= np.linalg.norm(right)
up = np.cross(to_cam, right)
focal = 0.5 * (disp_pts.min(axis=0) + disp_pts.max(axis=0))
rel = disp_pts - focal
w = np.ptp(rel @ right)
h = np.ptp(rel @ up)
# width-fit (half a tick label spare at either end), the row height follows:
# colour-bar-title band + wall + tick labels under it
m_per_px = w / (IMG_W - BAR_PX - 2 * 50)
tick_px = 1.6 * LABEL_FONT
wall_px = int(np.ceil(h / m_per_px))
row_px = ROW_TOP_PX + wall_px + int(tick_px) + 20
view_h = [TITLE_PX + row_px] + [row_px] * (len(rows) - 1)
view_h[-1] += NOTE_PX
img = [IMG_W, sum(view_h)]

_DisableFirstRenderCameraReset()
layout = CreateLayout("shell")
views = [CreateView("RenderView")]
AssignViewToLayout(view=views[0], layout=layout, hint=0)
if len(rows) > 1:
    layout.SplitVertical(0, view_h[0] / img[1])
    views.append(CreateView("RenderView"))
    AssignViewToLayout(view=views[1], layout=layout, hint=2)
layout.SetSize(*img)          # = ImageResolution: fonts stay unscaled

for row, (key, tname, caption, path) in enumerate(rows):
    view = views[row]
    view.Background = [1.0, 1.0, 1.0]
    view.UseColorPaletteForBackground = 0
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1
    vh = view_h[row]
    top_band = (TITLE_PX if row == 0 else 0) + ROW_TOP_PX
    bottom_band = vh - top_band - wall_px

    _, vals = load(path)
    lo, hi = float(vals.min()), float(vals.max())
    if hi - lo < 1e-6:                        # uniform field: still a bar
        lo, hi = lo - 0.005, hi + 0.005

    src = XMLPolyDataReader(FileName=[path])
    smooth = CellDatatoPointData(Input=src)   # smooth like the slices
    squash = Transform(Input=smooth)
    squash.Transform = "Transform"
    squash.Transform.Scale = [-1.0, 1.0, sz]
    disp = Show(squash, view)
    ColorBy(disp, ("POINTS", FIELD))
    lut = GetColorTransferFunction(f"{FIELD}_{key}")
    lut.ApplyPreset("Rainbow Uniform", True)
    lut.RescaleTransferFunction(lo, hi)
    lut.AutomaticRescaleRangeMode = "Never"
    disp.LookupTable = lut
    disp.Ambient = 0.25          # lift the shading on the steep flank

    bar = GetScalarBar(lut, view)
    bar.Visibility = 1
    bar.Title = "Ts (K)"
    bar.ComponentTitle = ""
    bar.TitleColor = [0.0, 0.0, 0.0]
    bar.LabelColor = [0.0, 0.0, 0.0]
    bar.TitleFontSize = 28
    bar.TitleBold = 1
    bar.LabelFontSize = 24
    bar.ScalarBarThickness = 24
    bar.HorizontalTitle = 1
    bar.AutomaticLabelFormat = 0
    bar.LabelFormat = "%-#.2f"
    bar.RangeLabelFormat = "%-#.2f"
    bar.AddRangeLabels = 1
    bar.UseCustomLabels = 1
    bar.CustomLabels = list(np.linspace(lo, hi, 5)[1:-1])
    bar.WindowLocation = "Any Location"
    bar.Orientation = "Vertical"
    bar.Position = [0.02, (bottom_band - 10) / vh]
    bar.ScalarBarLength = (wall_px + 10) / vh

    # true distance from the inlet, on a line just below the wall; the
    # billboard anchors the text's lower-left corner, so hang it by its height
    below = (rel @ up).min() - (0.25 + 1.3) * LABEL_FONT * m_per_px
    band = 0.02 * length * sz
    for i in range(N_TICKS):
        zf = i / (N_TICKS - 1)
        text = f"{1000.0 * (z0 + zf * length):.0f} mm"
        zd = (z0 + zf * length) * sz
        anchor = disp_pts[np.abs(disp_pts[:, 2] - zd) <= band].mean(axis=0)
        anchor = anchor + (below - (anchor - focal) @ up) * up
        half_w = 0.5 * 0.55 * LABEL_FONT * len(text) * m_per_px
        lbl = Text(Text=text)
        ld = Show(lbl, view)
        ld.TextPropMode = "Billboard 3D Text"
        ld.FontSize = LABEL_FONT
        ld.Color = [0.0, 0.0, 0.0]
        ld.BillboardPosition = list(anchor - half_w * right)

    if row == 0:                  # one title for the whole image
        head = Text(Text="Shell temperature: "
                         + ", ".join(r[2] for r in rows))
        hd = Show(head, view)
        hd.WindowLocation = "Upper Center"
        hd.FontSize = 34
        hd.Bold = 1
        hd.Color = [0.0, 0.0, 0.0]
    if row == len(rows) - 1:
        note = Text(Text=(f"length compressed {1.0 / sz:.0f}x for display"
                          if sz < 1.0 else "true scale"))
        nd = Show(note, view)
        nd.WindowLocation = "Lower Right Corner"
        nd.FontSize = 16
        nd.Color = [0.3, 0.3, 0.3]

    # wall placed in its band: right of the colour bar, between the top
    # bands and the tick labels
    centre = (focal
              - 0.5 * BAR_PX * m_per_px * right
              + 0.5 * (top_band - bottom_band) * m_per_px * up)
    view.CameraFocalPoint = list(centre)
    view.CameraPosition = list(centre + 10.0 * max(length * sz, ey) * to_cam)
    view.CameraViewUp = list(up)
    view.CameraParallelScale = 0.5 * vh * m_per_px
    Render(view)

SaveScreenshot(OUT, layout, SaveAllViews=1, ImageResolution=img,
               TransparentBackground=0)
print(f"[pv_shell] wrote {OUT} ({', '.join(r[1] for r in rows)})")
