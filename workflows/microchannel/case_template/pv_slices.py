"""thermalShell micro-channel -- cross-section slices, one PNG per phase.

Run by images.render_images() with pvbatch from the case dir:

    pvbatch --force-offscreen-rendering pv_slices.py \
        --T4 <time> --3T4 <time>

Reads postProcessing/slices/<time>/slice*.vtp (controlDict FO "slices":
z-planes at 0.2/0.4/0.6/0.8 l) and writes fluid_slices_<key>.png. Each image
is one window, two rows of the same four slices at their true z, in an
isometric view:
  top    |U| [m/s]
  bottom T   [K]
Each slice is labelled above its apex with its distance from the inlet
(z = 0). The domain is half a channel (symmetryPlane at x = 0), shown
as simulated (not mirrored).
Colour ranges are per image (the four slices of that time).
"""
import argparse
import glob
import os

from paraview.simple import (
    AssignViewToLayout, ColorBy, CreateLayout, CreateView, GetColorTransferFunction,
    GetScalarBar, Render, ResetSession, SaveScreenshot, Show, Text,
    XMLPolyDataReader, _DisableFirstRenderCameraReset,
)
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--T4", default=None, help="time dir of the T/4 slices")
ap.add_argument("--3T4", dest="T34", default=None,
                help="time dir of the 3T/4 slices")
args = ap.parse_args()

PHASES = [("T4", args.T4, "T/4 (forward pass)"),
          ("3T4", args.T34, "3T/4 (reverse pass)")]
FIELDS = [("U", "Magnitude", "U [m/s]"), ("T", None, "T [K]")]
IMG_W = 1400          # image width [px]; the height follows the content
CONTENT_W = 0.85      # share of the width for the slices (rest: colour bar)
TITLE_PX = 60         # title band above the top row [px]
LABEL_PX = 70         # right margin for the last slice's distance label [px]
LABEL_FONT = 24       # distance labels [pt]
MIN_ROW_PX = 320      # keeps the colour bar readable on squat geometries
MAX_ROW_PX = 700      # caps the image height on tall / short geometries
# Camera: down the channel axis from the inlet side, yawed to the side and
# pitched from above so the slices stand in a gently rising staggered row.
# The angles adapt to the geometry (cross-section size vs slice spacing),
# so every parameter set gets the same look:
STEP_W = 2.0          # projected step between slices / slice width
MIN_STEP_W = 1.15     # .. packed closest (wide, flat sections)
MAX_STEP_W = 12.0     # .. spread furthest (narrow, tall sections)
RISE_H = 0.17         # projected rise per slice / slice height
MAX_YAW = np.radians(65.0)   # beyond this the slices turn edge-on

def load(path):
    """Points (m) and the U-magnitude / T point arrays of one slice."""
    r = vtkXMLPolyDataReader()
    r.SetFileName(path)
    r.Update()
    d = dsa.WrapDataObject(r.GetOutput())
    U = np.asarray(d.PointData["U"])
    return (np.asarray(d.Points),
            {"U": np.linalg.norm(U, axis=1), "T": np.asarray(d.PointData["T"])})


def render(key, tname, caption):
    files = sorted(glob.glob(os.path.join("postProcessing", "slices", tname,
                                          "*.vtp")))
    if not files:
        print(f"[pv_slices] no slices in {tname}, skipping {key}")
        return
    data = [load(f) for f in files]
    data_files = sorted(zip(data, files), key=lambda df: df[0][0][:, 2].mean())
    pts_all = np.vstack([d[0][0] for d in data_files])
    height = pts_all[:, 1].max() - pts_all[:, 1].min()
    t = float(tname)

    zs = [pts[:, 2].mean() for (pts, _), _f in data_files]
    first = data_files[0][0][0]
    ex, ey = np.ptp(first[:, 0]), np.ptp(first[:, 1])   # one cross-section
    dz = (zs[-1] - zs[0]) / max(1, len(zs) - 1)          # slice spacing
    focal = 0.5 * (pts_all.min(axis=0) + pts_all.max(axis=0))
    rel = pts_all - focal

    def fit(step_w):
        """Camera axes for a projected slice step of step_w slice widths,
        and the scale [m/px] / row height [px] that fit the image."""
        yaw = min(np.arctan(step_w * ex / dz), MAX_YAW)
        pitch = np.arctan(RISE_H * ey / (dz * np.cos(yaw)))
        to_cam = np.array([-np.sin(yaw) * np.cos(pitch), np.sin(pitch),
                           -np.cos(yaw) * np.cos(pitch)])
        up = np.array([0.0, 1.0, 0.0]) - np.sin(pitch) * to_cam
        up /= np.linalg.norm(up)
        right = np.cross(-to_cam, up)             # screen right
        w = np.ptp(rel @ right)
        h = np.ptp(rel @ up) + 0.2 * height       # + the labels on top
        # Width-fit: the slices span CONTENT_W of the width (the rest is the
        # colour bar; the last label sticks out right of its apex, keep
        # LABEL_PX free). The row height follows, so no vertical slack --
        # clamped, and a tall row is height-fit instead.
        m_per_px = 1.03 * w / (CONTENT_W * IMG_W - LABEL_PX)
        natural_px = int(np.ceil(1.08 * h / m_per_px))
        row_px = natural_px
        if row_px > MAX_ROW_PX:
            row_px = MAX_ROW_PX
            m_per_px = 1.08 * h / row_px
        row_px = max(MIN_ROW_PX, row_px)
        step_px = dz * abs(np.array([0.0, 0.0, 1.0]) @ right) / m_per_px
        return to_cam, up, right, m_per_px, row_px, step_px, natural_px

    # Tune the step: spread the slices out when the row would be too tall
    # (narrow sections), pack them closer when it would be too flat (wide
    # sections) -- but never so close that neighbouring labels overlap.
    label_px = 0.55 * LABEL_FONT * len("00.0 mm")
    step_w = STEP_W
    for _ in range(8):
        to_cam, up, right, m_per_px, row_px, step_px, natural_px = fit(step_w)
        if natural_px > 1.05 * MAX_ROW_PX and step_w < MAX_STEP_W:
            step_w = min(MAX_STEP_W, step_w * natural_px / MAX_ROW_PX)
        elif step_px < 1.2 * label_px:
            step_w *= 1.2 * label_px / step_px
        elif (natural_px < 0.95 * MIN_ROW_PX and step_w > MIN_STEP_W
              and step_px > 1.3 * label_px):
            step_w = max(MIN_STEP_W, step_w * natural_px / MIN_ROW_PX,
                         step_w * 1.25 * label_px / step_px)
        else:
            break
    view_h = [row_px + TITLE_PX, row_px]   # only the top row has the title
    img = [IMG_W, sum(view_h)]

    _DisableFirstRenderCameraReset()     # keep the cameras set below
    # one window, split top (cell 1) / bottom (cell 2)
    layout = CreateLayout(f"slices_{key}")
    views = [CreateView("RenderView")]
    AssignViewToLayout(view=views[0], layout=layout, hint=0)
    layout.SplitVertical(0, view_h[0] / img[1])
    views.append(CreateView("RenderView"))
    AssignViewToLayout(view=views[1], layout=layout, hint=2)
    layout.SetSize(*img)          # = ImageResolution: fonts stay unscaled

    # content centred in each row (labels included), shifted right clear of
    # the colour bar; the top row's centre moves up by half the title band
    base_focal = (focal + 0.1 * height * up
                  - 0.5 * ((1.0 - CONTENT_W) * IMG_W - LABEL_PX)
                  * m_per_px * right)

    for row, (field, comp, title) in enumerate(FIELDS):
        view = views[row]
        view.Background = [1.0, 1.0, 1.0]
        view.UseColorPaletteForBackground = 0
        view.OrientationAxesVisibility = 0
        view.CameraParallelProjection = 1

        vals = np.concatenate([d[0][1][field] for d in data_files])
        lo, hi = float(vals.min()), float(vals.max())
        lut = GetColorTransferFunction(f"{field}_{key}")
        lut.ApplyPreset("Rainbow Uniform", True)
        lut.RescaleTransferFunction(lo, hi)
        lut.AutomaticRescaleRangeMode = "Never"
        if comp:
            lut.VectorMode = "Magnitude"

        for (pts, _), f in data_files:
            src = XMLPolyDataReader(FileName=[f])
            disp = Show(src, view)
            ColorBy(disp, ("POINTS", field, comp) if comp
                    else ("POINTS", field))
            disp.LookupTable = lut
            # distance from the inlet, above the apex (highest point)
            apex = pts[np.argmax(pts[:, 1])]
            z_mm = 1000.0 * pts[:, 2].mean()
            text = f"{z_mm:.1f} mm"
            lbl = Text(Text=text)
            ld = Show(lbl, view)
            ld.TextPropMode = "Billboard 3D Text"
            ld.FontSize = LABEL_FONT
            # the billboard anchors the text's lower-left corner: shift it
            # left by half the (estimated) text width to centre on the apex
            half_w = 0.5 * 0.55 * ld.FontSize * len(text) * m_per_px
            anchor = (np.array([apex[0], apex[1] + 0.06 * height, apex[2]])
                      - half_w * right)
            ld.BillboardPosition = list(anchor)
            ld.Color = [0.0, 0.0, 0.0]

        bar = GetScalarBar(lut, view)
        bar.Visibility = 1
        bar.Title = title
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
        # 5 labels: the two range ends + 3 evenly spaced in between
        bar.UseCustomLabels = 1
        bar.CustomLabels = list(np.linspace(lo, hi, 5)[1:-1])
        bar.WindowLocation = "Any Location"
        bar.Orientation = "Vertical"
        # same pixel size in both rows (the top one is taller by the title)
        bar.Position = [0.02, 0.08 * row_px / view_h[row]]
        bar.ScalarBarLength = 0.78 * row_px / view_h[row]

        if row == 0:                  # one title for the whole image
            head = Text(Text=f"Velocity and Temperature field at {caption}")
            hd = Show(head, view)
            hd.WindowLocation = "Upper Center"
            hd.FontSize = 36
            hd.Bold = 1
            hd.Color = [0.0, 0.0, 0.0]

        cam_focal = base_focal + 0.5 * (view_h[row] - row_px) * m_per_px * up
        view.CameraFocalPoint = list(cam_focal)
        view.CameraPosition = list(cam_focal
                                   + 10.0 * np.ptp(pts_all[:, 2]) * to_cam)
        view.CameraViewUp = list(up)
        view.CameraParallelScale = 0.5 * view_h[row] * m_per_px
        Render(view)

    Render()
    out = f"fluid_slices_{key}.png"
    SaveScreenshot(out, layout, SaveAllViews=1, ImageResolution=img,
                   TransparentBackground=0)
    print(f"[pv_slices] wrote {out} (t = {tname})")


for key, tname, caption in PHASES:
    if tname is None:
        continue
    render(key, tname, caption)
    ResetSession()
