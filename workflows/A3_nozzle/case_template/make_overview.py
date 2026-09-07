#!/usr/bin/env python3
"""A3_nozzle -- overview.png + overview.html aus den Run-Parametern.

Generiert pro Case einen Konfigurations-Snapshot:

  case_<name>/overview.png    matplotlib-Bild (Tabelle + 2D-Halbschnitt-Skizze)
  case_<name>/overview.html   info.md + embedded PNG + DOE-Code-Banner

Aufruf (vom Worker oder standalone):

  python make_overview.py <case_dir> [--info-md <path>] [--params-json <path>]

Wenn --params-json fehlt, wird `<case_dir>/.overview_params.json` erwartet.
Format:

  {
    "project": "A3_zhang_baseline",
    "doe_code": "A1B2C2D2",       # optional, leer = individuelle Parameter
    "geom": {Di, Ds, Ls, alpha, Do, Lo, Lk, n_inlet, L_exit, R_exit},
    "mesh": {N_Do, bl_cells},
    "flow": {u_inlet, end_time, n_proc},
    "toggles": {bl_passes, zhang_bcs, disable_prefill, enable_ambient},
    "ambient": {L_amb, R_amb} | null
  }
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.lines import Line2D


# ---------------------------------------------------------------------------
# Zhang DOE level tables (Section 4, Table 5)
# ---------------------------------------------------------------------------
DOE_LEVELS = {
    "theta": ("A", [0.0, 10.0, 20.0, 30.0]),    # expansion angle [deg]
    "Do":    ("B", [4.0,  5.0,  6.0,  7.0]),    # orifice diameter [mm]
    "alpha": ("C", [30.0, 45.0, 60.0, 75.0]),   # contraction angle [deg]
    "Ds":    ("D", [9.0, 10.0, 11.0, 12.0]),    # swirl chamber diameter [mm]
}


def doe_letter(var: str, value: float) -> str:
    """Return e.g. 'B2' for var='Do', value=5.0; '-' if not on a Zhang level."""
    if var not in DOE_LEVELS:
        return "-"
    letter, levels = DOE_LEVELS[var]
    for i, lv in enumerate(levels, start=1):
        if abs(value - lv) < 0.1:    # tolerate 0.1 deg / 0.1 mm rounding from derived theta
            return f"{letter}{i}"
    return f"{letter}?"   # value not on a Zhang level


def theta_from_geom(Do: float, Lk: float, R_exit: float) -> float:
    """Zhang's expansion angle [deg], FULL opening angle (Section 4.3):

        R_exit = Do/2 + Lk * tan(theta/2)   =>   theta = 2 * atan((R_exit - Do/2) / Lk)

    A4 (theta=30 deg) with Do=4 mm, Lk=5 mm gives R_exit = 3.34 mm.
    """
    if Lk <= 0:
        return 0.0
    return 2.0 * math.degrees(math.atan((R_exit - 0.5 * Do) / Lk))


# ---------------------------------------------------------------------------
# Parameter table
# ---------------------------------------------------------------------------
def build_param_rows(p: dict) -> list[tuple]:
    """Return list of (label, description, value_str, doe_letter) rows."""
    g = p["geom"]
    m = p["mesh"]
    f = p["flow"]
    t = p["toggles"]
    a = p.get("ambient")
    theta = theta_from_geom(g["Do"], g["Lk"], g["R_exit"])

    rows: list[tuple] = []
    rows.append(("Project",      "Run-/Case-Name",           p["project"],          "-"))
    if p.get("doe_code"):
        rows.append(("DOE-Code",     "Zhang Taguchi L16",        p["doe_code"],         "-"))

    rows.append(("Geometry (Zhang Fig. 1b)", "", "", ""))
    rows.append(("θ (theta)",    "Expansionswinkel",         f"{theta:.1f}°",        doe_letter("theta", theta)))
    rows.append(("Do",           "Bohrungs-Ø",               f"{g['Do']:.2f} mm",    doe_letter("Do", g["Do"])))
    rows.append(("α (alpha)",    "Kontraktions-Halbwinkel",  f"{g['alpha']:.1f}°",   doe_letter("alpha", g["alpha"])))
    rows.append(("Ds",           "Drallkammer-Ø",            f"{g['Ds']:.2f} mm",    doe_letter("Ds", g["Ds"])))
    rows.append(("Di",           "Tangential-Inlet-Ø",       f"{g['Di']:.2f} mm",    "-"))
    rows.append(("Ls",           "Drallkammer-Länge",        f"{g['Ls']:.2f} mm",    "-"))
    rows.append(("Lo",           "Bohrungs-Länge",           f"{g['Lo']:.2f} mm",    "-"))
    rows.append(("Lk",           "Expansions-Länge",         f"{g['Lk']:.2f} mm",    "-"))
    rows.append(("R_exit",       "Auslass-Radius",           f"{g['R_exit']:.2f} mm", "-"))
    rows.append(("L_exit",       "Gerade nach Expansion",    f"{g['L_exit']:.2f} mm", "-"))
    rows.append(("n_inlet",      "Anzahl Tangentialeinläufe", f"{int(g['n_inlet'])}", "-"))

    rows.append(("Mesh", "", "", ""))
    rows.append(("N_Do",         "Zellen über Do",            f"{m['N_Do']}",         "-"))
    rows.append(("bl_cells",     "BL-Band-Zellen",            f"{m['bl_cells']}",     "-"))
    rows.append(("BL passes",    "refineWallLayer-Pässe",     f"{t['bl_passes']}",    "-"))

    rows.append(("Flow & BC", "", "", ""))
    rows.append(("U_inlet",      "Inlet-Geschwindigkeit",     f"{f['u_inlet']:.2f} m/s", "-"))
    rows.append(("endTime",      "Solver-Endzeit",            f"{f['end_time']:.4f} s",  "-"))
    rows.append(("n_proc",       "MPI-Ranks",                 f"{f['n_proc']}",         "-"))
    rows.append(("Zhang outlet", "p_rgh=0 + alpha inletOutlet 0", "ja" if t["zhang_bcs"] else "nein", "-"))
    rows.append(("no-prefill",   "alpha=0 überall (no setFields)", "ja" if t["disable_prefill"] else "nein", "-"))

    if a:
        rows.append(("Ambient", "", "", ""))
        rows.append(("L_amb",   "Ambient-Länge",             f"{a['L_amb']:.1f} mm",  "-"))
        rows.append(("R_amb",   "Ambient-Radius",            f"{a['R_amb']:.1f} mm",  "-"))

    return rows


# ---------------------------------------------------------------------------
# Sketch (Zhang Fig. 1b style, 2D half-cut)
# ---------------------------------------------------------------------------
def draw_sketch(ax, p: dict) -> None:
    """2D-Halbschnitt: Drallkammer oben, Bohrung mittig, Expansion unten.

    Convention: r = radial, z = axial (down). Origin at the orifice top
    (= z=0 = chamber bottom / contraction-end). Drawing goes downward
    (negative z) as in Zhang Fig. 1b.
    """
    g = p["geom"]
    Ds = g["Ds"]; Ls = g["Ls"]; Do = g["Do"]; alpha = g["alpha"]
    Lo = g["Lo"]; Lk = g["Lk"]; L_exit = g["L_exit"]; R_exit = g["R_exit"]
    Di = g["Di"]
    Rs = 0.5 * Ds
    Ro = 0.5 * Do
    Lc = (Rs - Ro) * math.tan(math.radians(alpha))   # contraction axial length
    theta = theta_from_geom(Do, Lk, R_exit)

    # z-Koordinaten (downward = negative z in our convention; here we keep z
    # positive going down for plot simplicity, and invert the y-axis below).
    z_chamber_top   = -(Ls + Lc)           # top of swirl chamber
    z_chamber_bot   = -Lc                  # contraction starts here
    z_orifice_top   = 0.0                  # orifice begins
    z_orifice_bot   = Lo                   # orifice ends
    z_expansion_bot = Lo + Lk              # expansion ends
    z_exit          = Lo + Lk + L_exit     # exit straight ends

    # Outline polygon (right half: r>=0)
    outline = [
        (0.0,        z_chamber_top),
        (Rs,         z_chamber_top),
        (Rs,         z_chamber_bot),
        (Ro,         z_orifice_top),
        (Ro,         z_orifice_bot),
        (R_exit,     z_expansion_bot),
        (R_exit,     z_exit),
        (0.0,        z_exit),
    ]
    poly_right = Polygon(outline, closed=True, facecolor="#cfe5ff",
                          edgecolor="#1f4f8b", linewidth=1.6, zorder=1)
    # Mirror left half (r negative)
    outline_left = [(-x, y) for x, y in outline]
    poly_left = Polygon(outline_left, closed=True, facecolor="#cfe5ff",
                         edgecolor="#1f4f8b", linewidth=1.6, zorder=1)
    ax.add_patch(poly_right)
    ax.add_patch(poly_left)

    # Symmetry axis
    ax.plot([0, 0], [z_chamber_top, z_exit], "k:", linewidth=0.7, alpha=0.5)

    # ----- Dimensions (right side) -----
    fs = 9
    # Ds (chamber diameter)
    y_ds = z_chamber_top - 0.15 * Ls
    ax.annotate("", xy=(-Rs, y_ds), xytext=(Rs, y_ds),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax.text(0, y_ds - 0.05 * Ls, f"Ds = {Ds:.1f} mm  (D{doe_letter('Ds', Ds)[1:]})",
            ha="center", va="top", fontsize=fs)

    # Ls
    x_ls = Rs + 0.8
    ax.annotate("", xy=(x_ls, z_chamber_top), xytext=(x_ls, z_chamber_bot),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax.text(x_ls + 0.2, 0.5 * (z_chamber_top + z_chamber_bot),
            f"Ls = {Ls:.1f} mm", ha="left", va="center", fontsize=fs)

    # alpha (contraction half-angle, at right contraction wall)
    cx, cy = Rs, z_chamber_bot
    ax.text(cx + 0.3, cy + 0.5 * Lc + 0.1,
            f"α = {alpha:.0f}°  ({doe_letter('alpha', alpha)})",
            ha="left", va="center", fontsize=fs, color="#a02000")

    # Do (orifice diameter)
    y_do = 0.5 * Lo
    ax.annotate("", xy=(-Ro, y_do), xytext=(Ro, y_do),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax.text(Ro + 0.25, y_do, f"Do = {Do:.1f} mm  ({doe_letter('Do', Do)})",
            ha="left", va="center", fontsize=fs)

    # Lo
    x_lo = Ro + 0.8
    ax.annotate("", xy=(x_lo, z_orifice_top), xytext=(x_lo, z_orifice_bot),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax.text(x_lo + 0.2, 0.5 * Lo, f"Lo = {Lo:.1f} mm", ha="left", va="center", fontsize=fs)

    # Lk + theta
    x_lk = R_exit + 0.8
    ax.annotate("", xy=(x_lk, z_orifice_bot), xytext=(x_lk, z_expansion_bot),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.9))
    ax.text(x_lk + 0.2, 0.5 * (z_orifice_bot + z_expansion_bot),
            f"Lk = {Lk:.1f} mm", ha="left", va="center", fontsize=fs)

    # theta-Label (Expansionswinkel) am Lampenschirm
    if Lk > 0 and theta > 0.5:
        mx = 0.5 * (Ro + R_exit) + 0.2
        my = 0.5 * (z_orifice_bot + z_expansion_bot)
        ax.text(mx, my, f"θ = {theta:.1f}°  ({doe_letter('theta', theta)})",
                ha="left", va="center", fontsize=fs, color="#a02000")
    else:
        ax.text(R_exit + 0.3, z_expansion_bot - 0.1,
                f"θ = {theta:.1f}°  ({doe_letter('theta', theta)})",
                ha="left", va="top", fontsize=fs, color="#a02000")

    # L_exit (if > 0)
    if L_exit > 0.05:
        x_le = R_exit + 1.8
        ax.annotate("", xy=(x_le, z_expansion_bot), xytext=(x_le, z_exit),
                    arrowprops=dict(arrowstyle="<->", color="grey", lw=0.7))
        ax.text(x_le + 0.2, 0.5 * (z_expansion_bot + z_exit),
                f"L_exit = {L_exit:.1f}", ha="left", va="center", fontsize=fs - 1, color="grey")

    # Di (tangential inlet, drawn as arrow into chamber wall)
    z_inlet = z_chamber_top + 0.3 * Ls
    ax.add_patch(FancyArrowPatch((Rs + 1.6, z_inlet), (Rs, z_inlet),
                                  arrowstyle="-|>", mutation_scale=12,
                                  color="#1f78b4", lw=1.4))
    ax.text(Rs + 1.8, z_inlet, f"Di = {Di:.1f} mm\n× {int(g['n_inlet'])}",
            ha="left", va="center", fontsize=fs - 1, color="#1f78b4")
    # mirror inlet arrow
    ax.add_patch(FancyArrowPatch((-(Rs + 1.6), z_inlet), (-Rs, z_inlet),
                                  arrowstyle="-|>", mutation_scale=12,
                                  color="#1f78b4", lw=1.4))

    # R_exit-Label
    ax.text(R_exit + 0.25, z_expansion_bot, f"R_exit={R_exit:.2f}",
            ha="left", va="bottom", fontsize=fs - 1, color="#555")

    # Axes
    ax.set_xlim(-(Rs + 4.5), Rs + 4.5)
    ax.set_ylim(z_exit + 1.5, z_chamber_top - 0.5 * Ls)   # inverted: z down
    ax.set_aspect("equal")
    ax.set_xlabel("r [mm]")
    ax.set_ylabel("z [mm] (Strömungsrichtung →)")
    ax.set_title("Geometry (Zhang Fig. 1b half-cut)", fontsize=11)
    ax.grid(True, alpha=0.2)


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------
def render_overview_png(p: dict, out_png: Path) -> None:
    fig = plt.figure(figsize=(16, 11), constrained_layout=False)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0],
                          left=0.03, right=0.98, top=0.95, bottom=0.04, wspace=0.08)
    ax_tbl = fig.add_subplot(gs[0, 0])
    ax_sk  = fig.add_subplot(gs[0, 1])

    # --- table ---
    ax_tbl.set_axis_off()
    rows = build_param_rows(p)
    cell_text = []
    cell_colors = []
    for label, descr, value, doe in rows:
        if descr == "" and value == "":
            cell_text.append([label, "", "", ""])
            cell_colors.append(["#e6e6e6"] * 4)
        else:
            cell_text.append([label, descr, value, doe])
            cell_colors.append(["#ffffff"] * 4)
    tbl = ax_tbl.table(cellText=cell_text,
                       colLabels=["Parameter", "Beschreibung", "Wert", "DOE"],
                       cellLoc="left", colLoc="left", loc="upper left",
                       colWidths=[0.22, 0.44, 0.22, 0.10],
                       cellColours=cell_colors)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.0, 1.55)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#bbbbbb")
        if r == 0:
            cell.set_facecolor("#1f4f8b")
            cell.set_text_props(color="white", weight="bold")
        if c == 3:
            cell.set_text_props(ha="center")
    ax_tbl.set_title(f"A3 Pressure-Swirl Nozzle — Run «{p['project']}»",
                      fontsize=12, weight="bold", loc="left")

    # --- sketch ---
    draw_sketch(ax_sk, p)

    fig.savefig(out_png, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# HTML wrapper
# ---------------------------------------------------------------------------
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="de"><head><meta charset="utf-8">
<title>A3 Run Overview — {project}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         max-width: 1100px; margin: 1.5em auto; padding: 0 1em; color: #222; }}
  h1, h2 {{ color: #1f4f8b; }}
  .doe-banner {{ background: #1f4f8b; color: white; padding: 0.8em 1em;
                 border-radius: 6px; font-size: 1.1em; margin: 1em 0; }}
  pre {{ background: #f4f4f4; padding: 0.6em; border-radius: 4px; overflow-x: auto; }}
  img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 6px;
        margin: 1em 0; }}
  .info-md {{ border-top: 2px solid #1f4f8b; padding-top: 1em; margin-top: 2em; }}
  table {{ border-collapse: collapse; }} th, td {{ padding: 4px 8px; }}
</style></head>
<body>
<h1>A3 Pressure-Swirl Nozzle — Run «{project}»</h1>
<div class="doe-banner">
  <strong>DOE-Code:</strong> {doe_or_dash} &nbsp;•&nbsp;
  <strong>θ</strong>={theta:.1f}° ({doe_a}) &nbsp;
  <strong>Do</strong>={Do:.1f} mm ({doe_b}) &nbsp;
  <strong>α</strong>={alpha:.0f}° ({doe_c}) &nbsp;
  <strong>Ds</strong>={Ds:.1f} mm ({doe_d})
</div>
<img src="overview.png" alt="Konfigurations-Snapshot">
<div class="info-md">
<h2>Workflow-Beschreibung</h2>
{info_html}
</div>
</body></html>
"""


def md_to_html(md: str) -> str:
    """Tiny Markdown -> HTML converter (h1/h2/list/code/bold/italic/link/table).

    Avoids the python-markdown dependency. Good enough for our info.md.
    """
    out = []
    lines = md.split("\n")
    i = 0
    in_list = False
    in_table = False
    in_code = False
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("```"):
            if in_code:
                out.append("</pre>"); in_code = False
            else:
                out.append("<pre>"); in_code = True
            i += 1; continue
        if in_code:
            out.append(_esc(ln)); i += 1; continue
        if ln.startswith("# "):
            out.append(f"<h2>{_inline(ln[2:])}</h2>")
        elif ln.startswith("## "):
            out.append(f"<h3>{_inline(ln[3:])}</h3>")
        elif ln.startswith("### "):
            out.append(f"<h4>{_inline(ln[4:])}</h4>")
        elif ln.lstrip().startswith("- "):
            if not in_list: out.append("<ul>"); in_list = True
            out.append(f"<li>{_inline(ln.lstrip()[2:])}</li>")
        elif ln.startswith("|") and "|" in ln[1:]:
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            if all(set(c) <= set("-: ") for c in cells):
                pass    # separator
            else:
                tag = "th" if (not in_table) else "td"
                if not in_table:
                    out.append("<table border='1' cellspacing='0'>"); in_table = True
                out.append("<tr>" + "".join(f"<{tag}>{_inline(c)}</{tag}>" for c in cells) + "</tr>")
        elif ln.strip() == "":
            if in_list: out.append("</ul>"); in_list = False
            if in_table: out.append("</table>"); in_table = False
            out.append("")
        else:
            out.append(f"<p>{_inline(ln)}</p>")
        i += 1
    if in_list: out.append("</ul>")
    if in_table: out.append("</table>")
    if in_code: out.append("</pre>")
    return "\n".join(out)


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _inline(s: str) -> str:
    import re
    s = _esc(s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", s)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', s)
    return s


def render_overview_html(p: dict, info_md_path: Path | None, out_html: Path) -> None:
    g = p["geom"]
    theta = theta_from_geom(g["Do"], g["Lk"], g["R_exit"])
    info_html = ""
    if info_md_path and info_md_path.exists():
        info_html = md_to_html(info_md_path.read_text(encoding="utf-8"))
    out_html.write_text(HTML_TEMPLATE.format(
        project=p["project"],
        doe_or_dash=p.get("doe_code") or "(individuelle Parameter)",
        theta=theta, doe_a=doe_letter("theta", theta),
        Do=g["Do"], doe_b=doe_letter("Do", g["Do"]),
        alpha=g["alpha"], doe_c=doe_letter("alpha", g["alpha"]),
        Ds=g["Ds"], doe_d=doe_letter("Ds", g["Ds"]),
        info_html=info_html,
    ), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--info-md", default=None,
                    help="Pfad zu info.md (default: <case_dir>/../info.md oder <case_dir>/info.md)")
    ap.add_argument("--params-json", default=None,
                    help="Pfad zu params JSON (default: <case_dir>/.overview_params.json)")
    args = ap.parse_args(argv)

    case = Path(args.case_dir)
    params_json = Path(args.params_json) if args.params_json else case / ".overview_params.json"
    if not params_json.exists():
        raise SystemExit(f"params JSON not found: {params_json}")
    p = json.loads(params_json.read_text(encoding="utf-8"))

    info_md = Path(args.info_md) if args.info_md else None
    if info_md is None:
        for cand in (case / "info.md", case.parent / "info.md",
                     Path(__file__).parent.parent / "info.md"):
            if cand.exists():
                info_md = cand; break

    out_png = case / "overview.png"
    out_html = case / "overview.html"
    render_overview_png(p, out_png)
    render_overview_html(p, info_md, out_html)
    print(f"wrote {out_png}")
    print(f"wrote {out_html}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
