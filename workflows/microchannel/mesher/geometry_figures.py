"""Documentation figures of the channel geometry: the true fluid domain
against the CFD vertical cut (h_profile.png) and the full-pitch diagram
(full_pitch.png). Both are drawn from true_geometry's construction.

matplotlib is imported lazily (Agg) so the geometry gate and the mesher
stay usable without it.
"""
import numpy as np

from true_geometry import solve_cfd_outline, solve_true_geometry, true_wall

N_SAMPLES = 280          # points per fillet arc


def _foil_faces(g, s):
    """The foil's two faces as closed-polygon pieces, both from the true
    wall of the solved construction ``g``.

    :return: ``(trough, crest, fluid, far)`` -- the sampled fluid faces of
        the trough (floor contact ``(b, 0)`` -> ``T_t``) and the crest
        (``T_c`` -> apex); ``fluid = (x, y)``, the fluid face closed at the
        symmetry axis; ``far = (x, y)``, the face ``s`` out of the fluid,
        starting under the axis on the base sheet's far side.
    """
    trough_w, _, crest_w = true_wall(g)
    trough, crest = trough_w.sample(N_SAMPLES), crest_w.sample(N_SAMPLES)
    trough_o = trough_w.offset(-s).sample(N_SAMPLES)
    crest_o = crest_w.offset(-s).sample(N_SAMPLES)
    fluid = (np.r_[0.0, trough[:, 0], crest[:, 0]],
             np.r_[0.0, trough[:, 1], crest[:, 1]])
    far = (np.r_[0.0, g["Ct"][0], trough_o[:, 0], crest_o[:, 0]],
           np.r_[-s, -s, trough_o[:, 1], crest_o[:, 1]])
    return trough, crest, fluid, far


def plot_h_profile(a, b, s, R, out_path, h=0.05, dpi=150, suptitle=None):
    """True-fluid-vs-CFD figure for the vertical-cut approximation (README
    sec. 3.1), a parameterised port of the parametrization repo's
    study/05_domain_description/h_profile.py."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = solve_cfd_outline(a, b, s, R, h)
    T_c, T_t = g["T_c"], g["T_t"]
    Q, x_h = g["Q"], g["x_h"]
    trough, crest, (tx, ty), (ox, oy) = _foil_faces(g, s)

    # hypotenuse as one polyline, floor contact (b,0) -> apex (0,a_f);
    # y is strictly increasing along it, so y = h cuts it exactly once
    hyp_x = np.r_[trough[:, 0], T_c[0], crest[:, 0]]
    hyp_y = np.r_[trough[:, 1], T_c[1], crest[:, 1]]
    k = int(np.searchsorted(hyp_y, h))   # first hypotenuse sample at/above y=h

    # half thickness: the mid-surface of each foil, i.e. the s/2 the CFD
    # shell models. Base sheet and ribbon are separate foils (they only touch
    # at the trough), so their mid-lines are drawn as separate pieces.
    trough_w, _, crest_w = true_wall(g)
    crest_m = crest_w.offset(-0.5 * s).sample(N_SAMPLES)
    trough_m = trough_w.offset(-0.5 * s).sample(N_SAMPLES)
    mid_lines = [(np.array([0.0, g["Ct"][0]]), np.array([-0.5 * s, -0.5 * s])),
                 (np.r_[trough_m[:, 0], crest_m[:, 0]],
                  np.r_[trough_m[:, 1], crest_m[:, 1]])]

    # zoom window tracks the tangency height/depth so small-fillet and
    # large-fillet cases both stay framed (nominal keeps the original box)
    zx0 = min(b - 0.52, T_t[0] - 0.06)
    zy1 = max(0.42, T_t[1] + 0.10)
    views = [((-0.08, b + 0.16, -0.12, a + 0.14), "True fluid vs CFD domain"),
             ((zx0, b + 0.10, -0.06, zy1), "Zoomed view of the cut")]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6))
    for ax, ((x0, x1, y0, y1), ttl) in zip(axes, views):
        # the solid: everything between the outer edge and the true fluid
        ax.fill(np.r_[ox, tx[::-1]], np.r_[oy, ty[::-1]],
                color="#d6d6d6", zorder=1, label="Solid")
        ax.fill(tx, ty, color="#dbe9f6", zorder=2, label="True fluid domain")
        ax.plot(ox, oy, color="#8a8a8a", lw=2.0, zorder=4)
        for i, (mx, my) in enumerate(mid_lines):
            ax.plot(mx, my, color="#6b6b6b", lw=1.1, ls="-.", zorder=4,
                    label="Half thickness s/2 (CFD shell)" if i == 0 else None)
        ax.plot(tx, ty, color="#1f5fa8", lw=2.2, zorder=4)
        # approximated wall, in three pieces so each line starts its dash
        # phase with ink (one path would let a dash gap swallow the short
        # vertical). The crest is shared with the true wall, so the dashes
        # overlay it there.
        ax.plot([0.0, x_h], [0.0, 0.0],
                color="#b02418", lw=2.2, ls="--", zorder=5,
                label="CFD domain boundary")
        ax.plot([x_h, x_h], [0.0, h], color="#b02418", lw=2.2, ls="--", zorder=5)
        ax.plot(np.r_[Q[0], hyp_x[k:]], np.r_[h, hyp_y[k:]],
                color="#b02418", lw=2.2, ls="--", zorder=5)
        sx = np.r_[x_h, Q[0], hyp_x[:k][::-1], x_h]
        sy = np.r_[0.0, h, hyp_y[:k][::-1], 0.0]
        ax.fill(sx, sy, color="#f6d3a8", zorder=3, alpha=.95,
                label="Cut off for CFD")
        # construction: the y = h line and the two points it produces
        ax.plot([0.0, Q[0]], [h, h], color="#555", lw=1.0, ls=":", zorder=6,
                label="Cut off height h")
        ax.plot(*Q, marker="o", ms=4.5, color="#b02418", zorder=7)
        ax.plot([x_h], [0.0], marker="o", ms=4.5, color="#b02418", zorder=7)
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.grid(alpha=.25)
        ax.set_xlabel("x [mm]")
        ax.set_title(ttl, fontsize=11)
    axes[0].set_ylabel("y [mm]")
    axes[0].legend(fontsize=8.5, loc="upper right")
    axes[1].legend(fontsize=8.5, loc="upper right")

    fig.suptitle(suptitle or
                 "Domain Description: True fluid vs CFD approximation",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def half_pitch_polygons(a, b, s, R):
    """Polygons of the half-pitch section x in [0, b]: lower fluid, solid
    (base sheet + ribbon), upper sheet and the inverted (upper) fluid.
    Port of the parametrization repo's
    study/05_domain_description/overall_geomtry.py half_pitch."""
    g = solve_true_geometry(a, b, s, R)
    _, _, (fx, fy), (ox, oy) = _foil_faces(g, s)
    # the far face's arcs, (b, 0)-side first: the upper fluid runs them back
    far_x, far_y = ox[2:], oy[2:]

    # solid (base sheet + ribbon), closed: outer face forward + lower
    # fluid boundary reversed
    sx = np.r_[ox, fx[::-1]]
    sy = np.r_[oy, fy[::-1]]

    # upper sheet strip, y in [a, a + s]
    shx = np.array([0.0, b, b, 0.0])
    shy = np.array([a, a, a + s, a + s])

    # inverted (upper) fluid, closed polygon
    ux = np.r_[far_x[::-1], b]
    uy = np.r_[far_y[::-1], a]

    return dict(fx=fx, fy=fy, sx=sx, sy=sy, shx=shx, shy=shy,
                ux=ux, uy=uy, ox=ox, oy=oy, g=g)


def plot_full_pitch(a, b, s, R, out_path, dpi=150,
                    title="Full pitch geometry"):
    """Full-pitch diagram of the TRUE geometry (fluid + solid), port of the
    parametrization repo's study/05_domain_description/overall_geomtry.py
    draw_full_pitch.

    The half-pitch section is mirrored about x = b to complete one pitch,
    then about x = 0 for the neighbouring pitch: x in [-2b, 2b], two full
    pitches, three foil layers and two fluid channels."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    SOLID, FLUID = "#9e9e9e", "#dbe9f6"
    hp = half_pitch_polygons(a, b, s, R)
    a_f, a_s = hp["g"]["a_f"], a + s

    placements = [
        lambda x: x,          # original half, x in [0, b]
        lambda x: 2 * b - x,  # mirror about x = b
        lambda x: -x,         # mirror about x = 0
        lambda x: x - 2 * b,  # complete pitch mirrored about x = 0
    ]

    fig, ax = plt.subplots(figsize=(15, 6))
    for i, t in enumerate(placements):
        first = (i == 0)
        ax.fill(t(hp["sx"]), hp["sy"], color=SOLID, zorder=1,
                label="Solid" if first else None)
        ax.fill(t(hp["shx"]), hp["shy"], color=SOLID, zorder=1)
        ax.fill(t(hp["fx"]), hp["fy"], color=FLUID, zorder=2,
                label="Fluid" if first else None)
        ax.fill(t(hp["ux"]), hp["uy"], color=FLUID, zorder=2)

        # black boundary of the solid, physical edges only (no legend)
        ax.plot(t(hp["fx"][1:]), hp["fy"][1:], color="k", lw=1.3, zorder=4)
        ax.plot(t(hp["ox"][2:]), hp["oy"][2:], color="k", lw=1.3, zorder=4)

    # the ribbon end cuts at the figure edges, x = +/-2b
    for xe in (-2 * b, 2 * b):
        ax.plot([xe, xe], [a_f, a], color="k", lw=1.3, zorder=4)
    # base sheet and upper sheet rectangles
    for y0, y1 in ((-s, 0.0), (a, a_s)):
        ax.fill([-2 * b, 2 * b, 2 * b, -2 * b], [y0, y0, y1, y1],
                fill=False, edgecolor="k", lw=1.3, zorder=4)

    m = 0.06 * max(a, b)
    ax.set_xlim(-2 * b - m, 2 * b + m)
    ax.set_ylim(-s - 0.05 * a, a_s + 0.08 * a)
    ax.set_aspect("equal")
    ax.grid(alpha=.2)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.legend(fontsize=14, loc="lower left", framealpha=.95)
    ax.set_title(title, fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
