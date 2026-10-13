"""Geometry of the vertical-cut CFD domain: the cut height, the eight
outline/block points, and the band's inner edge.

Pure geometry -- no cell counts. The wall itself comes from
``true_geometry.solve_cfd_outline`` as a list of pieces (Q to the apex);
everything here is built from those pieces and their offsets.
"""
import math

import numpy as np

from true_geometry import Arc, flank_line, solve_cfd_outline

# The cut height as a FRACTION of the smaller of the channel height `a` and
# the half-wavelength `b`. 
H_CUT_FRAC = 0.05

# Backstop cap on p5's horizontal offset, as a multiple of the band depth
# h. Over the gate-passing box the offset (h/sin(theta) on a flank cut)
# never exceeds 6.36 h.
P5_MAX_OFFSET_RATIO = 8.0

# Hand-over clearance for the inner edge's corners I_t / I_c: each must
# clear its bounding line (y = h for I_t, the axis x = 0 for I_c) by this
# multiple of the band depth, or the edge has no room to turn there and
# feeds the core slivers. When it does not, the corner walks along the
# flank's offset in HANDOVER_STATIONS steps until it does. With p5 placed
# on the wall's own offset curve the walk should almost never fire, and
# when it does it gives up exact offset over the stretch it bridges.
HANDOVER_CLEARANCE = 0.5
HANDOVER_STATIONS = 40


def cut_height(a, b):
    """The cut height ``h = H_CUT_FRAC * min(a, b)``, in mm. One resolver
    for every consumer -- the mesher, build_case's per-case figure and the
    README figures -- so the cut the mesh is built on and the cut the
    figures draw cannot drift apart."""
    return H_CUT_FRAC * min(float(a), float(b))


def solve_cut_geometry(shape, h):
    """The eight outline/block points of the vertical-cut approximation::

        p0 (0, 0)      p1 (x_h, 0)    p2 (x_h, h)    p3 (0, a-s)
        p4 (0, h)      p5 on y = h, where the wall's offset h crosses it
        p6 (0, a-s-h)  p7 below p5 on the floor

    ``p5`` is NOT ``p2 - (h, 0)``: the band's depth is measured
    perpendicular to the wall, so p5 sits where the wall's own inward
    offset crosses ``y = h`` -- the concentric arc ``r_t + h`` on a trough
    cut, the flank's parallel line otherwise (an ``h/sin(theta)`` step).
    That puts p5 exactly on the inner band edge's offset curve, so the edge
    needs no bridging spline there. (Placing it by ``h/sin(theta)`` on a
    trough cut leaves it 0.0085 off the arc, and the bridge S-bends with
    62x the arc's curvature.)

    :param shape: ``a``, ``b``, ``s``, ``R`` [mm].
    :param h: the cut height (:func:`cut_height`) [mm].
    :return: geometry dict -- ``p0``..``p7``, ``h``, ``a_f``, ``x_h``,
        ``p5_offset``, the ``piece`` that owns the cut, ``wall`` (its
        pieces, Q to the apex) and ``solved`` (the full construction).
    """
    a, b = float(shape["a"]), float(shape["b"])
    s, R = float(shape["s"]), float(shape["R"])
    h = float(h)

    g = solve_cfd_outline(a, b, s, R, h)
    a_f, x_h = g["a_f"], g["x_h"]
    where = f"(a={a}, b={b}, s={s}, R={R})"

    # The block layout needs room for the band on both sides: p5 must stay
    # right of the axis, and p6 must stay above p4.
    if x_h <= 2.0 * h:
        raise ValueError(
            f"the cut lands too close to the axis for a band of depth "
            f"{h}: x_h={x_h:.4f} {where}")
    if a_f - h <= 2.0 * h:
        raise ValueError(
            f"the fluid apex a-s={a_f:.4f} leaves no room for a band of "
            f"depth {h} (a={a}, s={s})")

    try:
        p5_offset = x_h - g["wall"][0].offset(h).x_at_y(h)
    except ValueError as exc:
        raise ValueError(f"the wall's offset never reaches the cut height "
                         f"h={h} {where}: {exc}") from exc
    if p5_offset <= 0.0:
        raise ValueError(
            f"the wall's offset at the cut puts p5 at or right of p2 {where}")
    p5_offset = min(p5_offset, P5_MAX_OFFSET_RATIO * h, 0.5 * x_h)

    return dict(a=a, b=b, s=s, R=R, h=h, a_f=a_f, x_h=x_h,
                piece=g["piece"], solved=g, wall=g["wall"],
                p0=(0.0, 0.0), p1=(x_h, 0.0), p2=(x_h, h), p3=(0.0, a_f),
                p4=(0.0, h), p5=(x_h - p5_offset, h), p6=(0.0, a_f - h),
                p7=(x_h - p5_offset, 0.0),
                p5_offset=p5_offset)


def flank(geo):
    """The true flank ``T_t -> T_c`` (whether or not the cut kept all of it)."""
    g = geo["solved"]
    return flank_line(g, g["T_t"])


def inner_band_edge(geo, outer_fr, depth):
    """The band's inner edge ``p5 -> p6``: the TRUE offset of the wall.

    Built piece by piece from the wall's own offsets, so the band's
    PERPENDICULAR depth is ``depth`` wherever an offset exists:

    * ``p5 -> I_t``: the trough fillet's concentric arc ``r_t + depth``
      (on a trough cut; p5 is on it by construction).
    * ``I_t -> I_c``: the flank's parallel offset. ``I_t``/``I_c`` are the
      feet of the tangency normals, where the offset leaves the flank.
    * ``I_c -> p6``: the crest fillet's concentric arc ``r_c - depth``,
      whose top point is p6.

    A corner that crowds its bounding line walks along the flank offset
    (``HANDOVER_*``); the arc it then no longer lies on is replaced by a
    Hermite bridge. Same for a band deeper than the crest is round.

    The edge is sampled at the OUTER edge's fractions of its own total
    arclength -- deliberately one GLOBAL map, not piece by piece. A
    degenerate piece breaks the piecewise map: at a=2.375, b=3.0, R=0.2
    the cut leaves 0.0023 of trough arc on the wall against 0.0524 of its
    offset, and one station would have to jump 3.5% of the inner edge.

    :return: ``(len(outer_fr), 2)`` stations, ``p5`` first, ``p6`` last.
    """
    g = geo["solved"]
    p5 = np.asarray(geo["p5"], dtype=float)
    p6 = np.asarray(geo["p6"], dtype=float)
    Cc, Ct = np.asarray(g["Cc"]), np.asarray(g["Ct"])
    tdir = np.asarray(g["tdir"], dtype=float)     # the flank, DOWNHILL
    has_trough = geo["piece"] == "trough fillet"

    offset = flank(geo).offset(depth)
    I_t, I_c = offset.P0, offset.P1
    clear = HANDOVER_CLEARANCE * depth

    walked_t = False
    if has_trough:
        # I_t must clear the line y = h, or the edge crosses into block2.
        # Not rare: where the cut eats nearly the whole trough fillet, T_t
        # sits barely above the cut and its normal's foot lands BELOW the
        # line (a=1.75, b=3.0, R=0.4: -0.56 depths).
        k = 0
        while k < HANDOVER_STATIONS and I_t[1] - p5[1] < clear:
            k += 1
            I_t = offset.point_at([k / HANDOVER_STATIONS])[0]
            walked_t = True
    k = 0
    walked_c = False
    while k < HANDOVER_STATIONS and I_c[0] < clear:
        k += 1
        I_c = offset.point_at([1.0 - k / HANDOVER_STATIONS])[0]
        walked_c = True

    pieces = []
    if has_trough:
        if walked_t:
            # Leave p5 along the offset arc's tangent (so the seam with
            # block2 is as smooth as the arc's) and arrive along the flank.
            r5 = (p5 - Ct) / np.linalg.norm(p5 - Ct)
            t5 = np.array([-r5[1], r5[0]])
            if float(t5 @ (I_t - p5)) < 0:
                t5 = -t5
            t_end = tdir if float(tdir @ (I_t - p5)) > 0 else -tdir
            pieces.append(_hermite(p5, t5, I_t, t_end))
        else:
            th_5 = math.atan2(p5[1] - Ct[1], p5[0] - Ct[0])
            th_t = math.atan2(I_t[1] - Ct[1], I_t[0] - Ct[0])
            # Sweep the SHORT way: atan2 returns (-pi, pi], and an arc
            # straddling the branch cut otherwise runs the long way round
            # (a=3.0, b=0.5, R=0.3: +309.8 deg for a real -50.2).
            th_t = th_5 + (th_t - th_5 + math.pi) % (2.0 * math.pi) - math.pi
            arc = Arc(Ct, g["r_t"] + depth, th_5, th_t).sample(200)
            arc[0] = p5           # exact at the seam with block2
            pieces.append(arc)
    else:
        I_t = p5                  # the wall starts on the flank

    pieces.append(np.linspace(I_t, I_c, 200))      # exact flank offset

    r_off = g["r_c"] - depth
    if r_off <= 0.0 or walked_c:
        # No concentric offset to use (band deeper than the crest is round,
        # or I_c walked off it): bridge to p6, leaving I_c UPHILL along the
        # flank offset (tdir runs downhill, hence -tdir) and arriving
        # parallel to the apex. Leaving along +tdir doubled the edge back
        # down the flank and folded block3's cells under the crest.
        pieces.append(_hermite(I_c, -tdir, p6, np.array([-1.0, 0.0])))
    else:
        th_c = math.atan2(I_c[1] - Cc[1], I_c[0] - Cc[0])
        apex = math.pi / 2.0
        th_c = apex - ((apex - th_c + math.pi) % (2.0 * math.pi) - math.pi)
        arc = Arc(Cc, r_off, th_c, apex).sample(200)
        arc[-1] = p6              # exact at the axis
        pieces.append(arc)

    path = np.vstack([pieces[0]] + [q[1:] for q in pieces[1:]])
    path[0], path[-1] = p5, p6
    return resample(path, outer_fr)


def wall_points(wall, fractions):
    """Points at the given fractions of the wall's total arclength, placed
    on the exact pieces (a station never cuts a fillet's chord)."""
    L = np.array([p.length for p in wall], dtype=float)
    edges = np.concatenate([[0.0], np.cumsum(L)])
    s = np.asarray(fractions, dtype=float) * edges[-1]
    owner = np.clip(np.searchsorted(edges, s, side="right") - 1,
                    0, len(wall) - 1)
    out = np.empty((len(s), 2))
    for i, piece in enumerate(wall):
        mine = owner == i
        if mine.any():
            u = (s[mine] - edges[i]) / L[i] if L[i] > 0 else np.zeros(mine.sum())
            out[mine] = piece.point_at(np.clip(u, 0.0, 1.0))
    return out


def resample(pts, fractions):
    """Sample the polyline ``pts`` at the given fractions of its arclength."""
    pts = np.asarray(pts, dtype=float)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    s /= s[-1]
    fr = np.asarray(fractions, dtype=float)
    return np.stack([np.interp(fr, s, pts[:, 0]),
                     np.interp(fr, s, pts[:, 1])], axis=1)


def _hermite(P0, m0, P1, m1, n=200):
    """Cubic Hermite from ``P0`` to ``P1`` leaving along ``m0`` and
    arriving along ``m1``, both tangents scaled by the chord."""
    P0, P1 = np.asarray(P0, float), np.asarray(P1, float)
    m0, m1 = np.asarray(m0, float), np.asarray(m1, float)
    chord = float(np.linalg.norm(P1 - P0))
    m0 = m0 / np.linalg.norm(m0) * chord
    m1 = m1 / np.linalg.norm(m1) * chord
    t = np.linspace(0.0, 1.0, n)[:, None]
    return ((2 * t**3 - 3 * t**2 + 1) * P0 + (t**3 - 2 * t**2 + t) * m0
            + (-2 * t**3 + 3 * t**2) * P1 + (t**3 - t**2) * m1)
