"""Feasibility and construction of the TRUE fluid domain (the physical
foil geometry), and the vertical-cut CFD approximation derived from it.

The foil folds back on itself, so the fluid meets two fillets of unequal
fluid-facing radius -- the crest fillet at R (fluid on the concave side)
and the trough fillet at R+s (fluid on the convex side) -- joined by their
internal common tangent (the flank). That construction (documented in the
parametrization repo's study/05_domain_description/geometry_construction.md)
is what must be possible for the parameters to describe a real channel.

solve_true_geometry() solves that profile; solve_cfd_outline() then applies
the CFD approximation the mesher actually meshes -- a VERTICAL CUT at a
fixed height h, with its corner at (x_h, 0). mesher/CutMesh.py builds the
mesh on it.

The wall is modelled as a list of pieces -- :class:`Arc` and :class:`Line`
-- in the outline's counter-clockwise sense, so the fluid is always on a
piece's LEFT. :func:`true_wall` is the whole profile, floor contact to
apex; ``solve_cfd_outline(...)["wall"]`` is what survives the cut, Q to
apex. The figures (geometry_figures.py) and the mesher's band edges are
built from those pieces rather than from per-consumer branches on which fillet
the cut landed on.

check_parameter_values() (gate 1) and check_true_geometry() (gate 2) state
the feasibility conditions explicitly and return one message per violated
condition; py_microchannel.py runs both before building anything.

All inputs and outputs are in mm; everything here is
scale-invariant, so no unit conversion happens in this module.
"""
import math

import numpy as np


class Arc:
    """Circular arc about ``C`` of radius ``r``, from angle ``t0`` to ``t1``
    (either sense). Travelled counter-clockwise (``t1 > t0``) the fluid is
    inside the circle; clockwise, outside."""

    kappa = property(lambda self: 1.0 / self.r)

    def __init__(self, C, r, t0, t1):
        self.C, self.r = np.asarray(C, dtype=float), float(r)
        self.t0, self.t1 = float(t0), float(t1)

    @property
    def length(self):
        return self.r * abs(self.t1 - self.t0)

    def point_at(self, u):
        """Points at parameter ``u`` in [0, 1] (uniform in angle/arclength)."""
        t = self.t0 + np.asarray(u, dtype=float) * (self.t1 - self.t0)
        return np.column_stack([self.C[0] + self.r * np.cos(t),
                                self.C[1] + self.r * np.sin(t)])

    def sample(self, n):
        t = np.linspace(self.t0, self.t1, n)
        return np.column_stack([self.C[0] + self.r * np.cos(t),
                                self.C[1] + self.r * np.sin(t)])

    def offset(self, d):
        """The concentric arc ``d`` into the fluid (negative: out of it)."""
        ccw = 1.0 if self.t1 > self.t0 else -1.0
        return Arc(self.C, self.r - ccw * d, self.t0, self.t1)

    def x_at_y(self, y):
        """Where the full circle crosses ``y``, on the side ``t0`` is on."""
        dy = y - self.C[1]
        if abs(dy) >= self.r:
            raise ValueError(f"the circle (r={self.r:.6g}) never reaches "
                             f"y={y:.6g}")
        return float(self.C[0] + math.copysign(
            math.sqrt(self.r * self.r - dy * dy), math.cos(self.t0)))


class Line:
    """Straight segment ``P0 -> P1``; the fluid is on its left.

    ``direction`` may be given exactly: the flank's comes from the solved
    tangent, which stays defined when the fillets touch and the flank has
    zero length."""

    kappa = 0.0

    def __init__(self, P0, P1, direction=None):
        self.P0 = np.asarray(P0, dtype=float)
        self.P1 = np.asarray(P1, dtype=float)
        self._direction = (None if direction is None
                           else np.asarray(direction, dtype=float))

    @property
    def length(self):
        return float(np.linalg.norm(self.P1 - self.P0))

    @property
    def direction(self):
        if self._direction is not None:
            return self._direction
        return (self.P1 - self.P0) / self.length

    @property
    def normal(self):
        """Unit normal into the fluid (the left of the direction of travel)."""
        t = self.direction
        return np.array([-t[1], t[0]])

    def point_at(self, u):
        u = np.asarray(u, dtype=float)[:, None]
        return self.P0 + u * (self.P1 - self.P0)

    def sample(self, n):
        return np.linspace(self.P0, self.P1, n)

    def offset(self, d):
        """The parallel segment ``d`` into the fluid (negative: out of it)."""
        n = d * self.normal
        return Line(self.P0 + n, self.P1 + n, self._direction)

    def x_at_y(self, y):
        """Where the infinite line through the segment crosses ``y``."""
        dy = self.P1[1] - self.P0[1]
        if abs(dy) < 1e-12:
            raise ValueError("a horizontal line has no single crossing")
        return float(self.P0[0]
                     + (y - self.P0[1]) / dy * (self.P1[0] - self.P0[0]))


def solve_true_geometry(a, b, s, R):
    """Solve the true fluid-domain construction.

    Returns a dict with the fillet centres/radii, the flank tangent points
    and angles, and the flank length. Raises ValueError when the inner tangent does not exist
    (fillets overlap).
    """
    r_c, r_t = R, R + s          # fluid-facing radii (crest / trough)
    a_f = a - s                  # fluid apex; `a` is the outer envelope
    Cc = np.array([0.0, a_f - r_c])
    Ct = np.array([b, r_t])

    D = Ct - Cc
    L = float(np.linalg.norm(D))
    # Touching fillets (r_c + r_t == L) are allowed: a zero-length flank is
    # manufacturable, and the construction holds -- al = 90 deg, and both
    # tangency points collapse onto the internal homothetic centre P, the
    # point where the fillets touch.
    if r_c + r_t > L:
        raise ValueError(
            f"fillets overlap: r_c+r_t={r_c + r_t:.6g} > L={L:.6g}")

    phi = math.atan2(D[1], D[0])
    al = math.asin(min(1.0, (r_c + r_t) / L))
    tdir = np.array([math.cos(phi - al), math.sin(phi - al)])
    nrm = np.array([-tdir[1], tdir[0]])
    P = Cc + (r_c / (r_c + r_t)) * (Ct - Cc)   # internal homothetic centre
    T_c = Cc - nrm * np.dot(Cc - P, nrm)
    T_t = Ct - nrm * np.dot(Ct - P, nrm)
    flank = math.sqrt(max(0.0, L * L - (r_c + r_t) ** 2))
    a_c = math.atan2(T_c[1] - Cc[1], T_c[0] - Cc[0])

    # a_t: the trough tangency angle, expressed so the trough FILET is the
    # MINOR (clockwise-exterior) arc from the floor contact.
    #
    # The fluid wets the outside of the trough fold. Measured from the
    # floor-contact angle (-90 deg, the bottom of the trough circle), the
    # wetted arc runs CLOCKWISE -- from -90 deg down through the lower-left
    # quadrant -- to wherever the flank meets the fold:
    #     sweep = (-90 - a_t_raw) mod 360
    # and a_t = -90 - sweep, so it always lies in (-450, -90] deg.
    #
    # The raw atan2 angle jumps between branches exactly where the flank
    # turns overhang (case #105-family: the tangency rides just past the
    # fold's left horizontal): atan2 then reports ~+172 deg, and drawing
    # from -90 deg straight to +172 deg sweeps the MAJOR arc the wrong way
    # round (262 deg of fillet instead of the 98 deg the fluid actually
    # wets). Unwrapping to the clockwise-exterior sweep keeps the fillet the
    # minor arc at every tangency. For the gate-passing family (a_t_raw in
    # (-180, -90) deg, sweep <= 90 deg) this is the identity.
    a_t_raw = math.atan2(T_t[1] - Ct[1], T_t[0] - Ct[0])
    sweep_t = (-math.pi / 2 - a_t_raw) % (2 * math.pi)
    a_t = -math.pi / 2 - sweep_t

    return dict(r_c=r_c, r_t=r_t, a_f=a_f, Cc=Cc, Ct=Ct, T_c=T_c, T_t=T_t,
                tdir=tdir,
                a_c=a_c,
                a_t=a_t,
                sweep_t=sweep_t,
                flank=flank, L=L)


def flank_line(g, start):
    """The flank from ``start`` (T_t, or the cut point Q on it) up to T_c."""
    return Line(start, g["T_c"], direction=-np.asarray(g["tdir"]))


def true_wall(g):
    """The whole true profile, floor contact ``(b, 0)`` to the apex: the
    trough fillet (clockwise -- the fluid wets its outside), the flank, the
    crest fillet (counter-clockwise)."""
    return [Arc(g["Ct"], g["r_t"], -math.pi / 2, g["a_t"]),
            flank_line(g, g["T_t"]),
            Arc(g["Cc"], g["r_c"], g["a_c"], math.pi / 2)]


def solve_cfd_outline(a, b, s, R, h=0.05):
    """The CFD approximation: a vertical cut at fixed height h.

    Mirrors the construction documented in the parametrization repo's
    study/05_domain_description/h_profile.py: the horizontal y = h meets
    the true hypotenuse (trough fillet -> slant -> crest fillet) at Q;
    dropping straight down from Q gives the corner (x_h, 0), and everything
    right of x = x_h is handed to the solid.

    :param h: cut height above the floor, mm; must lie in (0, a_f).
    :return: :func:`solve_true_geometry`'s dict extended with ``h``, ``Q``,
        ``x_h``, ``piece`` (which hypotenuse branch owns the cut), ``th_h``
        -- the polar angle of Q about the owning fillet's centre (None on
        the slant branch), already mapped into (-pi, pi] so arc bookkeeping
        stays in one branch -- and ``wall``, the pieces that survive the
        cut, Q to the apex. Raises ValueError when the true geometry does
        not exist or h is out of range.
    """
    g = solve_true_geometry(a, b, s, R)
    r_c, r_t = g["r_c"], g["r_t"]
    Cc, Ct, T_c, T_t = g["Cc"], g["Ct"], g["T_c"], g["T_t"]
    crest = Arc(Cc, r_c, g["a_c"], math.pi / 2)

    if h <= 0.0 or h >= g["a_f"]:
        raise ValueError(f"cut height h={h} must lie in (0, a_f={g['a_f']})")

    if h <= T_t[1]:        # on the trough fillet
        # study's th_h = pi - asin(...) sits in (pi, 3pi/2); map into
        # (-pi, -pi/2) so it lives on the same branch as a_t.
        th_raw = math.pi - math.asin((h - r_t) / r_t)
        th_h = th_raw - 2.0 * math.pi
        Q = np.array([Ct[0] + r_t * math.cos(th_raw), h])
        piece = "trough fillet"
        wall = [Arc(Ct, r_t, th_h, g["a_t"]), flank_line(g, T_t), crest]
    elif h <= T_c[1]:      # on the slant flank
        th_h = None
        Q = T_t + (h - T_t[1]) / (T_c[1] - T_t[1]) * (T_c - T_t)
        piece = "slant flank"
        wall = [flank_line(g, Q), crest]
    else:                  # on the crest fillet
        th_h = math.asin((h - Cc[1]) / r_c)
        Q = np.array([Cc[0] + r_c * math.cos(th_h), h])
        piece = "crest fillet"
        wall = [Arc(Cc, r_c, th_h, math.pi / 2)]

    g.update(h=h, Q=Q, x_h=float(Q[0]), piece=piece, th_h=th_h, wall=wall)
    return g


# Every build parameter: a zero or negative value is meaningless for each
# of them (a,b,s,R shape the cross-section; l the extrusion length;
# N_rot the oscillation period T = 60/N_rot) and a non-number would crash
# the build with a less readable error downstream.
BUILD_PARAMETERS = ("l", "a", "s", "b", "R", "N_rot")


def check_parameter_values(shape):
    """One message per build parameter that is missing, non-numeric,
    non-finite, zero or negative. Empty means all six are usable.

    Complements check_true_geometry: this validates the parameter *values*
    themselves; the geometry gate then decides whether the values they
    form a possible channel with."""
    violations = []
    for name in BUILD_PARAMETERS:
        v = shape.get(name) if shape is not None else None
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            violations.append(f"{name} = {v!r} is missing or not a number")
        elif not math.isfinite(v):
            violations.append(f"{name} = {v!r} is not finite")
        elif v <= 0:
            violations.append(f"{name} = {v!r} must be positive")
    return violations


def check_true_geometry(a, b, s, R):
    """One message per way the true fluid domain fails to exist.

    Conditions, in the order checked:
      T1  the two fluid-facing fillets must not overlap, i.e. an internal
          common tangent (the flank) must exist:  2R + s <= L
          (equality: the fillets touch and the flank has zero length)
      T3  the flank's crest tangency must sit on the crest arc's right
          half (-90 < a_c < 90 deg), or the crest fillet folds over itself
      T4  the flank's trough tangency must sit in the trough arc's upper
          left quadrant (-180 < a_t < -90 deg), or the flank closes onto
          the trough fillet from the wrong side
      T5  the whole flank must stay above the floor y = 0

    T2 (minimum flank length) was removed 2026-09-21: a vanishing flank --
    the two-fillets-touching limit -- is an acceptable construction, down
    to and including zero length (2R + s == L, the fillets touching), which
    is manufacturable. T6 (the flank extended to the floor must land in
    0 < p2 < b) was removed 2026-10-06: it guarded the old "flank to the
    floor" outline, which the vertical-cut mesher does not use, and over
    200k random shapes in the box it never fired without T4. The remaining
    gate IDs are unchanged so the README's references keep their
    meaning.

    Expects values that passed check_parameter_values (gate 1): finite
    and positive. Returns a list of human-readable violations; empty means
    feasible."""
    violations = []
    try:
        g = solve_true_geometry(a, b, s, R)
    except ValueError as exc:
        return [f"T1 {exc.args[0]}: no straight flank exists"]

    a_c, a_t = math.degrees(g["a_c"]), math.degrees(g["a_t"])
    T_c, T_t = g["T_c"], g["T_t"]

    if not -90.0 < a_c < 90.0:
        violations.append(
            f"T3 crest tangency folds over: flank meets the crest fillet at "
            f"a_c = {a_c:.2f} deg (needs -90 < a_c < 90)")
    if not -180.0 < a_t < -90.0:
        violations.append(
            f"T4 trough tangency wrong side: flank meets the trough fillet at "
            f"a_t = {a_t:.2f} deg (needs -180 < a_t < -90, i.e. the fluid's "
            f"wetted trough arc must stay under 90 deg; here it is "
            f"{math.degrees(g['sweep_t']):.1f} deg); the flank cannot close "
            f"onto the trough fillet without over-hanging it")
    if T_c[1] <= 0 or T_t[1] <= 0:
        violations.append(
            f"T5 flank crosses the floor: tangent points at y = "
            f"{T_c[1]:.4f} (crest) / {T_t[1]:.4f} (trough)")
    return violations
