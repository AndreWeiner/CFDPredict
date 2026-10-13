"""Cell counts and point distributions for the vertical-cut mesher.

Three sizing rules, each matched to what the edge it sizes actually is:

* **across the band** -- a layer count set at refinement 1 from the band's
  depth ``h`` (:func:`band_layers`), scaled linearly by ``refinement``,
  graded fine at the wall. ``h`` is set by the cut, not by the domain, so
  the per-mm budget cannot size it: ``h * NY_PER_MM`` is under one cell at
  the nominal. One array (:attr:`Sizing.band_fr`) serves every edge that
  spans the depth.
* **along the floor** -- see ``CutMesh.build_cut_mesh``.
* **along the wall** -- the hypotenuse's own arclength against the global
  target cell (:func:`wall_fractions`), split piece by piece: each fillet
  by the angle it sweeps, the flank as the slack between them. Sizing the
  wall off the height's budget instead under-resolves it by the slant
  factor (1.55x at the nominal) and folds TFI cells at the crest.
"""
import math
from dataclasses import dataclass

import numpy as np

from LineDistribution import LineDistribution

# Per-mm cell densities at refinement 1: the refinement knob's basis.
NY_PER_MM = 12.0
NX_PER_MM = 12.0

# Hard floor on every block-edge cell count -- the TFI structural minimum.
MIN_CELL_COUNT = 3

# Across-band grading: TOTAL last/first cell ratio over the band's depth
# (blockMesh's simpleGrading convention), fine at the wall. Also used,
# reversed, on block2's corner edges. A steeper ratio on top of a stack
# sized by depth compounds: at 3.0 the floor ran from 0.56x to 1.67x the
# target cell; at 2.0 it is 0.71x to 1.41x.
BAND_EXPANSION = 2.0

# Layers across the band at refinement 1, for a band one target cell
# deep. The count scales with the band's depth as
# (h / ref_len_1) ** BAND_LAYER_EXPONENT -- sub-linear on purpose: the
# stack is graded, so its first cell is already a fraction of h, and sqrt
# lets a deep band coarsen a little instead of piling on layers. A flat
# count shredded a shallow band: at a=3.0, b=0.5 (0.30 target cells deep)
# 4 layers put the first one at 0.040 of a target cell.
BAND_LAYERS_PER_REF = 4
BAND_LAYER_EXPONENT = 0.5

# Floor on the layer count AT REFINEMENT 1, applied before scaling: the
# sqrt rule alone left the shallowest bands (h = 25 um) at 2 layers, too
# few to be a near-wall stack.
BAND_MIN_LAYERS_REF1 = 3

# Curvature weighting on the wall: the arcs claim a share of the cell
# budget larger than their share of the length,
#     share = (L_arc/L) * (1 + g) / (1 + g * L_arc/L)
# so g = 0 is purely length-proportional. At g = 1 the nominal's arcs get
# 28% of the cells for 15% of the length (crest 7.9 deg/cell); g = 3 gave
# a 5.3 deg/cell crest, 2.8x finer than the target cell. SAG_TOL below
# still guarantees the chord error, so coarsening here cannot starve an arc.
CURV_GAIN = 1.0

# Arc-chord SAGITTA tolerance, as a fraction of the channel width b: an
# arc drawn with n cells misses the true curve by R*(1 - cos(theta/2n)),
# which must stay under SAG_TOL * b. An absolute floor under the
# scale-free share above, which alone can starve a short arc.
SAG_TOL = 0.002

# A trough-fillet remnant shorter than this many times the band's last
# cell at p2 (where the wall's distribution starts) is merged into the
# flank. The cut can land almost tangent to the trough fillet and leave a
# sliver -- 0.0055 mm at a=1.125, b=2.375, R=0.5 -- whose one cell would
# set the flank's starting size and crowd ~11 stations into its first
# 0.1 mm. See merge_trough_sliver.
TROUGH_MERGE_FACTOR = 1.0

# A flank shorter than this many times its neighbouring arcs' cell gets NO
# cells: half its length goes to each arc, so the fillets meet in one wall
# station. The gates let the flank shrink to zero (touching fillets), and a
# flank that short is no feature -- given its own cells it only dictates
# sliver-sized cells to the arcs and core rows beside it.
FLANK_MERGE_FACTOR = 1.0

# A flank that keeps its cells gets enough of them to grade from its
# neighbouring arcs' cell sizes with no two adjacent wall cells -- the
# junction pairs included -- differing by more than this ratio. Without it
# a mid-length flank between two fine arcs got 3 cells and a 4.6x jump
# (49 / 199 / 44 um at a=1.125, b=0.5, R=0.3).
FLANK_MAX_RATIO = 1.3


def graded_fractions(n, length):
    """``n + 1`` fractions of an edge, graded by ``BAND_EXPANSION`` and fine
    at fraction 0."""
    return LineDistribution.multi_grading_fractions(
        [(1.0, 1.0, BAND_EXPANSION)], n, length)


def band_layers(h, refinement):
    """Number of layers across a band of depth ``h``::

        n_1 = max(BAND_MIN_LAYERS_REF1,
                  BAND_LAYERS_PER_REF * (h / ref_len_1) ** BAND_LAYER_EXPONENT)
        n   = max(2, round(refinement * n_1))

    Set at refinement 1 (``ref_len_1`` is the target cell there) and then
    scaled LINEARLY by the knob. Taking the depth ratio against the current
    target cell instead -- which already shrinks as 1/refinement -- made
    the count grow as refinement**1.5 and pinned the coarse levels of a
    ladder at the floor (2 / 2 / 2 / 3 layers over r = 0.5 / 0.7 / 1 / 1.3
    at h = 25 um; 2 / 2 / 3 / 4 now). 2 is the absolute floor: a band
    needs an inner and an outer layer to be one.
    """
    ref_len_1 = 2.0 / (NY_PER_MM + NX_PER_MM)
    n_1 = max(BAND_MIN_LAYERS_REF1,
              BAND_LAYERS_PER_REF
              * (float(h) / ref_len_1) ** BAND_LAYER_EXPONENT)
    return max(2, int(round(float(refinement) * n_1)))


@dataclass(frozen=True)
class Sizing:
    """Everything the mesher sizes from, resolved once per build so every
    consumer agrees on it: the refinement knob, the global target cell
    ``ref_len``, and the across-band distribution of a band ``h`` deep."""

    refinement: float
    h: float
    n_band: int
    band_fr: np.ndarray            # fine at the wall, fraction 0

    @classmethod
    def for_band(cls, refinement, h):
        n_band = band_layers(h, refinement)
        return cls(refinement=float(refinement), h=float(h),
                   n_band=n_band, band_fr=graded_fractions(n_band, h))

    @property
    def nx_per_mm(self):
        return NX_PER_MM * self.refinement

    @property
    def ref_len(self):
        """The global target cell size."""
        return 2.0 / (NY_PER_MM * self.refinement + self.nx_per_mm)

    @property
    def band_last_cell(self):
        """The band's coarsest cell, at its inner end -- on ``p1 -> p2``
        the cell at ``p2``, where the wall's distribution starts."""
        return float(np.diff(self.band_fr)[-1] * self.h)


# --- the wall -------------------------------------------------------------

def merge_trough_sliver(spans, last_cell):
    """Fold a too-short trough-fillet remnant into the flank.

    ``spans`` is the wall as ``dict(length, kappa)`` per piece, in wall
    order. The remnant is an artefact of where the cut fell, not a feature:
    below ``TROUGH_MERGE_FACTOR * last_cell`` its length joins the flank,
    which then starts at a size set by its own budget. Only ever the FIRST
    piece, and only an arc followed by a flank.
    """
    if len(spans) < 2:
        return spans
    first, second = spans[0], spans[1]
    if not (first["kappa"] > 0.0 and second["kappa"] == 0.0):
        return spans
    if first["length"] >= TROUGH_MERGE_FACTOR * last_cell:
        return spans
    merged = dict(second)
    merged["length"] = first["length"] + second["length"]
    return [merged] + list(spans[2:])


def wall_fractions(spans, wall_len, ref_len, width, block_rows=None):
    """Arclength fractions (0..1) of the wall's stations, ``p2 -> p3``.

    The count is the wall's own arclength against the target cell, raised
    so the arcs' sagitta demand cannot eat the flank's share; it is then
    split piece by piece (:func:`_piecewise_fractions`).

    :param width: the channel half-width ``b``; the arcs' sagitta
        tolerance is ``SAG_TOL * width``.
    :param block_rows: ``(lo, hi)`` bounds on the structured core's row
        count, or None when there is no structured core.
    """
    tol_abs = SAG_TOL * float(width)
    n_wall = max(MIN_CELL_COUNT, int(round(wall_len / ref_len)))
    n_wall = required_wall_cells(spans, n_wall, tol_abs)
    return _piecewise_fractions(spans, n_wall, tol_abs, block_rows=block_rows)


def _bulge_sizes(h0, h1, L, n):
    """``n`` cell sizes spanning ``L`` that start at ``h0``, grow
    geometrically to a peak, then shrink geometrically to ``h1``.

    The flank's distribution: it has no feature of its own, so its job is
    to get from the trough arc's cell size to the crest arc's without a
    jump at either junction. A single GP is monotone and cannot; two
    back-to-back legs at ratios ``r`` and ``1/r`` can. ``r`` is bisected
    until the legs sum to ``L``; when even uniform overshoots (a short
    flank between coarse arcs) the sizes are simply scaled to fit.

    :return: ``(n,)`` cell sizes summing to exactly ``L``.
    """
    n = int(n)
    if n <= 0:
        return np.zeros(0)
    if n == 1:
        return np.array([L])

    def sizes_for(r):
        """The n sizes implied by growth ratio r, before any rescale."""
        if abs(r - 1.0) < 1e-12:
            return np.linspace(h0, h1, n)
        # k = how many cells are spent on the rising leg: where
        # h0*r**k == h1*r**(n-1-k).
        k = 0.5 * (n - 1 + math.log(h1 / h0) / math.log(r))
        k = float(np.clip(k, 0.0, n - 1))
        i = np.arange(n)
        up = h0 * r ** np.minimum(i, k)
        down = r ** (-np.maximum(i - k, 0.0))
        return up * down

    # Bisect r so the sizes sum to L. sizes_for is increasing in r.
    lo, hi = 1.0, 1.0
    if sizes_for(1.0).sum() < L:
        hi = 1.0
        for _ in range(60):
            hi *= 1.2
            if sizes_for(hi).sum() >= L:
                break
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if sizes_for(mid).sum() < L:
                lo = mid
            else:
                hi = mid
        r = 0.5 * (lo + hi)
    else:
        r = 1.0   # even uniform overshoots; just rescale below

    sizes = sizes_for(r)
    return sizes * (L / sizes.sum())


def _flank_min_cells(h_lo, h_hi, L, max_ratio=FLANK_MAX_RATIO):
    """The fewest cells for which :func:`_bulge_sizes` grades a flank of
    length ``L`` from the arc cell ``h_lo`` to the arc cell ``h_hi`` with
    no two adjacent cells -- the junction pairs included -- differing by
    more than ``max_ratio``.

    ``h_lo`` / ``h_hi`` is ``None`` when no arc sits on that side (the cut
    is on the flank): that end is then free, sized ``L / n``. Scans ``n``
    up to the count at which the flank's cells are as small as the smaller
    arc cell; should none meet the limit, the count with the smallest
    worst ratio is returned.
    """
    ends = [h for h in (h_lo, h_hi) if h is not None]
    if not ends:
        return 1
    n_cap = max(1, int(math.ceil(L / min(ends))))
    best = (math.inf, n_cap)
    for n in range(1, n_cap + 1):
        a = L / n if h_lo is None else h_lo
        b = L / n if h_hi is None else h_hi
        seq = np.concatenate([[] if h_lo is None else [h_lo],
                              _bulge_sizes(a, b, L, n),
                              [] if h_hi is None else [h_hi]])
        worst = float(np.max(np.maximum(seq[1:] / seq[:-1],
                                        seq[:-1] / seq[1:])))
        if worst <= max_ratio:
            return n
        best = min(best, (worst, n))
    return best[1]


def _sagitta_min_cells(length, kappa, tol_abs, cap=60):
    """Fewest cells an arc needs for its chord error to meet ``tol_abs``."""
    if kappa <= 0.0 or tol_abs <= 0.0:
        return 1
    r = 1.0 / kappa
    theta = kappa * length
    n = 1
    while n < cap and r * (1.0 - math.cos(theta / (2.0 * n))) > tol_abs:
        n += 1
    return n


def required_wall_cells(spans, n_total, tol_abs):
    """Raise the wall's cell budget so the flank isn't starved by the arcs.

    The sagitta floor raises any arc whose chord error exceeds ``tol_abs``
    and takes the cells off the flank, the only donor. On a small-``b`` shape
    the tolerance is tight while the budget is small: at
    ``a=0.5, b=0.5, R=0.1`` the crest went 1 -> 4 of 8 cells and left the
    flank 3 for 81% of the wall, one of them 5.5x the target. So the budget
    GROWS instead, to at least (cells the arcs demand) + (cells the flank's
    own length buys at the target size).

    :return: the wall cell count to distribute, ``>= n_total``.
    """
    n_total = int(n_total)
    lengths = np.array([sp["length"] for sp in spans], dtype=float)
    kappas = np.array([sp["kappa"] for sp in spans], dtype=float)
    is_arc = kappas > 0.0
    if not is_arc.any() or is_arc.all():
        return n_total

    n_arcs = sum(_sagitta_min_cells(lengths[i], kappas[i], tol_abs)
                 for i in np.nonzero(is_arc)[0])
    total_len = float(lengths.sum())
    flank_len = float(lengths[~is_arc].sum())
    n_flank = int(math.ceil(n_total * flank_len / total_len))

    return max(n_total, n_arcs + max(n_flank, int((~is_arc).sum())))


def _piecewise_fractions(spans, n_total, tol_abs, block_rows=None):
    """Arclength fractions for the wall, distributed PIECE BY PIECE.

    * **arcs** share a budget proportional to their own ``kappa * length``
      (the angle they turn through) and are divided UNIFORMLY within
      themselves -- nothing prefers one end of a constant-curvature piece.
    * **the flank** gets the rest, distributed by :func:`_bulge_sizes` from
      the trough arc's cell size up and back down to the crest arc's.

    Then, in order: the sagitta floor raises any arc that needs it (taking
    from the flank); a flank shorter than ``FLANK_MERGE_FACTOR`` of its
    neighbouring arcs' cells is merged into them; a kept flank is raised to
    :func:`_flank_min_cells`; and ``block_rows`` holds the structured
    core's row count (every station below the crest arc) by moving the
    flank's count, never below that minimum.

    :return: the fractions; ``n_total + 1`` of them unless the merge or the
        row bounds moved the flank's count.
    """
    n_total = max(int(n_total), len(spans))
    lengths = np.array([sp["length"] for sp in spans], dtype=float)
    kappas = np.array([sp["kappa"] for sp in spans], dtype=float)
    total_len = float(lengths.sum())

    is_arc = kappas > 0.0
    counts = np.zeros(len(spans), dtype=int)

    if not is_arc.any():
        counts[:] = _split_evenly(n_total, len(spans))
    else:
        turn = kappas * lengths
        arc_len = float(lengths[is_arc].sum())
        share = (arc_len / total_len) * (1.0 + CURV_GAIN) / (
            1.0 + CURV_GAIN * (arc_len / total_len))
        n_arcs = int(round(share * n_total))
        n_arcs = int(np.clip(n_arcs, int(is_arc.sum()),
                             n_total - int((~is_arc).sum())))
        counts[is_arc] = _split_weighted(n_arcs, turn[is_arc])
        counts[~is_arc] = _split_evenly(n_total - n_arcs, int((~is_arc).sum()))

        # --- sagitta floor: raise an arc, taking cells off the flank ------
        for i in np.nonzero(is_arc)[0]:
            extra = (_sagitta_min_cells(lengths[i], kappas[i], tol_abs)
                     - counts[i])
            for j in np.nonzero(~is_arc)[0]:
                if extra <= 0:
                    break
                give = min(extra, max(0, counts[j] - 1))
                counts[j] -= give
                counts[i] += give
                extra -= give

    # --- vanishing flank: the two fillets meet --------------------------
    # Its length goes half-and-half to the arcs, so the junction station
    # lands at its midpoint (still on the true wall).
    lengths = lengths.copy()
    for j in np.nonzero(~is_arc)[0]:
        nb = [k for k in (j - 1, j + 1)
              if 0 <= k < len(spans) and is_arc[k] and counts[k] > 0]
        if nb and lengths[j] < FLANK_MERGE_FACTOR * min(
                lengths[k] / counts[k] for k in nb):
            for k in nb:
                lengths[k] += lengths[j] / len(nb)
            lengths[j], counts[j] = 0.0, 0

    def arc_cell(k):
        """Cell size of the arc at index ``k``, or None if there is none."""
        return (lengths[k] / counts[k]
                if 0 <= k < len(spans) and is_arc[k] and counts[k] > 0
                else None)

    # --- the flank's own minimum: grade smoothly between the arcs -------
    n_min = np.ones(len(spans), dtype=int)
    for j in np.nonzero(~is_arc)[0]:
        if counts[j] > 0:
            n_min[j] = _flank_min_cells(arc_cell(j - 1), arc_cell(j + 1),
                                        lengths[j])
            counts[j] = max(counts[j], n_min[j])

    # --- the structured core's row bounds, absorbed by the flank --------
    if block_rows is not None:
        lo, hi = block_rows
        rows = int(counts[:-1].sum())
        for j in np.nonzero(~is_arc)[0]:
            if counts[j] > 0:
                counts[j] = max(n_min[j],
                                counts[j] + int(np.clip(rows, lo, hi)) - rows)

    # --- each piece's own sizes -----------------------------------------
    sizes = []
    for i in range(len(spans)):
        n_i = int(counts[i])
        if n_i <= 0:
            continue
        if is_arc[i]:
            sizes.append(np.full(n_i, lengths[i] / n_i))
        else:
            fallback = lengths[i] / n_i
            h_lo, h_hi = arc_cell(i - 1), arc_cell(i + 1)
            sizes.append(_bulge_sizes(fallback if h_lo is None else h_lo,
                                      fallback if h_hi is None else h_hi,
                                      lengths[i], n_i))
    cum = np.concatenate([[0.0], np.cumsum(np.concatenate(sizes))])
    return cum / cum[-1]


def gp_fractions(h0, L, n):
    """``n + 1`` fractions (0..1) of an edge of length ``L`` split into ``n``
    cells in a geometric progression whose FIRST cell is ``h0``.

    The count and the first cell are both fixed, so the ratio is bisected
    (``h0 * (r**n - 1) / (r - 1) = L`` has no closed form). ``r < 1`` when
    ``n`` cells of ``h0`` would overshoot. Uniform when ``h0`` alone is not
    shorter than the edge.
    """
    n = int(n)
    if n <= 1 or h0 <= 0.0 or h0 >= L:
        return np.linspace(0.0, 1.0, n + 1)

    def total(r):
        return h0 * n if abs(r - 1.0) < 1e-12 else h0 * (r ** n - 1.0) / (r - 1.0)
    lo, hi = 1e-9, 1.0
    while total(hi) < L:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if total(mid) < L:
            lo = mid
        else:
            hi = mid
    sizes = h0 * (0.5 * (lo + hi)) ** np.arange(n)
    cum = np.concatenate([[0.0], np.cumsum(sizes)])
    return cum / cum[-1]


def _split_evenly(n, k):
    """``k`` integers >= 1 summing to ``n``, as equal as possible."""
    if k <= 0:
        return np.zeros(0, dtype=int)
    n = max(int(n), k)
    base = np.full(k, n // k, dtype=int)
    base[: n % k] += 1
    return base


def _split_weighted(n, weights, floor=1):
    """``len(weights)`` integers >= ``floor`` summing to ``n``, in
    proportion to ``weights`` (largest-remainder, so the total is exact)."""
    k = len(weights)
    n = max(int(n), floor * k)
    w = np.asarray(weights, dtype=float)
    if w.sum() <= 0:
        return _split_evenly(n, k)
    exact = w / w.sum() * n
    out = np.maximum(floor, np.floor(exact).astype(int))
    diff = n - out.sum()
    order = np.argsort(-(exact - np.floor(exact)))
    i = 0
    while diff > 0:
        out[order[i % k]] += 1
        diff -= 1
        i += 1
    while diff < 0:
        j = order[i % k]
        if out[j] > floor:
            out[j] -= 1
            diff += 1
        i += 1
    return out
