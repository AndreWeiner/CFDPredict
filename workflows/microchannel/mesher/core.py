"""The core of the vertical-cut mesh (block4, ``p4 p5 p6``).

A structured block of horizontal rows -- from the axis to each point of
the band's inner edge, ``p5``'s row up to the crest arc's start -- and a
small gmsh cap under the crest fillet, where the rows would narrow to
nothing. A shape with no room for the block is filled by gmsh whole.
Both reuse the neighbouring blocks' own point sequences, so the core is
conformal with them by construction.
"""
import math

import numpy as np

from cut_geometry import flank
from distribution import MIN_CELL_COUNT

# The block's cell count ACROSS its rows, n_x, from the width W, height H
# and area A of the region it covers:
#     c = A / max(W, H) / (CORE_N_SHORT * refinement)
#     n_x = ceil(W / (c * clamp(W/H, 1, CORE_ASPECT_MAX)))
# i.e. c resolves the mean gap across the shorter side, and each direction
# may be coarser than c by its own length ratio, up to CORE_ASPECT_MAX.
# First values from the parametrization repo's core_refinement_test
# prototype (study/03); the
# tuning across the box is still open.
CORE_N_SHORT = 6
CORE_ASPECT_MAX = 4.0
# Upper bound on the block's row count: rows no flatter than
# 1/CORE_ROW_ASPECT_MAX of the widest row's cell width. Separate from
# CORE_ASPECT_MAX so it moves only the row cap, never n_x.
CORE_ROW_ASPECT_MAX = 6.0
# The block stops below the crest so the gmsh cap above it is at least this
# many across-cells c tall.
CORE_CAP_MIN_CELLS = 2.0
# Cell size of the whole-core gmsh fill, as a ratio to the target cell.
CORE_CELL_SIZE_RATIO = 1.0

# gmsh's log messages, collected off the terminal by every gmsh_fill call.
# The caller drains them with take_gmsh_log() (build_case writes log.gmsh).
_GMSH_LOG = []


def take_gmsh_log():
    """The gmsh messages collected since the last call; clears the buffer."""
    lines = list(_GMSH_LOG)
    _GMSH_LOG.clear()
    return lines


def plan_core(geo, refinement):
    """The structured core's ``n_x``, its across-cell size ``c``, and the
    bounds on its row count.

    Sized from the region the block actually covers: ``p5``'s row up to
    the band's offset of the crest arc's start ``T_c``. (Sizing it off the
    flank's own offset band measured the wrong region: as the flank
    shrinks that band's height goes to zero while the block's does not,
    and n_x grew as 1/flank -- 406 at a 1 um flank, where 8 is needed.)

    :return: ``dict(n_x, c, n_y_min, n_y_max)``, or ``None`` when there is
        no wall below the crest to build rows along, or no region.
    """
    if len(geo["wall"]) < 2:
        return None
    top = flank(geo).offset(geo["h"]).P1
    y_b, W_b = float(geo["p5"][1]), float(geo["p5"][0])
    y_t, W_t = min(float(top[1]), geo["p6"][1]), max(0.0, float(top[0]))
    H = y_t - y_b
    if H <= 0.0 or W_b <= 0.0:
        return None
    area = 0.5 * (W_b + W_t) * H
    c = area / max(W_b, H) / (CORE_N_SHORT * refinement)

    def stretch(ratio):
        return min(max(ratio, 1.0), CORE_ASPECT_MAX)
    n_x = max(MIN_CELL_COUNT, math.ceil(W_b / (c * stretch(W_b / H))))
    n_y = math.ceil(H / (c * stretch(H / W_b)))
    # Row bounds -- the wall's own arclength rule is kept between them:
    # at least n_y (a small shape's short wall otherwise gets ~6 rows), at
    # most rows no flatter than 1/CORE_ROW_ASPECT_MAX of the widest row's
    # cell (a shallow flank's long wall otherwise stacks ~34 thin rows).
    n_y_max = math.floor(H * CORE_ROW_ASPECT_MAX * n_x / W_b)
    return dict(n_x=n_x, c=c,
                n_y_min=max(MIN_CELL_COUNT, n_y),
                n_y_max=max(MIN_CELL_COUNT, n_y, n_y_max))


def build_structured_core(bottom, inner, bottom_fr, c, j_top):
    """block4 as a structured block of horizontal rows plus a gmsh cap.

    Row ``j`` runs from the axis ``(0, y_j)`` to the inner-edge point
    ``inner[j]``. Row 0 is block1's upper edge ``bottom`` (``p4 -> p5``),
    laid out by ``bottom_fr``; the top row is uniform; between them each
    row blends the two by height -- what TFI with these boundaries gives.
    The block ends at row ``j_top`` (the crest tangency), pulled down until
    the cap above it is at least ``CORE_CAP_MIN_CELLS * c`` tall; the cap
    is filled by gmsh with every boundary segment pinned.

    :return: the core patch in ``Connect``'s unstructured form, or ``None``
        when there is no room for a block of at least two rows.
    """
    E = np.asarray(inner, dtype=float)               # p5 .. p6
    y6 = E[-1, 1]
    jt = int(j_top)
    while jt > 0 and y6 - E[jt, 1] < CORE_CAP_MIN_CELLS * c:
        jt -= 1
    if jt < 2 or np.any(np.diff(E[:jt + 1, 1]) <= 0.0):
        return None

    fr_b = np.asarray(bottom_fr, dtype=float)
    n_x = len(fr_b) - 1
    fr_t = np.linspace(0.0, 1.0, n_x + 1)
    y0 = E[0, 1]

    def row(j):
        eta = (E[j, 1] - y0) / (E[jt, 1] - y0)
        fr = (1.0 - eta) * fr_b + eta * fr_t
        return np.column_stack([fr * E[j, 0], np.full(n_x + 1, E[j, 1])])
    rows = [row(j) for j in range(jt + 1)]
    if (len(bottom) != n_x + 1
            or not np.allclose(np.asarray(bottom, dtype=float), rows[0],
                               atol=1e-12)):
        raise ValueError("block1's upper edge p4 -> p5 is not the core "
                         "block's first row")

    P = np.vstack(rows)
    quads = [[r * (n_x + 1) + i, r * (n_x + 1) + i + 1,
              (r + 1) * (n_x + 1) + i + 1, (r + 1) * (n_x + 1) + i]
             for r in range(jt) for i in range(n_x)]

    # the cap: top row (axis -> wall), the inner edge on to p6, then down
    # the axis at the top row's own cell width; sizes grown from the loop
    y_t, dx_t = E[jt, 1], E[jt, 0] / n_x
    n_ax = max(1, int(round((y6 - y_t) / dx_t)))
    axis = [(0.0, float(y)) for y in np.linspace(y6, y_t, n_ax + 1)[1:-1]]
    loop = as_tuples(rows[-1][:-1]) + as_tuples(E[jt:]) + axis
    Lp = np.asarray(loop, dtype=float)
    sizes = 0.5 * (np.linalg.norm(np.roll(Lp, -1, 0) - Lp, axis=1)
                   + np.linalg.norm(Lp - np.roll(Lp, 1, 0), axis=1))
    cap = gmsh_fill(loop, sizes)

    # one patch: the block, then the cap (Connect merges the shared row)
    off = len(P)
    P_all = np.vstack([P, np.asarray(cap["P"], dtype=float)])
    return {
        "P": as_tuples(P_all),
        "tri_connectivity": [[v + off for v in t]
                             for t in cap["tri_connectivity"]],
        "quad_connectivity": (clockwise(P_all, quads)
                              + [[v + off for v in q]
                                 for q in cap["quad_connectivity"]]),
    }


def gmsh_fill(boundary_pts, point_sizes):
    """Unstructured quad-dominant fill of one closed loop.

    ``boundary_pts`` reuses the structured blocks' own point sequences
    verbatim; ``point_sizes`` (one per loop point) is grown in from the
    boundary.

    Two guards: the loop must not self-intersect (gmsh does not reject a
    bad loop, it spins on it), and every boundary segment is pinned with
    ``setTransfiniteCurve(tag, 2)`` (without it gmsh re-discretises the
    boundary against its own size field and opens a hole against the
    structured blocks that every topological check still passes).
    """
    import gmsh  # imported lazily: only the core needs it installed

    crossings = _count_crossings(np.asarray(boundary_pts, dtype=float))
    if crossings:
        raise ValueError(
            f"the core's boundary loop self-intersects ({crossings} "
            f"crossing{'s' if crossings > 1 else ''}); gmsh would hang "
            "on it. Usually the inner band edge has folded -- check the "
            "band depth h against the domain's own size.")

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)  # log goes to _GMSH_LOG
        gmsh.logger.start()
        gmsh.model.add("thermalShell_cut_core")
        occ = gmsh.model.occ
        pt_tags = [occ.addPoint(float(x), float(y), 0.0, float(sz))
                   for (x, y), sz in zip(boundary_pts, point_sizes)]
        n = len(pt_tags)
        line_tags = [occ.addLine(pt_tags[i], pt_tags[(i + 1) % n])
                     for i in range(n)]
        occ.addPlaneSurface([occ.addCurveLoop(line_tags)])
        occ.synchronize()
        for tag in line_tags:
            gmsh.model.mesh.setTransfiniteCurve(tag, 2)

        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
        gmsh.model.mesh.generate(2)

        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        index = {int(t): i for i, t in enumerate(node_tags)}
        P = np.asarray(node_coords, dtype=float).reshape(-1, 3)[:, :2]
        cells = {2: [], 3: []}                       # gmsh tri / quad
        for etype, nodes in zip(*gmsh.model.mesh.getElements(dim=2)[::2]):
            if etype in cells:
                k = 3 if etype == 2 else 4
                cells[etype].extend([index[int(v)] for v in c]
                                    for c in np.asarray(nodes).reshape(-1, k))
    finally:
        _GMSH_LOG.extend(gmsh.logger.get())
        gmsh.logger.stop()
        gmsh.finalize()

    return {
        "P": as_tuples(P),
        "tri_connectivity": clockwise(P, cells[2]),
        "quad_connectivity": clockwise(P, cells[3]),
    }


def clockwise(P, cells):
    """``cells`` re-ordered clockwise. gmsh emits counter-clockwise cells,
    the structured blocks come out clockwise (BlockMesh's U-line order);
    mixed winding makes extrudeMesh produce open 3D cells even though
    every 2D cell is individually valid."""
    P = np.asarray(P, dtype=float)
    out = []
    for cell in cells:
        x, y = P[list(cell), 0], P[list(cell), 1]
        twice_area = float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
        out.append(list(reversed(cell)) if twice_area > 0 else list(cell))
    return out


def as_tuples(arr):
    """``(N, 2)`` array -> the list of ``(x, y)`` tuples BlockMesh wants."""
    return [(float(x), float(y)) for x, y in np.asarray(arr, dtype=float)]


def _count_crossings(P):
    """Proper crossings between non-adjacent segments of the closed loop."""
    def orient(a, b, c):
        v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        return 0 if abs(v) < 1e-15 else (1 if v > 0 else -1)

    n, count = len(P), 0
    for i in range(n):
        a, b = P[i], P[(i + 1) % n]
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            c, d = P[j], P[(j + 1) % n]
            o = (orient(a, b, c), orient(a, b, d),
                 orient(c, d, a), orient(c, d, b))
            if 0 not in o and o[0] != o[1] and o[2] != o[3]:
                count += 1
    return count
