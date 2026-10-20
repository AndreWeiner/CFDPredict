import numpy as np


class LineDistribution:
    """Parametric point distributions along straight block edges."""

    @staticmethod
    def multi_grading_fractions(divisions, n_total, length):
        """Parametric fractions (0..1) for a line built from several
        geometric-progression pieces stitched together -- the same
        "multi-grading" convention OpenFOAM blockMesh uses for
        ``simpleGrading`` (e.g. a ``centerGrad`` entry in a blockMeshDict):
        each division is ``(length_fraction, cell_fraction, expansion_ratio)``
        where ``expansion_ratio`` is that piece's last-cell/first-cell size
        ratio, solved the same closed-form way blockMesh solves it
        (``ratio ** (1/(n-1))`` for the piece's own cell count).

        The total cell count ``n_total`` is a hard input, never derived
        from a target cell size -- so this degrades gracefully at small
        ``n_total``
        (unlike a single large-ratio GP, which can produce a near-empty or
        exponentially skewed sequence when forced to very few cells) and
        lets total cell count be chosen independently, e.g. proportional to
        the line's physical length so resolution stays consistent as the
        geometry is resized.

        :param divisions: list of ``(length_fraction, cell_fraction,
            expansion_ratio)`` tuples. ``length_fraction``s must sum to 1.0,
            as must ``cell_fraction``s.
        :type divisions: list[tuple[float, float, float]]
        :param n_total: total cell count across all divisions.
        :type n_total: int
        :param length: physical length of the line being distributed.
        :type length: float
        :return: array of fractions from 0.0 to 1.0.
        :rtype: np.ndarray
        """
        divisions = list(divisions)
        length_fracs = np.array([d[0] for d in divisions], dtype=float)
        cell_fracs = np.array([d[1] for d in divisions], dtype=float)
        ratios = np.array([d[2] for d in divisions], dtype=float)

        if not np.isclose(length_fracs.sum(), 1.0):
            raise ValueError("division length_fractions must sum to 1.0")
        if not np.isclose(cell_fracs.sum(), 1.0):
            raise ValueError("division cell_fractions must sum to 1.0")
        if length <= 0:
            raise ValueError("length must be positive")

        counts = LineDistribution.division_counts(cell_fracs, n_total)

        cumulative = [0.0]
        running = 0.0
        for (length_frac, _, ratio), n in zip(divisions, counts):
            piece_length = length_frac * length
            if n == 1:
                sizes = np.array([piece_length])
            elif np.isclose(ratio, 1.0):
                sizes = np.full(n, piece_length / n)
            else:
                r = ratio ** (1.0 / (n - 1))
                first = piece_length * (r - 1.0) / (r**n - 1.0)
                sizes = first * r ** np.arange(n)
            for s in sizes:
                running += s
                cumulative.append(running)

        cumulative = np.array(cumulative)
        return cumulative / cumulative[-1]

    @staticmethod
    def division_counts(cell_fracs, n_total):
        """Integer cell counts per division, at least 1 each, summing to
        exactly ``n_total``.

        Naive rounding can overshoot/undershoot n_total (e.g. cell_fracs
        (.25,.5,.25) at n_total=3 rounds to (1,2,1)=4); rather than dumping
        the whole remainder onto the last division -- which can drive it
        below 1 and raise even though a valid allocation exists -- the
        remainder is distributed one cell at a time, always taking from /
        giving to the currently-largest division, never below the floor
        of 1.

        :param cell_fracs: each division's share of the cells (sum 1).
        :param n_total: total cell count across all divisions.
        :return: ``np.ndarray`` of ints, one per division.
        """
        cell_fracs = np.asarray(cell_fracs, dtype=float)
        n_total = int(n_total)
        if n_total < len(cell_fracs):
            raise ValueError(
                f"n_total={n_total} too small for {len(cell_fracs)} divisions "
                "(each division needs at least 1 cell)"
            )
        counts = np.maximum(1, np.round(cell_fracs * n_total).astype(int))
        diff = n_total - counts.sum()
        order = np.argsort(-counts)  # largest first
        i = 0
        guard = 0
        while diff != 0 and guard < 10 * len(counts) + 10:
            j = order[i % len(counts)]
            if diff > 0:
                counts[j] += 1
                diff -= 1
            elif counts[j] > 1:
                counts[j] -= 1
                diff += 1
            i += 1
            guard += 1
        if diff != 0:
            raise ValueError(
                f"n_total={n_total} too small for these cell_fractions"
            )
        return counts

    @staticmethod
    def apply_fractions(p_start, p_end, fractions):
        """Map parametric fractions (0..1) onto points along the straight line
        from ``p_start`` to ``p_end``.

        :param p_start: line start point.
        :type p_start: tuple[float, float]
        :param p_end: line end point.
        :type p_end: tuple[float, float]
        :param fractions: monotonic fractions from 0.0 to 1.0.
        :type fractions: array-like
        :return: list of ``(x, y)`` points, one per fraction.
        :rtype: list[tuple[float, float]]
        """
        p_start = np.asarray(p_start, dtype=float)
        p_end = np.asarray(p_end, dtype=float)
        fractions = np.asarray(fractions, dtype=float)
        pts = p_start[None, :] + fractions[:, None] * (p_end - p_start)[None, :]
        return [tuple(pt) for pt in pts]
