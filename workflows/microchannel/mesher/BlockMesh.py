import numpy as np
from scipy.interpolate import splprep

"""
This code is taken and adapted according to requirement from PyAero repository.
GitHub: https://github.com/chiefenne/PyAero
"""


class BlockMesh:
    """Container class for block-structured 2D meshes.

    The mesh is represented as a list of *U-lines* (lines in u-direction), where
    each U-line is a list of `(x, y)` tuples.
    """

    def __init__(self):
        self.ULines = list()

    def getDivUV(self):
        """Return the number of divisions in u and v directions.

        For a structured grid, if there are `U+1` points along each U-line and
        `V+1` U-lines total, then the number of divisions is `(U, V)`.

        :return: `(U, V)` where U is divisions in u-direction and V in v-direction.
        :rtype: tuple[int, int]
        """
        u = len(self.getULines()[0]) - 1
        v = len(self.getULines()) - 1
        return u, v

    def getULines(self):
        """Return the stored U-lines.

        :return: list of U-lines.
        :rtype: list[list[tuple[float, float]]]
        """
        return self.ULines

    def transfinite(self, boundary):
        """Fill the block by transfinite interpolation (TFI) from its four
        boundary curves.

        Reference:
        http://en.wikipedia.org/wiki/Transfinite_interpolation

                       upper
                --------------------
                |                  |
                |                  |
           left |                  | right
                |                  |
                |                  |
                --------------------
                       lower

        :param boundary: `[lower, upper, left, right]`, each a list of
            `(x, y)` tuples. `lower`/`upper` must run in the same direction,
            as must `left`/`right`.
        :type boundary: list[list[tuple[float, float]]]
        """
        lower, upper, left, right = (np.array(b) for b in boundary)

        # convert the block boundary curves into parametric form
        # as curves need to be between 0 and 1
        # interpolate B-spline through data points
        # here, a linear interpolant is derived "k=1"
        # splprep returns:
        # tck ... tuple (t,c,k) containing the vector of knots,
        #         the B-spline coefficients, and the degree of the spline.
        #   u ... array of the parameters for each given point (knot)
        tck_lower, u_lower = splprep(lower.T, s=0, k=1)
        tck_upper, u_upper = splprep(upper.T, s=0, k=1)
        tck_left, u_left = splprep(left.T, s=0, k=1)
        tck_right, u_right = splprep(right.T, s=0, k=1)

        nodes = np.zeros((len(left) * len(lower), 2))

        # corner points
        c1 = lower[0]
        c2 = upper[0]
        c3 = lower[-1]
        c4 = upper[-1]

        for i, xi in enumerate(u_lower):
            for j, eta in enumerate(u_left):

                node = i * len(u_left) + j

                point = (
                    (1.0 - xi) * left[j]
                    + xi * right[j]
                    + (1.0 - eta) * lower[i]
                    + eta * upper[i]
                    - (
                        (1.0 - xi) * (1.0 - eta) * c1
                        + (1.0 - xi) * eta * c2
                        + xi * (1.0 - eta) * c3
                        + xi * eta * c4
                    )
                )

                nodes[node, 0] = point[0]
                nodes[node, 1] = point[1]

        vlines = list()
        vline = list()
        i = 0
        for node in nodes:
            i += 1
            vline.append(node)
            if i % len(left) == 0:
                vlines.append(vline)
                vline = list()

        vlines.reverse()
        self.ULines = self.makeUfromV(vlines)

    @staticmethod
    def makeUfromV(vlines):
        """Convert a list of V-lines into a list of U-lines.

        This is a helper to switch the internal representation after generating
        nodes in V-line order.

        :param vlines: list of V-lines.
        :type vlines: list[list[tuple[float, float]]]
        :return: list of U-lines.
        :rtype: list[list[tuple[float, float]]]
        """
        ulines = list()
        uline = list()
        for i in range(len(vlines[0])):
            for vline in vlines:
                x, y = vline[i][0], vline[i][1]
                uline.append((x, y))
            ulines.append(uline[::-1])
            uline = list()
        return ulines
