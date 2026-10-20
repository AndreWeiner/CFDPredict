import numpy as np
import os
import logging

logger = logging.getLogger(__name__)
"""
This code is taken in part and adapted according to requirement from PyAero repository.
GitHub: https://github.com/chiefenne/PyAero 
"""

class Mesh:
    """Lightweight mesh container with simple export utilities."""

    def __init__(self, vertices, connectivity):
        """Construct a mesh container.

        :param vertices: vertex coordinates with shape ``(N, 3)`` or ``(N, 2)``
        :type vertices: np.ndarray
        :param connectivity: element-to-vertex indices (0-based), e.g. TRI3 as ``(M, 3)``
        :type connectivity: np.ndarray
        """
        self.mesh = vertices, connectivity

    def write_obj(self, filename):
        """Write mesh to a OBJ file.

        Vertices are written as ``v x y z``. If vertices are 2D (shape ``(N, 2)``),
        a zero z-coordinate is appended. Faces are written as ``f i j k ...`` using
        1-based indexing as required by the OBJ specification.

        :param filename: output
        :type filename: str
        :return: None
        :rtype: None
        :raises OSError: if the file cannot be created or written
        """
        with open(filename, "w") as f:
            # Write vertices
            vertices, connectivity = self.mesh
            for v in vertices:
                if len(v) == 2:
                    f.write(f"v {v[0]} {v[1]} 0.0\n")  # Add z=0 for 2D vertices
                else:
                    f.write(f"v {v[0]} {v[1]} {v[2]}\n")

            # Write faces (convert 0-based to 1-based indices)
            for face in connectivity:
                indices_str = " ".join(str(i + 1) for i in face)
                f.write(f"f {indices_str}\n")
