"""ARVO python api to compute molecular volumes and surfaces."""
from __future__ import annotations

from ctypes import byref, c_double, c_int, create_string_buffer

import numpy as np

from libarvo.lib import lib
from libarvo.mytyping import Array1D, Array2D, ArrayLike1D, ArrayLike2D


def molecular_vs(
    coordinates: ArrayLike2D,
    radii: ArrayLike1D,
    probe_radius: float = 0,
) -> tuple[float, float]:
    """Calculates molecular volume and surface.
    Args:
        coordinates: Coordinates as n x 3 matrix (Å)
        radii: vdW radii (Å)
        probe_radius (optional): radius of the probe for SASA (Å)
    Returns:
        V: Molecular volume (Å^3)
        S: Molecular surface (Å^2)
    Raises:
        ValueError: If libarvo exits with non-zero error code or if input validation fails.
    """
    # Set up input
    coordinates: Array2D = np.ascontiguousarray(coordinates)
    radii: Array1D = np.ascontiguousarray(radii, dtype=np.float64)

    if len(coordinates) != len(radii):
        raise ValueError(
            f"Length of coordinates, {len(coordinates)}, "
            f"and radii, {len(radii)}, must be the same."
        )

    V_ = c_double()
    S_ = c_double()
    stat_ = c_int()
    errmsg_ = create_string_buffer(100)
    n_atoms = c_int(coordinates.shape[0])

    # Run libconeangle
    lib.arvo(
        n_atoms,
        coordinates,
        radii,
        V,
        S,
        byref(stat_),
        errmsg_,
    )

    # Check for non-zero exit code
    stat = int(stat_.value)
    if stat != 0:
        raise ValueError(
            f"libarvo exited with non-zero exit code: {stat}. "
            f"Error message: {errmsg_.value.decode()}"
        )

    # Take out results
    V = float(V_.value)
    S = float(S_.value)

    return V, S
