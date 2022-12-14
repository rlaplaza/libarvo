"""ARVO python api to compute molecular volumes and surfaces."""
from __future__ import annotations

from ctypes import byref, c_char_p, c_double, c_int, create_string_buffer, POINTER

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
    coordinates: Array2D = np.ascontiguousarray(coordinates, dtype=np.float64)
    radii: Array1D = np.ascontiguousarray(radii, dtype=np.float64)

    if coordinates.ndim != 1:
        if coordinates.shape[0] != radii.shape[0]:
            raise ValueError(
                f"Length of coordinates, {len(coordinates)}, "
                f"and radii, {len(radii)}, must be the same."
            )
    else:
        if coordinates.shape[0] != 3:
            raise ValueError(
                "For a single atom, three (x,y,z) coordinates are required."
            )
        if radii.shape[0] != 1:
            raise ValueError("For a single atom, a single radius is required.")
        coordinates = np.expand_dims(coordinates, axis=0)

    volume_ = c_double()
    area_ = c_double()
    n_v: Array1D = np.zeros(coordinates.shape[0], dtype=np.float64)
    n_s: Array1D = np.zeros(coordinates.shape[0], dtype=np.float64)
    stat_ = c_int()
    errmsg_ = create_string_buffer(100)
    probe_radius = c_double(probe_radius)
    n_atoms = c_int(coordinates.shape[0])

    # Run libarvo
    lib.arvo(
        n_atoms,
        coordinates,
        radii,
        probe_radius,
        byref(volume_),
        byref(area_),
        n_v,
        n_s,
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
    V = float(volume_.value)
    S = float(area_.value)

    return V, S


def atomic_vs(
    coordinates: ArrayLike2D,
    radii: ArrayLike1D,
    probe_radius: float = 0,
) -> tuple[Array1D, Array1D]:
    """Calculates atomic volumes and surfaces.
    Args:
        coordinates: Coordinates as n x 3 matrix (Å)
        radii: vdW radii (Å)
        probe_radius (optional): radius of the probe for SASA (Å)
    Returns:
        n_v: Array of atomic volumes (Å^3)
        n_s: Array of atomic surfaces (Å^2)
    Raises:
        ValueError: If libarvo exits with non-zero error code or if input validation fails.
    """
    # Set up input
    coordinates: Array2D = np.ascontiguousarray(coordinates, dtype=np.float64)
    radii: Array1D = np.ascontiguousarray(radii, dtype=np.float64)

    if coordinates.ndim != 1:
        if coordinates.shape[0] != radii.shape[0]:
            raise ValueError(
                f"Length of coordinates, {len(coordinates)}, "
                f"and radii, {len(radii)}, must be the same."
            )
    else:
        if coordinates.shape[0] != 3:
            raise ValueError(
                "For a single atom, three (x,y,z) coordinates are required."
            )
        if radii.shape[0] != 1:
            raise ValueError("For a single atom, a single radius is required.")
        coordinates = np.expand_dims(coordinates, axis=0)

    volume_ = c_double()
    area_ = c_double()
    n_v: Array1D = np.zeros(coordinates.shape[0], dtype=np.float64)
    n_s: Array1D = np.zeros(coordinates.shape[0], dtype=np.float64)
    stat_ = c_int()
    errmsg_ = create_string_buffer(100)
    probe_radius = c_double(probe_radius)
    n_atoms = c_int(coordinates.shape[0])

    # Run libarvo
    lib.arvo(
        n_atoms,
        coordinates,
        radii,
        probe_radius,
        byref(volume_),
        byref(area_),
        n_v,
        n_s,
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
    V = float(volume_.value)
    S = float(area_.value)

    return n_v, n_s
