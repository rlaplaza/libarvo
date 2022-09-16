"""Test libarvo."""

import numpy as np
from numpy.testing import assert_allclose
from numpy.typing import NDArray
import pytest

from libarvo import molecular_vs


def test_onesphere():
    coords = np.asarray([0, 0, 0], dtype=float)
    radii = np.asarray([1.7], dtype=float)
    radius_probe = 0.0
    n_atoms = 1
    v, s = molecular_vs(coords, radii, radius_probe)
    assert_allclose(v, 20.5795, rtol=1e-5)
    assert_allclose(s, 36.3168, rtol=1e-5)


def test_threespheres():
    coords = np.asarray([[0, 0, 0], [0, 0, 3.4], [3.4, 0, 0]], dtype=float)
    radii = np.asarray([1.7, 1.7, 1.7], dtype=float)
    radius_probe = 0.0
    n_atoms = 3
    v, s = molecular_vs(coords, radii, radius_probe)
    assert_allclose(v, 61.7386, rtol=1e-5)


def test_twoospheres():
    rs = [1.7, 1.7, 1.7]
    ds = [3.4, 2.8, 1.7]
    for r, d in zip(rs, ds):
        coords = np.asarray([[0, 0, -d / 2], [0, 0, d / 2]], dtype=float)
        radii = np.asarray([r, r], dtype=float)
        radius_probe = 0.0
        n_atoms = 2
        v, s = molecular_vs(coords, radii, radius_probe)
        # Analytical expression for the overlap area of two spheres
        av = 2 * 20.5795 - (1 / 12) * np.pi * (4 * r + d) * (2 * r - d) ** 2
        assert_allclose(v, av, rtol=1e-5)


if __name__ == "__main__":
    test_onesphere()
    test_threespheres()
    test_twoospheres()
