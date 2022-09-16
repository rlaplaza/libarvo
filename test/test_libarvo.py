"""Test libarvo."""

import faulthandler

import numpy as np
from numpy.testing import assert_allclose
from numpy.typing import NDArray
import pytest

from libarvo import molecular_vs


def test_sphere():
    coords = np.asarray([0, 0, 0], dtype=float)
    radii = np.asarray([1.7], dtype=float)
    radius_probe = 0.0
    n_atoms = 1
    faulthandler.enable()
    v, s = molecular_vs(coords, radii, radius_probe)
    assert_allclose(v, 20.5795, rtol=1e-5)
    assert_allclose(s, 36.3168, rtol=1e-5)


if __name__ == "__main__":
    test_sphere()
