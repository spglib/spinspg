"""Independent operation-set oracles for the public input-cell API."""

from itertools import permutations, product

import numpy as np
import pytest

from spinspg import get_spin_symmetry
from spinspg.spin import SpinOnlyGroupType


@pytest.mark.parametrize("backend", ["spglib", "moyopy"])
@pytest.mark.parametrize("kind", list(SpinOnlyGroupType))
def test_fcc_operation_set(fcc, backend, kind):
    lattice, positions, numbers, _ = fcc
    if kind == SpinOnlyGroupType.NONMAGNETIC:
        moments = np.zeros((4, 3))
    elif kind == SpinOnlyGroupType.COLLINEAR:
        moments = np.array([[0, 0, 1], [0, 0, 1], [0, 0, -1], [0, 0, -1]])
    elif kind == SpinOnlyGroupType.COPLANAR:
        moments = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]])
    else:
        moments = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])

    # The cubic point group consists of all signed coordinate permutations.
    # The four FCC sites are also its input-cell centering translations.
    expected = {}
    for axes in permutations(range(3)):
        for signs in product([-1, 1], repeat=3):
            rotation = np.eye(3, dtype=int)[list(axes)] * np.array(signs)[:, None]
            for translation in positions:
                transformed = (positions @ rotation.T + translation) % 1
                site_map = [
                    next(i for i, site in enumerate(positions) if np.array_equal(site, target))
                    for target in transformed
                ]
                targets = moments[site_map]
                # An orthogonal spin map exists exactly when all inner products agree.
                if not np.array_equal(targets @ targets.T, moments @ moments.T):
                    continue
                spin_rotation = np.eye(3)
                if kind == SpinOnlyGroupType.COLLINEAR:
                    spin_rotation[2, 2] = np.dot(targets[:, 2], moments[:, 2]) / 4
                elif kind == SpinOnlyGroupType.COPLANAR:
                    spin_rotation[:2, :2] = targets[:, :2].T @ moments[:, :2] / 2
                    spin_rotation[2, 2] = np.linalg.det(spin_rotation[:2, :2])
                elif kind == SpinOnlyGroupType.NONCOPLANAR:
                    spin_rotation = targets.T @ moments / 4
                key = (tuple(rotation.ravel()), tuple(translation))
                expected[key] = spin_rotation

    sog, rotations, translations, spin_rotations = get_spin_symmetry(
        lattice, positions, numbers, moments, backend=backend
    )
    assert sog.spin_only_group_type == kind
    if sog.axis is not None:
        np.testing.assert_allclose(np.abs(sog.axis), [0, 0, 1])
    assert len(rotations) == len(expected)
    for rotation, translation, spin_rotation in zip(rotations, translations, spin_rotations):
        # All translations here are exact half-integers; rounding only removes
        # backend noise when constructing dictionary keys.
        key = (tuple(rotation.ravel()), tuple(np.round(translation % 1, 8) % 1))
        np.testing.assert_allclose(spin_rotation, expected.pop(key), atol=1e-12)
    assert not expected


@pytest.mark.parametrize("backend", ["spglib", "moyopy"])
@pytest.mark.parametrize("factor,order", [(0.99, 192), (1.01, 32)])
def test_magnetic_tolerance_boundary(fcc, backend, factor, order):
    lattice, positions, numbers, _ = fcc
    tolerance = 1e-3
    moments = np.zeros((4, 3))
    moments[:, 2] = 1
    moments[positions[:, 0] == 0, 2] += factor * tolerance
    sog, rotations, _, _ = get_spin_symmetry(
        lattice, positions, numbers, moments, mag_symprec=tolerance, backend=backend
    )
    assert sog.spin_only_group_type == SpinOnlyGroupType.COLLINEAR
    # Above the threshold only 16 x-axis-preserving rotations and two
    # translations preserving each moment magnitude remain.
    assert len(rotations) == order
