import numpy as np
import pytest
from spglib import get_symmetry_dataset

from spinspg.permutation import Permutation, get_symmetry_permutations, is_overlap_with_origin


def test_symmetry_permutations(fcc):
    lattice, positions, numbers, _ = fcc
    symprec = 1e-5
    dataset = get_symmetry_dataset((lattice, positions, numbers), symprec)
    rotations = dataset.rotations
    translations = dataset.translations
    permutations = get_symmetry_permutations(
        lattice, positions, numbers, rotations, translations, symprec
    )

    assert len(permutations) == len(rotations)
    for permutation in permutations:
        assert np.all(np.sort(permutation.permutation) == np.arange(len(positions)))


def test_permutation_composition():
    p = Permutation(np.array([2, 0, 3, 1]))
    q = Permutation(np.array([0, 2, 1, 3]))
    product = p * q
    assert [product(i) for i in range(4)] == [p(q(i)) for i in range(4)]
    assert not np.array_equal(product.permutation, (q * p).permutation)
    product.permutation[:] = 0
    np.testing.assert_array_equal(p.permutation, [2, 0, 3, 1])


def test_species_matching_and_first_available_site():
    # Both nearby sites of a species are within tolerance. Greedy matching must
    # keep the first unused index, even when a later index is nearer.
    lattice = np.array([[2, 0, 0], [0.3, 1, 0], [0.1, 0.2, 3]])
    positions = np.array([[0.002, 0, 0], [0, 0, 0], [0.001, 0, 0], [0, 0, 0]])
    numbers = np.array([0, 1, 0, 1])
    rotations = np.array([np.eye(3, dtype=int), -np.eye(3, dtype=int)])
    translations = np.array([[1, 0, 0], [0, 0, 0]])
    tolerance = 0.02
    actual = get_symmetry_permutations(
        lattice, positions, numbers, rotations, translations, tolerance
    )
    for rotation, translation, permutation in zip(rotations, translations, actual):
        expected = []
        for site, number in zip(positions, numbers):
            target = rotation @ site + translation
            match = next(
                j
                for j, candidate in enumerate(positions)
                if j not in expected
                and numbers[j] == number
                and is_overlap_with_origin(lattice, target - candidate, tolerance)
            )
            expected.append(match)
        np.testing.assert_array_equal(permutation.permutation, expected)


@pytest.mark.parametrize("distance", [np.nextafter(0.125, 0), 0.125, np.nextafter(0.125, 1)])
def test_site_matching_strict_cartesian_threshold(distance):
    lattice = np.diag([2.0, 3.0, 4.0])
    arguments = (
        lattice,
        np.zeros((1, 3)),
        np.array([0]),
        np.array([np.eye(3, dtype=int)]),
        np.array([[distance / 2, 0, 0]]),
        0.125,
    )
    if distance < 0.125:
        (permutation,) = get_symmetry_permutations(*arguments)
        np.testing.assert_array_equal(permutation.permutation, [0])
    else:
        with pytest.raises(ValueError, match="does not map all sites"):
            get_symmetry_permutations(*arguments)
