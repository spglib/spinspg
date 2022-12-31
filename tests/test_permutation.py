import numpy as np
from spglib import get_symmetry_dataset

from spinspg.permutation import get_symmetry_permutations


def test_symmetry_permutations(fcc):
    lattice, positions, numbers, _ = fcc
    symprec = 1e-5
    dataset = get_symmetry_dataset((lattice, positions, numbers), symprec)
    rotations = dataset["rotations"]
    translations = dataset["translations"]
    permutations = get_symmetry_permutations(
        lattice, positions, numbers, rotations, translations, symprec
    )

    assert len(permutations) == len(rotations)
    for permutation in permutations:
        assert np.all(np.sort(permutation.permutation) == np.arange(len(positions)))
