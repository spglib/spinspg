"""Permutations from action of symmetry operation on sites."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from typing_extensions import Self

if TYPE_CHECKING:
    from spinspg.utils import NDArrayFloat, NDArrayInt


@dataclass
class Permutation:
    """Permutation of list."""

    permutation: NDArrayInt

    def __call__(self, idx: int) -> int:
        """Return a permuted index of ``idx``."""
        return self.permutation[idx]

    def __mul__(self, rhs: Self) -> Permutation:
        """Return product with a given permutation ``rhs``.

        (self * rhs)(i) = self(rhs(i))
        """
        assert len(rhs.permutation) == len(self.permutation)
        return Permutation(np.asarray(self.permutation)[rhs.permutation])


def get_symmetry_permutations(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    symprec: float,
) -> list[Permutation]:
    """Return site permutations in the same order as the supplied operations.

    Raise ``ValueError`` if an operation cannot match every site. Matching is
    greedy in input-site order with a strict Cartesian distance threshold.
    """
    num_sites = len(positions)
    sites_by_species = {
        number: np.flatnonzero(numbers == number).tolist() for number in np.unique(numbers)
    }

    permutations = []
    for rot, trans in zip(rotations, translations):
        new_positions = positions @ rot.T + trans[None, :]
        perm = np.empty(num_sites, dtype=np.int_)
        found = [False for _ in range(num_sites)]
        for i in range(num_sites):
            for j in sites_by_species[numbers[i]]:
                if found[j]:
                    continue
                if is_overlap_with_origin(lattice, new_positions[i] - positions[j], symprec):
                    # Keep the first available matching site, not the nearest site.
                    perm[i] = j
                    found[j] = True
                    break
            else:
                raise ValueError("Symmetry operation does not map all sites within symprec")
        permutations.append(Permutation(perm))

    return permutations


def is_overlap_with_origin(lattice, frac_coords, symprec) -> bool:
    """Return true iff ``frac_coords`` is overlapped with the origin up to lattice translations."""
    diff = lattice.T @ (frac_coords - np.rint(frac_coords))
    return np.linalg.norm(diff) < symprec
