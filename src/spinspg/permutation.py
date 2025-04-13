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
        n = len(self.permutation)
        assert len(rhs.permutation) == n
        mul = np.zeros(n, dtype=np.int_)
        for i in range(n):
            mul[i] = self(rhs(i))
        return Permutation(mul)


def get_symmetry_permutations(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    rotations: NDArrayInt,
    translations: NDArrayFloat,
    symprec: float,
) -> list[Permutation]:
    """Return permutations of sites from given symmetry operations."""
    num_sites = len(positions)

    permutations = []
    for rot, trans in zip(rotations, translations):
        new_positions = positions @ rot.T + trans[None, :]
        perm = [-1 for _ in range(num_sites)]
        found = [False for _ in range(num_sites)]
        for i in range(num_sites):
            for j in range(num_sites):
                if found[j] or (numbers[i] != numbers[j]):
                    continue
                if is_overlap_with_origin(lattice, new_positions[i] - positions[j], symprec):
                    perm[i] = j
                    found[j] = True
                    break

        if np.all(perm != -1):
            permutations.append(Permutation(perm))

    return permutations


def is_overlap_with_origin(lattice, frac_coords, symprec) -> bool:
    """Return true iff ``frac_coords`` is overlapped with the origin up to lattice translations."""
    diff = lattice.T @ (frac_coords - np.rint(frac_coords))
    return np.linalg.norm(diff) < symprec
