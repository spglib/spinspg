from __future__ import annotations

import numpy as np

from spinspg.group import get_primitive_spin_symmetry, get_symmetry_with_cell
from spinspg.spin import SpinOnlyGroup
from spinspg.utils import NDArrayFloat, NDArrayInt


def get_spin_symmetry(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    magmoms: NDArrayFloat,
    symprec: float = 1e-5,
    angle_tolerance: float = -1.0,
) -> tuple[SpinOnlyGroup, NDArrayInt, NDArrayFloat, NDArrayFloat]:
    """Return spin symmetry operations of given spin arrangement.

    Parameters
    ----------
    lattice: array, (3, 3)
    positions: array, (num_sites, 3)
    numbers: array[int], (num_sites, )
    magmoms: array, (num_sites, 3)
        in Cartesian coordinates
    symprec: float
    angle_tolerance: float

    Returns
    -------
    spin_only_group: SpinOnlyGroup
    rotations: array[int], (num_sym, 3, 3)
        Rotation parts of spin symmetry operations
    translations: array, (num_sym, 3)
        Translation parts of spin symmetry operations
    spin_rotations: array, (num_sym, 3, 3)
        Spin rotation parts of spin symmetry operations
    """
    ns = get_symmetry_with_cell(lattice, positions, numbers, symprec, angle_tolerance)
    ssg = get_primitive_spin_symmetry(ns, magmoms, symprec)

    spin_only_group = ssg.spin_only_group
    tmat = ssg.transformation
    invtmat = np.linalg.inv(tmat)
    rotations = []
    translations = []
    spin_rotations = []
    for ops in ssg.nontrivial_coset:
        # Transform to primitive to input cell
        new_rotation = np.around(invtmat @ ops.rotation @ tmat).astype(np.int_)
        for centering in ssg.prim_centerings:
            new_translation = np.remainder(invtmat @ (ops.translation + centering), 1)
            rotations.append(new_rotation)
            translations.append(new_translation)
            spin_rotations.append(ops.spin_rotation)

    return spin_only_group, np.array(rotations), np.array(translations), np.array(spin_rotations)
