"""Core APIs."""
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
    """Return spin symmetry operations of a given spin arrangement.

    See :ref:`Spglib's document <spglib:py_variables_crystal_structure>` for how to specify the spin arrangement by ``lattice``, ``positions``, ``numbers``, and ``magmoms`` in details.

    Parameters
    ----------
    lattice: array, (3, 3)
        ``lattice[i, :]`` is the ``i``-th basis vector of a lattice
    positions: array, (num_sites, 3)
        ``positions[i, :]`` is a fractional coordinates of the ``i``-th site w.r.t. ``lattice``.
    numbers: array[int], (num_sites, )
        ``numbers[i]`` specifies a specie at the ``i``-th site.
    magmoms: array, (num_sites, 3)
        ``magmoms[i, :]`` is a magnetic moments at the ``i``-th site in Cartesian coordinates.
    symprec: float, default=1e-5
        See :ref:`spglib:variables_symprec`.
    angle_tolerance: float, default=-1
        See :ref:`spglib:variables_angle_tolerance`.

    Returns
    -------
    spin_only_group: :class:`spin.SpinOnlyGroup`
    rotations: array[int], (num_sym, 3, 3)
        Rotation parts of spin symmetry operations w.r.t. ``lattice``.
    translations: array, (num_sym, 3)
        Translation parts of spin symmetry operations w.r.t. ``lattice``.
    spin_rotations: array, (num_sym, 3, 3)
        Spin rotation parts of spin symmetry operations in Cartesian coordinates.
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
