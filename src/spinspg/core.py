"""Core APIs."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from spinspg.group import get_primitive_spin_symmetry, get_symmetry_with_cell

if TYPE_CHECKING:
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
    With the returned spin symmetry operations, the ``i``-th spin symmetry operation maps point coordinates ``x`` to ``rotations[i] @ x + translations[i]`` and magnetic moment ``m`` to ``spin_rotations[i] @ m``.

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
        Spin-only group, which preserve the given spin arrangement with identity spatial operation.
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

    # Products of "translations in cell", "nontrivial spin translation group's coset", and "nontrivial spin space group's coset"
    for ops in ssg.nontrivial_coset:
        # Transform to primitive to input cell
        new_rotation = np.around(invtmat @ ops.rotation @ tmat).astype(np.int_)
        for ops_st in ssg.spin_translation_coset:
            for centering in ssg.prim_centerings:
                new_translation = np.remainder(
                    invtmat @ (ops.translation + ops_st.translation + centering), 1
                )
                rotations.append(new_rotation)
                translations.append(new_translation)
                spin_rotations.append(ops_st.spin_rotation @ ops.spin_rotation)

    return spin_only_group, np.array(rotations), np.array(translations), np.array(spin_rotations)
