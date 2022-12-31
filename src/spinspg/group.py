from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from spglib import get_symmetry_dataset

from spinspg.permutation import Permutation, get_symmetry_permutations
from spinspg.spin_only import SpinOnlyGroup
from spinspg.utils import NDArrayFloat, NDArrayInt, ndarray2d_to_integer_tuple


@dataclass
class NonmagneticSymmetry:
    """Crystal structure with symmetry operations

    Attributes
    ----------
    prim_lattice: array, (3, 3)
        Primitive basis vectors w.r.t. nonmagnetic symmetry
    prim_rotations: array[int], (order, 3, 3)
        w.r.t. ``prim_lattice``
    prim_translations: array, (order, 3)
        w.r.t. ``prim_lattice``
    prim_permutations: array[int], (order, num_sites)
        ``num_sites`` is a number of sites in input cell.
        ``(prim_rotations[p], prim_translations[p])`` moves the ``i``-th site to the ``prim_permutations[p][i]``
    centerings: array, (nc, 3)
        Centering translations w.r.t. ``prim_lattice``
    centering_permutations: (nc, N)
    """

    prim_lattice: NDArrayFloat
    prim_rotations: NDArrayInt
    prim_translations: NDArrayFloat
    prim_permutations: list[Permutation]
    prim_centerings: NDArrayFloat
    prim_centering_permutations: list[Permutation]


def get_symmetry_with_cell(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    symprec: float,
    angle_tolerance: float,
) -> NonmagneticSymmetry:
    dataset = get_symmetry_dataset((lattice, positions, numbers), symprec, angle_tolerance)
    rotations = dataset["rotations"]
    translations = dataset["translations"]
    prim_lattice = dataset["primitive_lattice"]

    # Unique by rotation parts
    uniq_rotations = []
    uniq_translations = []
    centerings = []
    found_rotations = set()
    for rot, trans in zip(rotations, translations):
        if np.allclose(rot, np.eye(3)):
            centerings.append(trans)

        rot_int = ndarray2d_to_integer_tuple(rot)
        if rot_int in found_rotations:
            continue
        uniq_rotations.append(rot)
        uniq_translations.append(trans)
        found_rotations.add(rot_int)

    # Primitive transformation
    tmat = np.linalg.inv(prim_lattice.T) @ lattice.T
    assert np.allclose(np.abs(np.linalg.det(tmat)), len(centerings))

    # Permutations of sites
    prim_permutations = get_symmetry_permutations(
        lattice,
        positions,
        numbers,
        rotations=uniq_rotations,
        translations=uniq_translations,
        symprec=symprec,
    )
    prim_centering_permutations = get_symmetry_permutations(
        lattice,
        positions,
        numbers,
        rotations=[np.eye(3) for _ in range(len(centerings))],
        translations=centerings,
        symprec=symprec,
    )

    # To primitive basis (never take modulus!)
    prim_rotations = []
    prim_translations = []
    invtmat = np.linalg.inv(tmat)
    for rot, trans in zip(uniq_rotations, uniq_translations):
        prim_rotations.append(np.around(tmat @ rot @ invtmat).astype(np.int_))
        prim_translations.append(tmat @ trans)

    prim_centerings = []
    for centering in centerings:
        prim_centerings.append(tmat @ centering)

    return NonmagneticSymmetry(
        prim_lattice=prim_lattice,
        prim_rotations=np.array(prim_rotations),
        prim_translations=np.array(prim_translations),
        prim_permutations=prim_permutations,
        prim_centerings=np.array(prim_centerings),
        prim_centering_permutations=prim_centering_permutations,
    )


@dataclass
class SpinSymmetryOperation:
    """Spin symmetry operation.

    Attributes
    ----------
    rotations: array[int], (3, 3)
    translation: array, (3, )
    spin_rotation: array, (3, 3)
    """

    rotation: NDArrayInt
    translation: NDArrayFloat
    spin_rotation: NDArrayFloat


@dataclass
class SpinSpaceGroup:
    """Spin space group.

    Attributes
    ----------
    primitive_lattice: array, (3, 3)
        primitive_lattice[i] is the i-th primitive basis vector for maximal space subgroup
    spin_only_group: SpinOnlyGroup
    spin_translation_coset: list[SpinSymmetryOperation]
        N.B. translation parts are distinct
    nontrivial_coset: list[SpinSymmetryOperation]
        N.B. rotation parts are distinct
    """

    prim_lattice: NDArrayFloat
    spin_only_group: SpinOnlyGroup
    spin_translation_coset: list[SpinSymmetryOperation]
    nontrivial_coset: list[SpinSymmetryOperation]
