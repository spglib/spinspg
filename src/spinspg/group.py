"""Group for spin symmetry operations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from hsnf import column_style_hermite_normal_form
from spglib import get_symmetry_dataset

from spinspg.permutation import Permutation, get_symmetry_permutations
from spinspg.spin import SpinOnlyGroup, SpinOnlyGroupType, get_spin_only_group, solve_procrustes
from spinspg.utils import (
    NDArrayFloat,
    NDArrayInt,
    is_integer_array,
    ndarray2d_to_integer_tuple,
)


@dataclass
class NonmagneticSymmetry:
    """Crystal structure with symmetry operations.

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
    prim_centerings: array[int], (nc, 3)
        Centering translations w.r.t. ``prim_lattice``
    prim_centering_permutations: (nc, N)
    transformation: array[int], (3, 3)
        Transformation matrix from primitive to given cell
    """

    prim_lattice: NDArrayFloat
    prim_rotations: NDArrayInt
    prim_translations: NDArrayFloat
    prim_permutations: list[Permutation]
    prim_centerings: NDArrayInt
    prim_centering_permutations: list[Permutation]
    transformation: NDArrayInt


def get_symmetry_with_cell(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    symprec: float,
    angle_tolerance: float,
) -> NonmagneticSymmetry:
    """Find spatial symmetry operations from nonmagnetic crystal structure."""
    dataset = get_symmetry_dataset((lattice, positions, numbers), symprec, angle_tolerance)
    rotations = dataset.rotations
    translations = dataset.translations
    prim_lattice = dataset.primitive_lattice

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
    assert np.isclose(np.abs(np.linalg.det(tmat)), len(centerings))

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
        prim_centering = tmat @ centering
        assert is_integer_array(prim_centering)
        prim_centerings.append(np.around(prim_centering).astype(int))

    return NonmagneticSymmetry(
        prim_lattice=prim_lattice,
        prim_rotations=np.array(prim_rotations),
        prim_translations=np.array(prim_translations),
        prim_permutations=prim_permutations,
        prim_centerings=np.array(prim_centerings),
        prim_centering_permutations=prim_centering_permutations,
        transformation=np.around(tmat).astype(np.int_),
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
    centerings: array, (?, 3)
    nontrivial_coset: list[SpinSymmetryOperation]
        N.B. rotation parts are distinct
    transformation: array[int], (3, 3)
        Transformation matrix from primitive to given cell
    """

    prim_lattice: NDArrayFloat
    spin_only_group: SpinOnlyGroup
    spin_translation_coset: list[SpinSymmetryOperation]
    prim_centerings: NDArrayFloat
    nontrivial_coset: list[SpinSymmetryOperation]
    transformation: NDArrayInt


def get_primitive_spin_symmetry(
    nonmagnetic_symmetry: NonmagneticSymmetry, magmoms: NDArrayFloat, mag_symprec: float
) -> SpinSpaceGroup:
    """Return spin space group symmetry.

    Parameters
    ----------
    nonmagnetic_symmetry : NonmagneticSymmetry
    magmoms : array, (num_sites, 3)
    mag_symprec : float

    Returns
    -------
    SpinSpaceGroup
    """
    # Centerings for maximal space subgroup of spin space group
    stg_centerings = []
    stg_centering_permutations = []
    for centering, perm in zip(
        nonmagnetic_symmetry.prim_centerings, nonmagnetic_symmetry.prim_centering_permutations
    ):
        if np.max(np.linalg.norm(magmoms[perm.permutation] - magmoms, axis=1)) < mag_symprec:
            stg_centerings.append(centering)
            stg_centering_permutations.append(perm)
    assert len(nonmagnetic_symmetry.prim_centerings) % len(stg_centerings) == 0

    # Transformation matrix to primitive cell of maximal space subgroup
    stg_vectors = np.concatenate(
        [
            nonmagnetic_symmetry.transformation,  # (3, 3)
            np.array(stg_centerings).T,  # (3, ?)
        ],
        axis=1,
    )
    tmat_stg, _ = column_style_hermite_normal_form(stg_vectors)
    tmat_stg = tmat_stg[:, :3]  # (3, 3)
    invtmat_stg = np.linalg.inv(tmat_stg)

    # Spin only group
    spin_only_group = get_spin_only_group(magmoms, mag_symprec)

    prim_centerings = [invtmat_stg @ centering for centering in stg_centerings]
    # prim_spin_lattice.T @ transformation == nonmagnetic_symmetry.prim_lattice.T @ nonmagnetic_symmetry.transformation
    transformation = invtmat_stg @ nonmagnetic_symmetry.transformation
    assert is_integer_array(transformation)
    transformation = np.around(transformation).astype(np.int_)

    # Spin translation group search
    spin_translation_coset = []
    found_stg_centerings = []  # type: ignore
    for centering, perm in zip(
        nonmagnetic_symmetry.prim_centerings, nonmagnetic_symmetry.prim_centering_permutations
    ):
        # Two centerings are equivalent in primitive cell of spin translation group if they
        # are translated to each other by lattice translations in `transformation`.
        is_new_centering = True
        for other in found_stg_centerings:
            residual = invtmat_stg @ (centering - other)
            residual -= np.rint(residual)
            if np.allclose(residual, 0):
                is_new_centering = False
                break
        if not is_new_centering:
            continue
        found_stg_centerings.append(centering)

        # Search W in O(3) s.t. new_magmoms @ W.T = magmoms[perm.permutation]
        new_magmoms = magmoms.copy()
        perm_magmoms = magmoms[perm.permutation]
        W = solve_procrustes(new_magmoms, perm_magmoms)
        W = purify_spin_rotation(W, spin_only_group)

        new_magmoms = new_magmoms @ W.T
        if np.max(np.linalg.norm(new_magmoms - perm_magmoms, axis=1)) < mag_symprec:
            # w.r.t. primitive cell of spin space group
            reduced_centering = invtmat_stg @ centering
            spin_translation_coset.append(
                SpinSymmetryOperation(
                    rotation=np.eye(3, dtype=np.int_),
                    translation=reduced_centering,
                    spin_rotation=W,
                )
            )

    assert len(nonmagnetic_symmetry.prim_centerings) % len(found_stg_centerings) == 0

    # Spin space group search
    nontrivial_coset = []
    for rot, trans, perm in zip(
        nonmagnetic_symmetry.prim_rotations,
        nonmagnetic_symmetry.prim_translations,
        nonmagnetic_symmetry.prim_permutations,
    ):
        # Point group symmetry compatible with the primitive cell
        rot_prim = invtmat_stg @ rot @ tmat_stg
        if not is_integer_array(rot_prim):
            continue

        # Need to consider centerings for subgroup
        for centering, centering_perm in zip(
            nonmagnetic_symmetry.prim_centerings, nonmagnetic_symmetry.prim_centering_permutations
        ):
            new_perm: Permutation = centering_perm * perm
            new_magmoms = magmoms.copy()
            perm_magmoms = magmoms[new_perm.permutation]
            W = solve_procrustes(new_magmoms, perm_magmoms)
            W = purify_spin_rotation(W, spin_only_group)

            new_magmoms = new_magmoms @ W.T
            if np.max(np.linalg.norm(new_magmoms - perm_magmoms, axis=1)) < mag_symprec:
                # w.r.t. primitive cell of spin space group
                new_trans = centering + trans
                nontrivial_coset.append(
                    SpinSymmetryOperation(
                        rotation=rot_prim,
                        translation=invtmat_stg @ new_trans,
                        spin_rotation=W,
                    )
                )
                break

    # Transform centerings to primitive cell of spin space group
    prim_spin_lattice = tmat_stg.T @ nonmagnetic_symmetry.prim_lattice

    return SpinSpaceGroup(
        prim_lattice=prim_spin_lattice,
        spin_only_group=spin_only_group,
        spin_translation_coset=spin_translation_coset,
        prim_centerings=prim_centerings,
        nontrivial_coset=nontrivial_coset,
        transformation=transformation,
    )


def purify_spin_rotation(
    W: NDArrayFloat, spin_only_group: SpinOnlyGroup, atol: float = 1e-5
) -> NDArrayFloat:
    """Purify spin rotation matrix `W` by eliminating the spin-only group part."""
    identity = np.eye(3, dtype=np.float64)
    if spin_only_group.spin_only_group_type == SpinOnlyGroupType.NONMAGNETIC:
        # spin_only_group = O(3): All spin rotations are equivalent to identity
        return identity
    elif spin_only_group.spin_only_group_type == SpinOnlyGroupType.COLLINEAR:
        # normalizer of spin_only_group: infty/mm = (infty m) x m (ITA Table 3.5.4.2)
        axis = spin_only_group.axis
        assert isinstance(axis, np.ndarray)
        if np.allclose(W @ axis, axis, atol=atol):
            return identity
        elif np.allclose(W @ axis, -axis, atol=atol):
            # mirror along axis
            return _get_mirror_along_axis(axis)
        else:
            raise ValueError("Given spin rotation does not normalize the spin-only group.")
    elif spin_only_group.spin_only_group_type == SpinOnlyGroupType.COPLANAR:
        # normalizer of spin_only_group: infty/mm
        mirror = _get_mirror_along_axis(spin_only_group.axis)
        if np.linalg.det(W) > 0:
            return W
        else:
            return W @ mirror
    elif spin_only_group.spin_only_group_type == SpinOnlyGroupType.NONCOPLANAR:
        return W


def _get_mirror_along_axis(axis: NDArrayFloat) -> NDArrayFloat:
    """Return mirror matrix along the given axis."""
    return np.eye(3, dtype=np.float64) - 2 * axis[:, None] @ axis[None, :]
