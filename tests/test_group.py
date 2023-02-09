import numpy as np
import pytest
from spglib import get_magnetic_symmetry

from spinspg.core import get_spin_symmetry
from spinspg.group import get_primitive_spin_symmetry, get_symmetry_with_cell
from spinspg.spin import SpinOnlyGroupType


def test_get_symmetry_with_cell(fcc):
    lattice, positions, numbers, _ = fcc
    symmetry = get_symmetry_with_cell(lattice, positions, numbers, 1e-5, -1)
    assert symmetry.prim_rotations.shape == (48, 3, 3)
    assert symmetry.prim_centerings.shape == (4, 3)


def test_spin_space_group_fcc(fcc):
    lattice, positions, numbers, magmoms = fcc
    symprec = 1e-5
    ns = get_symmetry_with_cell(lattice, positions, numbers, symprec, -1)
    ssg = get_primitive_spin_symmetry(ns, magmoms, symprec)

    assert ssg.spin_only_group.spin_only_group_type == SpinOnlyGroupType.COPLANAR
    assert len(ssg.spin_translation_coset) == 2


def test_spin_space_group_kagome(layer_triangular_kagome):
    lattice, positions, numbers, magmoms = layer_triangular_kagome
    symprec = 1e-5
    ns = get_symmetry_with_cell(lattice, positions, numbers, symprec, -1)
    ssg = get_primitive_spin_symmetry(ns, magmoms, symprec)

    assert ssg.spin_only_group.spin_only_group_type == SpinOnlyGroupType.COPLANAR
    assert len(ssg.spin_translation_coset) == 1

    kernel_pointgroup = []
    for ops in ssg.nontrivial_coset:
        if np.allclose(ops.spin_rotation, np.eye(3)):
            kernel_pointgroup.append(ops.rotation)
    assert len(kernel_pointgroup) == 4  # 2/m
    assert len(ssg.nontrivial_coset) == 24  # 6/mmm


def test_spin_space_group_rutile(rutile):
    lattice, positions, numbers, magmoms = rutile
    symprec = 1e-5
    ns = get_symmetry_with_cell(lattice, positions, numbers, symprec, -1)
    ssg = get_primitive_spin_symmetry(ns, magmoms, symprec)

    assert ssg.spin_only_group.spin_only_group_type == SpinOnlyGroupType.COLLINEAR
    assert np.allclose(np.cross(ssg.spin_only_group.axis, [0, 0, 1]), 0)
    assert len(ssg.spin_translation_coset) == 1
    assert len(ssg.prim_centerings) == 1
    assert np.isclose(np.linalg.det(ssg.transformation), 1)

    kernel_pointgroup = []
    for ops in ssg.nontrivial_coset:
        if np.allclose(ops.spin_rotation, np.eye(3)):
            kernel_pointgroup.append(ops.rotation)

    assert len(ssg.nontrivial_coset) == 16  # 4/mmm
    assert len(kernel_pointgroup) == 8  # 4/m


@pytest.mark.parametrize(
    "testcase,spin_only_group_type,axis",
    [
        ("Cr_in_Cr2O3", SpinOnlyGroupType.COLLINEAR, [1, 0, 0]),
        ("Pr_in_PrScSb", SpinOnlyGroupType.COLLINEAR, [0, 0, 1]),
        ("Mn_in_Mn3ReO6", SpinOnlyGroupType.COPLANAR, [0, 0, 1]),
        ("Ni_in_NiTa2O6", SpinOnlyGroupType.COLLINEAR, [1, -1, 0]),
    ],
)
def test_spin_space_groups(request, testcase, spin_only_group_type, axis):
    lattice, positions, numbers, magmoms = request.getfixturevalue(testcase)
    symprec = 1e-5
    ns = get_symmetry_with_cell(lattice, positions, numbers, symprec, -1)
    ssg = get_primitive_spin_symmetry(ns, magmoms, symprec)

    assert ssg.spin_only_group.spin_only_group_type == spin_only_group_type
    assert np.allclose(np.cross(ssg.spin_only_group.axis, axis), 0)

    num_sym = (
        len(ssg.nontrivial_coset)
        * len(ssg.spin_translation_coset)
        * np.around(np.abs(np.linalg.det(ssg.transformation))).astype(int)
    )
    mag_symmetry = get_magnetic_symmetry((lattice, positions, numbers, magmoms), symprec=symprec)
    assert num_sym >= len(mag_symmetry["rotations"])


def test_get_spin_symmetry(rutile):
    lattice, positions, numbers, magmoms = rutile
    sog, rotations, translations, spin_rotations = get_spin_symmetry(
        lattice, positions, numbers, magmoms
    )

    assert sog.spin_only_group_type == SpinOnlyGroupType.COLLINEAR
    assert np.allclose(np.cross(sog.axis, [0, 0, 1]), 0)  # parallel to [0, 0, 1]

    identity = np.eye(3)
    mz = np.diag([1, 1, -1])
    expects = [
        (
            # (1) x, y, z; 1
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [0, 0, 0],
            identity,
        ),
        (
            # (2) -x, -y, z; 1
            [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
            [0, 0, 0],
            identity,
        ),
        (
            # (3) -y + 1/2, x + 1/2, z + 1/2; mz
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
            [0.5, 0.5, 0.5],
            mz,
        ),
        (
            # (4) y + 1/2, -x + 1/2, z + 1/2; mz
            [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
            [0.5, 0.5, 0.5],
            mz,
        ),
        (
            # (5) -x + 1/2, y + 1/2, -z + 1/2; mz
            [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
            [0.5, 0.5, 0.5],
            mz,
        ),
        (
            # (6) x + 1/2, -y + 1/2, -z + 1/2; mz
            [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
            [0.5, 0.5, 0.5],
            mz,
        ),
        (
            # (7) y, x, -z; 1
            [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
            [0, 0, 0],
            identity,
        ),
        (
            # (8) -y, -x, -z; 1
            [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],
            [0, 0, 0],
            identity,
        ),
        (
            # (9) -x, -y, -z; 1
            [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
            [0, 0, 0],
            identity,
        ),
        (
            # (10) x, y, -z; 1
            [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
            [0, 0, 0],
            identity,
        ),
        (
            # (11) y + 1/2, -x + 1/2, -z + 1/2; mz
            [[0, 1, 0], [-1, 0, 0], [0, 0, -1]],
            [0.5, 0.5, 0.5],
            mz,
        ),
        (
            # (12) -y + 1/2, x + 1/2, -z + 1/2; mz
            [[0, -1, 0], [1, 0, 0], [0, 0, -1]],
            [0.5, 0.5, 0.5],
            mz,
        ),
        (
            # (13) x + 1/2, -y + 1/2, z + 1/2; mz
            [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
            [0.5, 0.5, 0.5],
            mz,
        ),
        (
            # (14) -x + 1/2, y + 1/2, z + 1/2; mz
            [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [0.5, 0.5, 0.5],
            mz,
        ),
        (
            # (15) -y, -x, z; 1
            [[0, -1, 0], [-1, 0, 0], [0, 0, 1]],
            [0, 0, 0],
            identity,
        ),
        (
            # (16) y, x, z; 1
            [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
            [0, 0, 0],
            identity,
        ),
    ]

    found = [False for _ in expects]
    assert len(rotations) == len(expects)
    for rot, trans, srot in zip(rotations, translations, spin_rotations):
        for i in range(len(expects)):
            if (not found[i]) and np.allclose(rot, expects[i][0]):
                assert np.allclose(trans, expects[i][1])
                assert np.allclose(srot, expects[i][2])
                found[i] = True

    assert all(found)
