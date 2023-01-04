import numpy as np

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
