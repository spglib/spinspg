from spinspg.group import get_spin_space_group, get_symmetry_with_cell
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
    ssg = get_spin_space_group(ns, magmoms, symprec)

    assert ssg.spin_only_group.spin_only_group_type == SpinOnlyGroupType.COPLANAR
    assert len(ssg.spin_translation_coset) == 2
