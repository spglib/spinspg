from spinspg.group import get_symmetry_with_cell


def test_get_symmetry_with_cell(fcc):
    lattice, positions, numbers, _ = fcc
    symmetry = get_symmetry_with_cell(lattice, positions, numbers, 1e-5, -1)
    assert symmetry.prim_rotations.shape == (48, 3, 3)
    assert symmetry.prim_centerings.shape == (4, 3)
