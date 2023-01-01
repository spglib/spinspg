import numpy as np

from spinspg.spin import SpinOnlyGroupType, get_spin_only_group, solve_procrustes


def test_spin_only_group():
    mag_symprec = 1e-5

    # Nonmagnetic
    magmoms_nonmagnetic = np.array(
        [
            [0, 0, 0],
            [0, 0, 0],
        ],
        dtype=np.float_,
    )
    so_nonmagnetic = get_spin_only_group(magmoms_nonmagnetic, mag_symprec)
    assert so_nonmagnetic.spin_only_group_type == SpinOnlyGroupType.NONMAGNETIC

    # Collinear
    magmoms_collinear = np.array(
        [
            [0, 0, 1],
            [0, 0, -2],
        ],
        dtype=np.float_,
    )
    so_collinear = get_spin_only_group(magmoms_collinear, mag_symprec)
    assert so_collinear.spin_only_group_type == SpinOnlyGroupType.COLLINEAR
    assert np.allclose(np.cross(so_collinear.axis, np.array([0, 0, 1])), 0)

    # Coplanar
    magmoms_collinear = np.array(
        [
            [0, 0, 1],
            [0, 1, 0],
        ],
        dtype=np.float_,
    )
    so_coplanar = get_spin_only_group(magmoms_collinear, mag_symprec)
    assert so_coplanar.spin_only_group_type == SpinOnlyGroupType.COPLANAR
    assert np.allclose(np.cross(so_coplanar.axis, np.array([1, 0, 0])), 0)

    # Non-coplanar
    magmoms_noncoplanar = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        dtype=np.float_,
    )
    so_noncoplanar = get_spin_only_group(magmoms_noncoplanar, mag_symprec)
    assert so_noncoplanar.spin_only_group_type == SpinOnlyGroupType.NONCOPLANAR


def test_solve_procrustes():
    A = np.array(
        [
            [1, 0, 0],
            [1, 1, 0],
            [1, 1, 1],
            [0, 0, 1],
        ]
    )
    # Random orthogonal matrix:
    #   a = np.random.random((3, 3))
    #   R = sp.linalg.polar(a)[0]
    R = np.array(
        [
            [-0.03831202, 0.99710779, 0.06563723],
            [0.28470653, -0.05207084, 0.95719947],
            [0.95784883, 0.05535959, -0.28188816],
        ]
    )
    B = A @ R.T

    R_actual = solve_procrustes(A, B)
    assert np.allclose(R_actual, R)
