import numpy as np
import pytest

from spinspg.spin import SpinOnlyGroupType, get_spin_only_group, solve_procrustes


@pytest.fixture
def nonmagnetic():
    magmoms_nonmagnetic = np.array(
        [
            [0, 0, 0],
            [0, 0, 0],
        ],
        dtype=np.float_,
    )
    return magmoms_nonmagnetic


@pytest.fixture
def collinear():
    magmoms_collinear = np.array(
        [
            [0, 0, 1],
            [0, 0, -2],
        ],
        dtype=np.float_,
    )
    return magmoms_collinear


@pytest.fixture
def coplanar():
    magmoms_coplanar = np.array(
        [
            [0, 0, 1],
            [0, 1, 0],
        ],
        dtype=np.float_,
    )
    return magmoms_coplanar


@pytest.fixture
def noncoplanar():
    magmoms_noncoplanar = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        dtype=np.float_,
    )
    return magmoms_noncoplanar


@pytest.mark.parametrize(
    "magmoms_type,sog_type,axis",
    [
        ("nonmagnetic", SpinOnlyGroupType.NONMAGNETIC, None),
        ("collinear", SpinOnlyGroupType.COLLINEAR, [0, 0, 1]),
        ("coplanar", SpinOnlyGroupType.COPLANAR, [1, 0, 0]),
        ("noncoplanar", SpinOnlyGroupType.NONCOPLANAR, None),
    ],
)
def test_spin_only_group(request, magmoms_type, sog_type, axis):
    magmoms = request.getfixturevalue(magmoms_type)
    sog = get_spin_only_group(magmoms, mag_symprec=1e-5)
    assert sog.spin_only_group_type == sog_type
    if axis is not None:
        assert np.allclose(np.cross(sog.axis, axis), 0)

    # Test __str__
    assert str(sog)


FOURFOLD_Z = [
    [0, -1, 0],
    [1, 0, 0],
    [0, 0, 1],
]
MIRROR_X = [
    [-1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
]


@pytest.mark.parametrize(
    "magmoms_type,linear,expect",
    [
        ("nonmagnetic", FOURFOLD_Z, True),
        ("nonmagnetic", MIRROR_X, True),
        ("collinear", FOURFOLD_Z, True),
        ("collinear", MIRROR_X, True),
        ("coplanar", FOURFOLD_Z, False),
        ("coplanar", MIRROR_X, True),
        ("noncoplanar", FOURFOLD_Z, False),
        ("noncoplanar", MIRROR_X, False),
    ],
)
def test_contain(request, magmoms_type, linear, expect):
    magmoms = request.getfixturevalue(magmoms_type)
    sog = get_spin_only_group(magmoms, mag_symprec=1e-5)
    assert sog.contain(linear) == expect


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
