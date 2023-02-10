from __future__ import annotations

import numpy as np
import pytest
from spglib import get_pointgroup, get_symmetry_from_database
from spgrep.pointgroup import pg_dataset
from spgrep.utils import is_integer_array

from spinspg.pointgroup import (
    POINT_GROUP_GENERATORS,
    POINT_GROUP_REPRESENTATIVES,
    SPIN_POINT_GROUP_TYPES,
    get_integer_point_group,
    get_pointgroup_representative,
    get_pointgroup_representative_from_symbol,
    traverse_spin_operations,
)


def test_point_group_representatives():
    assert set(POINT_GROUP_REPRESENTATIVES.keys()) == set(pg_dataset.keys())


def test_point_group_generators():
    identity = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for symbol, generators in POINT_GROUP_GENERATORS.items():
        group = get_pointgroup_representative_from_symbol(symbol)
        generated = traverse_spin_operations(
            [(identity, g) for g in [group[i] for i in generators]]
        )
        assert len(generated) == len(group)


def test_spin_point_group_table():
    founds = set()
    for symbol_R, datum_R in SPIN_POINT_GROUP_TYPES.items():
        R = get_pointgroup_representative_from_symbol(symbol_R)
        generator_R = [R[i] for i in POINT_GROUP_GENERATORS[symbol_R]]
        for symbol_r, datum_R_r in datum_R.items():  # type: ignore
            for symbol_B, datum_R_r_B in datum_R_r.items():
                B = get_pointgroup_representative_from_symbol(symbol_B)
                for number, mapping in datum_R_r_B:
                    print(f"#{number}: R={symbol_R}, r={symbol_r}, B={symbol_B}")
                    spg_generators = [(B[idx], generator_R[i]) for i, idx in enumerate(mapping)]
                    spg = traverse_spin_operations(spg_generators)

                    R2 = list({rot for _, rot in spg})
                    B2 = list({srot for srot, _ in spg})
                    r2 = [rot for srot, rot in spg if srot == ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
                    symbol_r2, _, _ = get_pointgroup(r2)

                    assert len(spg) == len(R)
                    assert set(R2) == set(R)
                    assert set(B2) == set(B)
                    assert symbol_r2 == symbol_r
                    assert spg not in founds
                    founds.add(spg)

    assert len(founds) == 598


@pytest.mark.parametrize(
    "symbol,idx",
    [
        ("mm2", 0),  # unique axis a -> unique axis c
        ("mm2", 1),  # unique axis b -> unique axis c
        ("mm2", 2),  # unique axis c -> unique axis c
        ("-42m", 0),  # -42m -> -42m
        ("-42m", 1),  # -4m2 -> -42m
        ("32", 0),  # 312 -> 312
        ("32", 1),  # 321 -> 312
        ("3m", 0),  # 3m1 -> 3m1
        ("3m", 1),  # 31m -> 3m1
        ("-3m", 0),  # -31m -> -31m
        ("-3m", 1),  # -3m1 -> -31m
        ("-6m2", 0),  # -6m2 -> -6m2
        ("-6m2", 1),  # -62m -> -6m2
    ],
)
def test_get_pointgroup_representative(symbol, idx):
    pg = np.array(pg_dataset[symbol][idx])[::-1]
    symbol_actual, P, mapping = get_pointgroup_representative(pg)

    assert symbol_actual == symbol

    idx_std = POINT_GROUP_REPRESENTATIVES[symbol]
    pg_std = np.array(pg_dataset[symbol][idx_std])
    assert len(pg_std) == len(pg)

    for i in range(len(pg)):
        ri_actual = np.linalg.inv(P) @ pg[mapping[i]] @ P
        assert is_integer_array(ri_actual)
        ri_actual = np.around(ri_actual).astype(int)
        assert np.allclose(ri_actual, pg_std[i])


def test_get_integer_point_group_3m():
    # 3m
    cart_rotations = np.array(
        [
            np.eye(3),
            [
                [np.cos(np.pi * 2 / 3), -np.sin(np.pi * 2 / 3), 0],
                [np.sin(np.pi * 2 / 3), np.cos(np.pi * 2 / 3), 0],
                [0, 0, 1],
            ],
            [
                [np.cos(np.pi * 4 / 3), -np.sin(np.pi * 4 / 3), 0],
                [np.sin(np.pi * 4 / 3), np.cos(np.pi * 4 / 3), 0],
                [0, 0, 1],
            ],
            np.diag([1, -1, 1]),
            [
                [np.cos(np.pi * 2 / 3), -np.sin(np.pi * 2 / 3), 0],
                [-np.sin(np.pi * 2 / 3), -np.cos(np.pi * 2 / 3), 0],
                [0, 0, 1],
            ],
            [
                [np.cos(np.pi * 4 / 3), -np.sin(np.pi * 4 / 3), 0],
                [-np.sin(np.pi * 4 / 3), -np.cos(np.pi * 4 / 3), 0],
                [0, 0, 1],
            ],
        ],
        dtype=np.float_,
    )
    P, rotations = get_integer_point_group(cart_rotations)

    for rot, cart_rot in zip(rotations, cart_rotations):
        rot2 = np.linalg.inv(P) @ cart_rot @ P
        assert is_integer_array(rot2)
        assert np.allclose(rot, rot2)


@pytest.mark.parametrize(
    "hall_number,pointgroup_international",
    [
        (1, "1"),
        (2, "-1"),
        (3, "2"),
        (18, "m"),
        (57, "2/m"),
        (108, "222"),
        (125, "mm2"),
        (227, "mmm"),
        (349, "4"),
        (355, "-4"),
        (357, "4/m"),
        (366, "422"),
        (376, "4mm"),
        (388, "-42m"),
        (400, "4/mmm"),
        (430, "3"),
        (435, "-3"),
        (438, "32"),
        (446, "3m"),
        (454, "-3m"),
        (462, "6"),
        (468, "-6"),
        (469, "6/m"),
        (471, "622"),
        (477, "6mm"),
        (481, "-6m2"),
        (485, "6/mmm"),
        (489, "23"),
        (494, "m-3"),
        (503, "432"),
        (511, "-43m"),
        (517, "m-3m"),
    ],
)
def test_get_integer_point_group(hall_number, pointgroup_international):
    symmetry = get_symmetry_from_database(hall_number)
    rotations_expect = symmetry["rotations"]

    # Generate random orthogonal matrix
    np.random.seed(0)
    P, _ = np.linalg.qr(np.random.rand(3, 3))

    # A hexagonal lattice
    P_hex = np.array(
        [
            [1, 0, 0],
            [-0.5, np.sqrt(3) / 2, 0],
            [0, 0, 1],
        ]
    ).T

    cart_rotations = []
    for rot in rotations_expect:
        if pointgroup_international in [
            "3",
            "-3",
            "32",
            "3m",
            "-3m",
            "6",
            "-6",
            "6/m",
            "622",
            "6mm",
            "-6m2",
            "6/mmm",
        ]:
            rot = P_hex @ rot @ np.linalg.inv(P_hex)

        cart_rot = P @ rot @ np.linalg.inv(P)
        assert np.allclose(cart_rot @ cart_rot.T, np.eye(3))
        cart_rotations.append(cart_rot)

    P_actual, rotations_actual = get_integer_point_group(cart_rotations)

    for rot, cart_rot in zip(rotations_actual, cart_rotations):
        rot2 = np.linalg.inv(P_actual) @ cart_rot @ P_actual
        assert is_integer_array(rot2)
        assert np.allclose(rot, rot2)
