from __future__ import annotations

import numpy as np
import pytest
from spglib import get_pointgroup
from spgrep.pointgroup import pg_dataset
from spgrep.utils import is_integer_array

from spinspg.pointgroup import (
    POINT_GROUP_GENERATORS,
    POINT_GROUP_REPRESENTATIVES,
    SPIN_POINT_GROUP_TYPES,
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
