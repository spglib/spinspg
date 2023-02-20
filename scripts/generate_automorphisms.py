import json
from pathlib import Path

import numpy as np
from spgrep.utils import is_integer_array, ndarray2d_to_integer_tuple

from spinspg.pointgroup import (
    POINT_GROUP_REPRESENTATIVES,
    get_pointgroup_representative_from_symbol,
)


def get_automorphism_permutation(group, g):
    conjugated_group = []
    for h in group:
        conj = np.linalg.inv(g) @ h @ g
        assert is_integer_array(conj)
        conjugated_group.append(ndarray2d_to_integer_tuple(conj))

    assert set(conjugated_group) == set(group)
    mapping = [-1 for _ in range(len(group))]
    for i, hi in enumerate(group):
        j = conjugated_group.index(hi)  # type: ignore
        mapping[j] = i

    return mapping


AUT_D4h = np.array(
    [
        # x, y, z
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ],
        # -y, x, z
        [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
        # -x, y, -z
        [
            [-1, 0, 0],
            [0, 1, 0],
            [0, 0, -1],
        ],
        # y, x, -z
        [
            [0, 1, 0],
            [1, 0, 0],
            [0, 0, -1],
        ],
    ]
)

AUT_D6h = np.array(
    [
        # x, y, z
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        # -y, x-y, z
        [
            [0, -1, 0],
            [1, -1, 0],
            [0, 0, 1],
        ],
        # -x+y, -x, z
        [
            [-1, 1, 0],
            [-1, 0, 0],
            [0, 0, 1],
        ],
        # y, x, -z
        [
            [0, 1, 0],
            [1, 0, 0],
            [0, 0, -1],
        ],
        # -x+y, y, -z
        [
            [-1, 1, 0],
            [0, 1, 0],
            [0, 0, -1],
        ],
        # x, x-y, -z
        [
            [1, 0, 0],
            [1, -1, 0],
            [0, 0, -1],
        ],
    ]
)


def main():
    symbols = POINT_GROUP_REPRESENTATIVES.keys()
    all_datum = {}
    for symbol in symbols:
        print(symbol)
        if symbol in ["222", "mmm"]:
            # Aut:
            # Aut: | m-3m / mmm | = 6
            aut = np.array(
                [
                    # x, y, z
                    [
                        [1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1],
                    ],
                    # z, x, y
                    [
                        [0, 0, 1],
                        [1, 0, 0],
                        [0, 1, 0],
                    ],
                    # y, z, x
                    [
                        [0, 1, 0],
                        [0, 0, 1],
                        [1, 0, 0],
                    ],
                    # y, x, -z
                    [
                        [0, 1, 0],
                        [1, 0, 0],
                        [0, 0, -1],
                    ],
                    # x, z, -y
                    [
                        [1, 0, 0],
                        [0, 0, 1],
                        [0, -1, 0],
                    ],
                    # z, y, -x
                    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
                ]
            )
        elif symbol in ["23", "m-3", "432", "-43m", "m-3m"]:
            # Aut: 432
            aut = np.array(get_pointgroup_representative_from_symbol("432"))
        elif symbol == "mm2":
            # Aut: | 4/mmm / mmm | = 2
            aut = np.array(
                [
                    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],  # 1
                    [
                        [0, -1, 0],
                        [1, 0, 0],
                        [0, 0, 1],
                    ],  # 4^+
                ]
            )
        elif symbol == "-42m":
            # Aut: | 4/mmm / 2/m | = 4
            aut = AUT_D4h
        elif symbol in ["422", "4mm", "4/mmm"]:
            # Aut: | 8/mmm / 2/m | = 8
            C8 = np.array(
                [
                    [np.cos(np.pi / 4), -np.sin(np.pi / 4), 0],
                    [np.sin(np.pi / 4), np.cos(np.pi / 4), 0],
                    [0, 0, 1],
                ]
            )
            aut = np.concatenate([AUT_D4h, [C8 @ g for g in AUT_D4h]])
        elif symbol in ["32", "3m", "-3m", "-6m2"]:
            # Aut: | 6/mmm / 2/m | = 6
            aut = AUT_D6h
        elif symbol in ["622", "6mm", "6/mmm"]:
            # Aut: 12 22
            C12_hex = np.array(
                [
                    [2 / np.sqrt(3), -1 / np.sqrt(3), 0],
                    [1 / np.sqrt(3), 1 / np.sqrt(3), 0],
                    [0, 0, 1],
                ]
            )
            aut = np.concatenate([AUT_D6h, [C12_hex @ g for g in AUT_D6h]])
        elif symbol in ["1", "-1", "2", "m", "2/m", "4", "-4", "4/m", "3", "-3", "6", "-6", "6/m"]:
            aut = np.array(get_pointgroup_representative_from_symbol("1"))
        else:
            raise ValueError("Unreachable!")

        group = get_pointgroup_representative_from_symbol(symbol)
        permutations = [get_automorphism_permutation(group, g) for g in aut]
        unique_perms = {tuple(p) for p in permutations}
        assert len(permutations) == len(unique_perms)

        all_datum[symbol] = [(g.tolist(), perm) for g, perm in zip(aut, permutations)]

    path = Path(__file__).parent.parent / "src" / "spinspg" / "automorphisms.json"
    with open(path, "w") as f:
        json.dump(all_datum, f)


if __name__ == "__main__":
    main()
