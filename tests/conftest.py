import numpy as np
import pytest


@pytest.fixture
def fcc():
    lattice = np.eye(3)
    positions = np.array(
        [
            [0, 0, 0],
            [0, 0.5, 0.5],
            [0.5, 0, 0.5],
            [0.5, 0.5, 0],
        ]
    )
    numbers = np.array([0, 0, 0, 0])
    magmoms = np.array(
        [
            [0, 0, 1],
            [0, 0, 1],
            [1, 0, 0],
            [1, 0, 0],
        ],
        dtype=np.float_,
    )
    return lattice, positions, numbers, magmoms
