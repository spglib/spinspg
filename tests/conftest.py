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


@pytest.fixture
def layer_triangular_kagome():
    # Kagome lattice with similar to the spin arrangement of the non-collinear antiferromagnetic Mn3Ge and Mn3Sn
    # Family space group: P6/mmm (No. 191)
    # Magnetic space group: Cmm'm'
    # Spin point group: ^{3z}6 / ^{1}m ^{2x}m ^{2xy}m (498)
    a = 5.2
    c = 20
    m = 1.0
    lattice = np.array(
        [
            [-0.5 * a, -np.sqrt(3) / 2 * a, 0],
            [a, 0, 0],
            [0, 0, c],
        ]
    )
    positions = np.array(
        [
            [0, 0, 0],  # Ge(1a)
            [0.5, 0, 0],  # Mn(3f)
            [0, 0.5, 0],  # Mn(3f)
            [0.5, 0.5, 0],  # Mn(3f)
        ]
    )
    numbers = np.array([0, 1, 1, 1])
    magmoms = np.array(
        [
            [0, 0, 0],
            [-0.5 * m, np.sqrt(3) / 2 * m, 0],
            [m, 0, 0],
            [-0.5 * m, -np.sqrt(3) / 2 * m, 0],
        ]
    )
    return lattice, positions, numbers, magmoms


@pytest.fixture
def triangular_kagome():
    # Non-collinear antiferromagnetic Mn3Ge and Mn3Sn
    # https://materialsproject.org/materials/mp-1078873?formula=Mn3Ge
    # Family space group: P6_3/mmc (No. 194)
    a = 5.2
    c = 4.2
    m = 1.0
    x_6h = -1 / 6
    lattice = np.array(
        [
            [-0.5 * a, -np.sqrt(3) / 2 * a, 0],
            [a, 0, 0],
            [0, 0, c],
        ]
    )
    positions = np.array(
        [
            [1 / 3, 2 / 3, 1 / 4],  # Ge(2c)
            [2 / 3, 1 / 3, 3 / 4],  # Ge(2c)
            [x_6h, 2 * x_6h, 1 / 4],  # Mn(6h)
            [-2 * x_6h, -x_6h, 1 / 4],  # Mn(6h)
            [x_6h, -x_6h, 1 / 4],  # Mn(6h)
            [-x_6h, -2 * x_6h, 3 / 4],  # Mn(6h)
            [2 * x_6h, x_6h, 3 / 4],  # Mn(6h)
            [-x_6h, x_6h, 3 / 4],  # Mn(6h)
        ]
    )
    numbers = np.array([0, 0, 1, 1, 1, 1, 1, 1])
    magmoms = np.array(
        [
            [0, 0, 0],  # Ge
            [0, 0, 0],  # Ge
            [-0.5 * m, np.sqrt(3) / 2 * m, 0],  # Mn-A
            [m, 0, 0],  # Mn-B
            [-0.5 * m, -np.sqrt(3) / 2 * m, 0],  # Mn-C
            [-0.5 * m, np.sqrt(3) / 2 * m, 0],  # Mn-A
            [m, 0, 0],  # Mn-B
            [-0.5 * m, -np.sqrt(3) / 2 * m, 0],  # Mn-C
        ]
    )
    return lattice, positions, numbers, magmoms
