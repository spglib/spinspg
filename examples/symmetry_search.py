import numpy as np

from spinspg import get_spin_symmetry

# Antiferromagnetic rutile structure
a = 4.87
c = 3.31
x_4f = 0.695169
lattice = np.diag([a, a, c])
positions = np.array(
    [  # Fractional coordinates
        [0, 0, 0],  # Mn(2a)
        [0.5, 0.5, 0.5],  # Mn(2a)
        [x_4f, x_4f, 0],  # F(4f)
        [-x_4f, -x_4f, 0],  # F(4f)
        [-x_4f + 0.5, x_4f + 0.5, 0.5],  # F(4f)
        [x_4f + 0.5, -x_4f + 0.5, 0.5],  # F(4f)
    ]
)
numbers = np.array([0, 0, 1, 1, 1, 1])
magmoms = np.array(
    [  # In Cartesian coordinates
        [0, 0, 2.5],
        [0, 0, -2.5],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
    ]
)

# Find spin symmetry operations
sog, rotations, translations, spin_rotations = get_spin_symmetry(
    lattice, positions, numbers, magmoms
)

print(f"Spin-only group: {sog}")  # COLLINEAR(axis=[0. 0. 1.])

# Some operations have nontrivial spin rotations
idx = 2
print(f"Rotation ({idx})\n{rotations[idx]}")
print(f"Translation ({idx})\n{translations[idx]}")
print(f"Spin rotation ({idx})\n{spin_rotations[idx]}")  # -> diag([1, 1, -1])
