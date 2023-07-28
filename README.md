# spinspg
[![testing](https://github.com/spglib/spinspg/actions/workflows/testing.yml/badge.svg)](https://github.com/spglib/spinspg/actions/workflows/testing.yml)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


Python package for detecting spin space group on top of spglib

- Document(latest): <https://spinspg.readthedocs.io/en/latest/>
- GitHub: <https://github.com/spglib/spinspg>

## Features

- Find spin symmetry operations from spin arrangements

## Usage

```python
import numpy as np
from spinspg import get_spin_symmetry

# Antiferromagnetic rutile structure
a = 4.87
c = 3.31
x_4f = 0.695169
lattice = np.diag([a, a, c])
positions = np.array([  # Fractional coordinates
    [0, 0, 0],  # Mn(2a)
    [0.5, 0.5, 0.5],  # Mn(2a)
    [x_4f, x_4f, 0],  # F(4f)
    [-x_4f, -x_4f, 0],  # F(4f)
    [-x_4f + 0.5, x_4f + 0.5, 0.5],  # F(4f)
    [x_4f + 0.5, -x_4f + 0.5, 0.5],  # F(4f)
])
numbers = np.array([0, 0, 1, 1, 1, 1])
magmoms = np.array([  # In Cartesian coordinates
    [0, 0, 2.5],
    [0, 0, -2.5],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
])

# Find spin symmetry operations
sog, rotations, translations, spin_rotations = get_spin_symmetry(lattice, positions, numbers, magmoms)

print(f"Spin-only group: {sog}")  # COLLINEAR(axis=[0. 0. 1.])

# Some operations have nontrivial spin rotations
idx = 2
print(f"Rotation ({idx})\n{rotations[idx]}")
print(f"Translation ({idx})\n{translations[idx]}")
print(f"Spin rotation ({idx})\n{spin_rotations[idx]}")  # -> diag([1, 1, -1])
```

See [API documents](docs/api/api.md) for more details.

## Installation

```shell
git clone git@github.com:spglib/spinspg.git
cd spinspg
pip install .
```

## How to cite spinspg

If you use spinspg in your research, please cite [Spglib](https://spglib.readthedocs.io/en/latest/) and the following paper:

```
@misc{spinspg,
    author = {Kohei Shinohara and Atsushi Togo and Hikaru Watanabe and Takuya Nomoto and Isao Tanaka and Ryotaro Arita},
    title = {Algorithm for spin symmetry operation search},
    year = {2023},
    eprint = {arXiv:2307.12228},
    howpublished = {\url{https://arxiv.org/abs/2307.12228}},
}
```

## Change log

See the [change log](docs/changelog.md) for recent changes.

## License

Spinspg is distributed under a BSD 3-clause license.
