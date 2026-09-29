# spinspg
[![testing](https://github.com/spglib/spinspg/actions/workflows/testing.yml/badge.svg)](https://github.com/spglib/spinspg/actions/workflows/testing.yml)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![PyPI version](https://badge.fury.io/py/spinspg.svg)](https://badge.fury.io/py/spinspg)
![PyPI - Downloads](https://img.shields.io/pypi/dm/spinspg)

`spinspg` is a Python package for detecting spin space group on top of `spglib`

- Document(latest): <https://spinspg.readthedocs.io/en/latest/>
- GitHub: <https://github.com/spglib/spinspg>
- PyPI: <https://pypi.org/project/spinspg/>

## Features

- Find spin symmetry operations from spin arrangements

## Usage

{func}`spinspg.get_spin_symmetry` returns spin symmetry operations of a given spin arrangement, analogous to Spglib's {ref}`spglib:py_get_magnetic_symmetry` for magnetic symmetry operations.
For comprehensive output details, refer to [API documents](docs/api/api.md).

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

For multiple moment fields on the same crystal, prepare its nonmagnetic symmetry
once:

```python
from spinspg import prepare_spin_symmetry

prepared = prepare_spin_symmetry(lattice, positions, numbers, symprec=1e-3)
for angle in np.linspace(0, 2 * np.pi, 256, endpoint=False):
    moments = np.zeros((len(positions), 3))
    moments[0] = 2.5 * np.array([np.sin(angle), 0, np.cos(angle)])
    moments[1] = -moments[0]
    sog, rotations, translations, spin_rotations = prepared.get_spin_symmetry(moments)
```

Moments must follow the original site order. Results have the same format and
input-cell coordinates as `get_spin_symmetry`. The preparation owns its geometry;
changing the original arrays does not change it. Prepare a new object for another
lattice, ordered positions, species, geometry tolerance, or backend. The magnetic
tolerance can vary per evaluation via `mag_symprec` and otherwise defaults to the
preparation's `symprec`.

## Installation

```shell
pip install spinspg
```

## How to cite spinspg

If you use `spinspg` in your research, please cite both [Spglib](https://spglib.readthedocs.io/en/latest/) and the subsequent paper:

```
@article{spinspg,
    author = "Shinohara, Kohei and Togo, Atsushi and Watanabe, Hikaru and Nomoto, Takuya and Tanaka, Isao and Arita, Ryotaro",
    title = "{Algorithm for spin symmetry operation search}",
    journal = "Acta Cryst. A",
    year = "2024",
    volume = "80",
    number = "1",
    pages = "94--103",
    month = "Jan",
    doi = {10.1107/S2053273323009257},
    url = {https://doi.org/10.1107/S2053273323009257},
}

```

## Change log

See the [change log](docs/changelog.md) for recent changes.

## License

`spinspg` is released under a BSD 3-clause license.
