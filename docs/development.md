# Development memo

## Installation

```shell
git clone git@github.com:spglib/spinspg.git
cd spinspg
conda create -y -n spinspg python=3.10 pip
conda activate spinspg
pip install -e ".[dev,docs]"
pre-commit install
```

## Reusing crystal preparation

Issue [#13](https://github.com/spglib/spinspg/issues/13) separates a fixed crystal
from the magnetic moment fields evaluated on it. The public contract is:

- `prepare_spin_symmetry(lattice, positions, numbers, symprec=1e-5,
  angle_tolerance=-1.0, backend="spglib")` returns a `PreparedSpinSymmetry`.
- Preparation owns the derived nonmagnetic symmetry and site permutations. It
  copies the geometry inputs, exposes no mutable geometry, and fixes the ordered
  sites, geometry tolerances, and backend for its lifetime. Changing any of these
  requires another preparation; there is no global cache.
- `prepared.get_spin_symmetry(magmoms, mag_symprec=None)` accepts moments in the
  original site order. It checks their `(num_sites, 3)` shape and defaults the
  magnetic tolerance to the preparation's `symprec`. Each call returns fresh
  results and retains no magnetic moments or previous results.
- Both entry points return the existing four-tuple, including operations in the
  input cell and the same operation order and spin-only representatives.

Before sharing the input-cell expansion, the rutile operation oracle and analytic
FCC operation sets characterize the existing behavior. Validation also covers
both backends, all four spin-only groups, nonprimitive cells, tolerance
boundaries, changed input arrays, and reuse without another backend call.

Performance changes must preserve species matching, the first available matching
site, strict Cartesian distance thresholds, and search early exits. Measure
unprofiled complete searches and separate call-tree profiles before deciding
whether a native kernel warrants an additional build dependency. Keep SVD and
floating-point tolerance decisions unchanged in the initial implementation.

## Compile documents

```shell
sphinx-autobuild docs docs_build
# open localhost:8000 in your browser
```

## Release

```shell
# Confirm the version number via `setuptools-scm`
python -m setuptools_scm

# Update changelog here
vim docs/changelog.md

# Push with tag
git tag <next-version>
git push origin main
git push origin <next-version>
```
