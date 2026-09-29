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

## Symmetry-search performance

Run the public README rutile example with changing spin frames:

```shell
uv run python scripts/benchmark_spin_symmetry.py \
  --output .cache/spin-symmetry.json --profile-dir .cache/spin-symmetry-profiles
```

The script limits BLAS threads before importing NumPy, times one representative
search, and reports the median of five unprofiled 256-frame passes. It also times
256 distinct cells obtained by uniformly scaling the rutile lattice. Both
workloads call `get_spin_symmetry` for each frame, including geometry preparation.
Optional cProfile passes run after the unprofiled measurements. Use
`--backend moyopy` to measure that backend.

The retained changes use indexed permutation composition, species-grouped site
candidates, exact integer-rotation identity checks, batched rotation-integrality
checks, and direct collinear comparisons with the original absolute and relative
tolerances. Batching distances over candidate sites was also measured, but made
site matching slower on this example; the scalar Cartesian norm was retained.

The rutile operation oracle and analytic FCC operation sets cover both backends,
all four spin-only groups, and nonprimitive cells. Boundary tests preserve strict
Cartesian site matching and both absolute and relative spin tolerances.

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
