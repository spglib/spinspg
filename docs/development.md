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

### Benchmark and native-kernel decision

Run the public README rutile example with changing spin frames:

```shell
uv run python scripts/benchmark_spin_symmetry.py \
  --output .cache/spin-symmetry.json --profile-dir .cache/spin-symmetry-profiles
```

The script limits BLAS threads before importing NumPy, times one representative
search, and reports the median of five unprofiled 256-frame passes. It also times
256 distinct cells obtained by uniformly scaling the rutile lattice. Preparation
is timed separately from prepared evaluations. Optional cProfile passes run
after the unprofiled measurements. Use `--backend moyopy` to measure that backend.
The script also works with versions that only provide the one-shot API.

Results measured on macOS ARM64, Python 3.13.5, NumPy 2.3.3, spglib 2.6.0,
and moyopy 0.7.1, using the spglib backend (seconds per 256-frame pass):

| Implementation | Repeated one-shot calls | Cold distinct cells | Prepared evaluations |
| --- | ---: | ---: | ---: |
| Original (`e0bd13f`) | 2.119 | 2.129 | — |
| Reuse API, before Python optimizations (`5603229`) | 2.031 | 2.111 | 0.869 |
| Reuse API and Python optimizations | 1.566 | 1.591 | 0.476 |

The final preparation took about 4 ms. Including that cost, the repeated prepared
workload was about 4.4 times faster than the original one-shot loop. Cold searches
improved by about 25%. These are local measurements on this small public example,
not performance guarantees or timing assertions for CI.

The retained changes use indexed permutation composition, species-grouped site
candidates, exact integer-rotation identity checks, batched rotation-integrality
checks, and direct collinear comparisons with the original absolute and relative
tolerances. Batching distances over candidate sites was also measured, but made
site matching slower on this example; the scalar Cartesian norm was retained.

The final prepared profile took 0.715 s, of which permutation composition used
0.013 s (about 2%) and SVD used 0.149 s (about 21%). Even removing SVD's entire
cost would cap the prepared speedup at about 1.26 times, before native-call
overhead. Moving permutation composition alone has little remaining benefit.
The remaining search time is spread across residual checks, spin purification,
group construction, and input-cell expansion.

No PyO3 dependency is added on this evidence. If larger-cell workloads motivate
native work, first compare a fused permutation/residual/comparison kernel against
complete prepared and cold searches using the same scientific regression tests.
Keep the existing SVD and search order until those measurements justify moving
more of the search.

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
