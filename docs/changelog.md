# Change Log

## v0.3.0 (28 Oct. 2025)

- Add moyopy backend support for symmetry search

## v0.2.0 (5 Oct. 2025)

- Drop support for Python 3.9
- Fix collinear case with almost-zero magnetic moments
- Eliminate spin rotations from spin-only group


## v0.1.5 (30 Oct. 2024)

- Drop support for Python 3.8
- Fix typo in scripts/generate_spin_point_group_table.py by @janpriessnitz

## v0.1.2 (28 Jul. 2023)

- First release version

## v0.1.1 (20 Feb. 2023)
- Fix spin symmetry search and add table of nontrivial spin point-group types [[#2]](https://github.com/spglib/spinspg/pull/2)
    - Tabulate 598 spin point group types
    - Fix eigenvectors for spin-only group search
    - Fix to correctly search complement in nontrivial spin space group search
    - Add note for spin point group and glossary

## v0.1.0 (18 Jan. 2023)
- Add functionality to find spin symmetry operations
    - {func}`spinspg.get_spin_symmetry` returns spin symmetry operations from a given spin arrangement
    - {class}`spinspg.spin.SpinOnlyGroup` represents a spin-only group.
