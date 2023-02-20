# Change Log

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
