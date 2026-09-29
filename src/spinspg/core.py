"""Core APIs."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from spinspg.group import (
    SYMMETRY_FINDER_BACKEND,
    get_primitive_spin_symmetry,
    get_symmetry_with_cell,
)

if TYPE_CHECKING:
    from spinspg.group import SpinSpaceGroup
    from spinspg.spin import SpinOnlyGroup
    from spinspg.utils import NDArrayFloat, NDArrayInt


def get_spin_symmetry(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    magmoms: NDArrayFloat,
    symprec: float = 1e-5,
    mag_symprec: float | None = None,
    angle_tolerance: float = -1.0,
    backend: SYMMETRY_FINDER_BACKEND = "spglib",
) -> tuple[SpinOnlyGroup, NDArrayInt, NDArrayFloat, NDArrayFloat]:
    """Return spin symmetry operations of a given spin arrangement.

    See :ref:`Spglib's document <spglib:py_variables_crystal_structure>` for how to specify the spin arrangement by ``lattice``, ``positions``, ``numbers``, and ``magmoms`` in details.
    With the returned spin symmetry operations, the ``i``-th spin symmetry operation maps point coordinates ``x`` to ``rotations[i] @ x + translations[i]`` and magnetic moment ``m`` to ``spin_rotations[i] @ m``.

    Parameters
    ----------
    lattice: array, (3, 3)
        ``lattice[i, :]`` is the ``i``-th basis vector of a lattice
    positions: array, (num_sites, 3)
        ``positions[i, :]`` is a fractional coordinates of the ``i``-th site w.r.t. ``lattice``.
    numbers: array[int], (num_sites, )
        ``numbers[i]`` specifies a specie at the ``i``-th site.
    magmoms: array, (num_sites, 3)
        ``magmoms[i, :]`` is a magnetic moments at the ``i``-th site in Cartesian coordinates.
    symprec: float, default=1e-5
        See :ref:`spglib:variables_symprec`.
    mag_symprec: float | None
        See :ref:`spglib:variables_mag_symprec`.
    angle_tolerance: float, default=-1
        See :ref:`spglib:variables_angle_tolerance`.
    backend: {"spglib", "moyopy"}, default="spglib"
        Backend for finding the nonmagnetic symmetry.

    Returns
    -------
    spin_only_group: :class:`spin.SpinOnlyGroup`
        Spin-only group, which preserve the given spin arrangement with identity spatial operation.
    rotations: array[int], (num_sym, 3, 3)
        Rotation parts of spin symmetry operations w.r.t. ``lattice``.
    translations: array, (num_sym, 3)
        Translation parts of spin symmetry operations w.r.t. ``lattice``.
    spin_rotations: array, (num_sym, 3, 3)
        Spin rotation parts of spin symmetry operations in Cartesian coordinates.

    """
    return prepare_spin_symmetry(
        lattice, positions, numbers, symprec, angle_tolerance, backend
    ).get_spin_symmetry(magmoms, mag_symprec)


class PreparedSpinSymmetry:
    """Reusable nonmagnetic symmetry for a fixed, ordered crystal.

    Construct with :func:`prepare_spin_symmetry`. Geometry inputs are copied
    during preparation, and only private derived geometry is retained. Mutating
    the original inputs cannot affect this object. Prepare another object when
    the lattice, ordered positions, species, geometry tolerances, or backend
    changes. Magnetic moments and returned results are never cached.
    """

    def __init__(
        self,
        lattice: NDArrayFloat,
        positions: NDArrayFloat,
        numbers: NDArrayInt,
        symprec: float = 1e-5,
        angle_tolerance: float = -1.0,
        backend: SYMMETRY_FINDER_BACKEND = "spglib",
    ) -> None:
        lattice = np.array(lattice, dtype=np.float64, copy=True)
        positions = np.array(positions, dtype=np.float64, copy=True)
        numbers = np.array(numbers, dtype=np.int_, copy=True)
        self._num_sites = len(positions)
        self._symprec = symprec
        self._angle_tolerance = angle_tolerance
        self._backend = backend
        self._nonmagnetic_symmetry = get_symmetry_with_cell(
            lattice, positions, numbers, symprec, angle_tolerance, backend=backend
        )

    @property
    def num_sites(self) -> int:
        """Number of sites in the input cell, including nonmagnetic sites."""
        return self._num_sites

    @property
    def symprec(self) -> float:
        """Fixed geometry tolerance and default magnetic tolerance."""
        return self._symprec

    @property
    def angle_tolerance(self) -> float:
        """Fixed angle tolerance in degrees; negative means backend default."""
        return self._angle_tolerance

    @property
    def backend(self) -> SYMMETRY_FINDER_BACKEND:
        """Backend used to prepare the nonmagnetic symmetry."""
        return self._backend

    def get_spin_symmetry(
        self, magmoms: NDArrayFloat, mag_symprec: float | None = None
    ) -> tuple[SpinOnlyGroup, NDArrayInt, NDArrayFloat, NDArrayFloat]:
        """Evaluate a moment field without repeating geometry preparation.

        Parameters
        ----------
        magmoms: array, (num_sites, 3)
            Cartesian magnetic moments in the original input site order.
        mag_symprec: float | None
            Magnetic tolerance. Defaults to this preparation's ``symprec``.

        Returns
        -------
        tuple
            The same spin-only group, rotations, translations, and spin rotations
            as :func:`get_spin_symmetry`, in the original input cell and operation
            order. Each call returns independent results.

        Raises
        ------
        ValueError
            If the moment array does not have shape ``(num_sites, 3)``.
        """
        magmoms = np.asarray(magmoms, dtype=np.float64)
        if magmoms.shape != (self.num_sites, 3):
            raise ValueError(f"magmoms must have shape ({self.num_sites}, 3)")
        ssg = get_primitive_spin_symmetry(
            self._nonmagnetic_symmetry,
            magmoms,
            mag_symprec=self.symprec if mag_symprec is None else mag_symprec,
        )
        return _expand_spin_symmetry(ssg)


def prepare_spin_symmetry(
    lattice: NDArrayFloat,
    positions: NDArrayFloat,
    numbers: NDArrayInt,
    symprec: float = 1e-5,
    angle_tolerance: float = -1.0,
    backend: SYMMETRY_FINDER_BACKEND = "spglib",
) -> PreparedSpinSymmetry:
    """Prepare a crystal once for evaluating multiple magnetic moment fields.

    Parameters
    ----------
    lattice: array, (3, 3)
        Lattice basis vectors as rows, in Cartesian coordinates.
    positions: array, (num_sites, 3)
        Fractional site coordinates with respect to ``lattice``.
    numbers: array[int], (num_sites,)
        Species identifiers in the same site order as ``positions``.
    symprec: float, default=1e-5
        Geometry tolerance; also the default magnetic tolerance for evaluations.
    angle_tolerance: float, default=-1
        Angle tolerance in degrees. A negative value uses the backend default.
    backend: {"spglib", "moyopy"}, default="spglib"
        Backend for finding the nonmagnetic symmetry.

    Returns
    -------
    PreparedSpinSymmetry
        An owned preparation for this geometry, site order, and tolerance/backend
        choice. Call its :meth:`PreparedSpinSymmetry.get_spin_symmetry` with each
        moment field. Changes to the original arrays do not alter the preparation.
    """
    return PreparedSpinSymmetry(lattice, positions, numbers, symprec, angle_tolerance, backend)


def _expand_spin_symmetry(
    ssg: SpinSpaceGroup,
) -> tuple[SpinOnlyGroup, NDArrayInt, NDArrayFloat, NDArrayFloat]:
    """Expand primitive spin symmetry operations into the input cell."""
    spin_only_group = ssg.spin_only_group
    tmat = ssg.transformation
    invtmat = np.linalg.inv(tmat)
    rotations = []
    translations = []
    spin_rotations = []

    # Products of "translations in cell", "nontrivial spin translation group's coset", and "nontrivial spin space group's coset"
    for ops in ssg.nontrivial_coset:
        # Transform to primitive to input cell
        new_rotation = np.around(invtmat @ ops.rotation @ tmat).astype(np.int_)
        for ops_st in ssg.spin_translation_coset:
            for centering in ssg.prim_centerings:
                new_translation = np.remainder(
                    invtmat @ (ops.translation + ops_st.translation + centering), 1
                )
                rotations.append(new_rotation)
                translations.append(new_translation)
                spin_rotations.append(ops_st.spin_rotation @ ops.spin_rotation)

    return spin_only_group, np.array(rotations), np.array(translations), np.array(spin_rotations)
