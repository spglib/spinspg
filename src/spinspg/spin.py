"""Spin rotation for magnetic moments."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from typing import TYPE_CHECKING

import numpy as np
from spgrep.spinor import get_rotation_angle_and_axis
from typing_extensions import Self

if TYPE_CHECKING:
    from spinspg.utils import NDArrayFloat


class SpinOnlyGroupType(Enum):
    """Type of spin only group."""

    NONMAGNETIC = auto()  # O(3)
    COLLINEAR = auto()  # inf/m
    COPLANAR = auto()  # m
    NONCOPLANAR = auto()  # 1

    def __str__(self) -> str:
        """Return string representation."""
        return self.name


@dataclass
class SpinOnlyGroup:
    """Spin only group.

    This class represents one of nonmagnetic, collinear, coplanar, and non-coplanar spin-only groups.

    Attributes
    ----------
    spin_only_group_type: :class:`spin.SpinOnlyGroupType`
        Nonmagnetic, collinear, coplanar, or non-coplanar
    axis: array or None, (3, )
        For collinear and coplanar, unit vector perpendicular to its rotation axis.
    """

    spin_only_group_type: SpinOnlyGroupType
    axis: NDArrayFloat | None

    def __str__(self) -> str:
        """Return string representation."""
        if self.spin_only_group_type in [SpinOnlyGroupType.COLLINEAR, SpinOnlyGroupType.COPLANAR]:
            return str(self.spin_only_group_type) + f"(axis={self.axis})"
        else:
            return str(self.spin_only_group_type)

    def contain(self, linear: NDArrayFloat, atol: float = 1e-5) -> bool:
        """Return if this spin only group contains ``linear``."""
        if self.spin_only_group_type == SpinOnlyGroupType.NONMAGNETIC:
            return True
        if np.allclose(linear, np.eye(3), atol=atol):
            # Group contains identity
            return True

        if np.linalg.det(linear) > 0:
            sign = 1
            rotation = np.array(linear).copy()
        else:
            sign = -1
            rotation = -np.array(linear)

        theta, axis = get_rotation_angle_and_axis(rotation)

        is_parallel = (self.axis is not None) and np.allclose(
            np.cross(axis, self.axis), 0, atol=atol
        )
        two_fold = np.isclose(theta, np.pi, atol=atol)
        if self.spin_only_group_type == SpinOnlyGroupType.COPLANAR:
            if (sign == -1) and is_parallel and two_fold:
                # ``linear`` is a mirror along self.axis
                return True
            else:
                return False
        elif self.spin_only_group_type == SpinOnlyGroupType.COLLINEAR:
            if (sign == 1) and is_parallel:
                # ``linear`` is a rotation along self.axis
                return True
            elif (sign == -1) and np.isclose(np.inner(axis, self.axis), 0, atol=atol) and two_fold:
                return True
            else:
                return False
        else:
            return False

    @classmethod
    def nonmagnetic(cls) -> Self:
        """Instantiate nonmagnetic spin-only group."""
        return cls(SpinOnlyGroupType.NONMAGNETIC, None)

    @classmethod
    def collinear(cls, axis: NDArrayFloat) -> Self:
        """Instantiate collinear spin-only group with the parallel axis."""
        return cls(SpinOnlyGroupType.COLLINEAR, axis)

    @classmethod
    def coplanar(cls, axis: NDArrayFloat) -> Self:
        """Instantiate coplanar spin-only group with the perpendicular axis."""
        return cls(SpinOnlyGroupType.COPLANAR, axis)

    @classmethod
    def noncoplanar(cls) -> Self:
        """Instantiate noncoplanar spin-only group."""
        return cls(SpinOnlyGroupType.NONCOPLANAR, None)


def get_spin_only_group(magmoms: NDArrayFloat, mag_symprec: float) -> SpinOnlyGroup:
    """Determine spin only group of given spin arrangement.

    Parameters
    ----------
    magmoms : array, (num_sites, 3)
        Magnetic moments in Cartesian coordinates

    Returns
    -------
    spin_only_group: SpinOnlyGroup
    """
    # Nonmagnetic
    if np.max(np.linalg.norm(magmoms, axis=1)) < mag_symprec:
        return SpinOnlyGroup.nonmagnetic()

    moment = np.einsum("ij,ik->jk", magmoms, magmoms, optimize="greedy")  # (3, 3), symmetric
    _, eigvecs = np.linalg.eigh(moment)  # eigenvalues in ascending order

    # Collinear
    parallel_axis = eigvecs[:, -1] / np.linalg.norm(eigvecs[:, -1])
    residual_collinear = (
        magmoms - (magmoms @ parallel_axis)[:, None] * parallel_axis[None, :]
    )  # (N, 3)
    if np.max(2 * np.linalg.norm(residual_collinear, axis=1)) < mag_symprec:
        return SpinOnlyGroup.collinear(parallel_axis)

    # Coplanar
    vertical_axis = eigvecs[:, 0] / np.linalg.norm(eigvecs[:, 0])
    residual_coplanar = (magmoms @ vertical_axis)[:, None] * vertical_axis[None, :]  # (N, 3)
    if np.max(2 * np.linalg.norm(residual_coplanar, axis=1)) < mag_symprec:
        return SpinOnlyGroup.coplanar(vertical_axis)

    # Non-coplanar
    return SpinOnlyGroup.noncoplanar()


def solve_procrustes(A: NDArrayFloat, B: NDArrayFloat) -> NDArrayFloat:
    """Solve orthogonal Procrustes problem.

        argmin_{ R in O(3) } || R A^T - B^T ||_{F}

    Parameters
    ----------
    A: array, (n, 3)
    B: array, (n, 3)

    Returns
    -------
    R: array, (3, 3)
        orthogonal matrix
    """
    u, s, vt = np.linalg.svd(B.T @ A)
    R = u @ vt
    return R
