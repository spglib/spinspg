from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto

import numpy as np

from spinspg.utils import NDArrayFloat


class SpinOnlyGroupType(Enum):
    NONMAGNETIC = auto()  # O(3)
    COLLINEAR = auto()  # inf/m
    COPLANAR = auto()  # m
    NONCOPLANAR = auto()  # 1


@dataclass
class SpinOnlyGroup:
    """Spin only group.

    Attributes
    ----------
    spin_only_group_type: Nonmagnetic, collinear, coplanar, or non-coplanar
    axis: array or None, (3, )
        For collinear and coplanar, unit vector perpendicular to its rotation axis.
    """

    spin_only_group_type: SpinOnlyGroupType
    axis: NDArrayFloat | None

    @classmethod
    def nonmagnetic(cls) -> SpinOnlyGroup:
        return SpinOnlyGroup(SpinOnlyGroupType.NONMAGNETIC, None)

    @classmethod
    def collinear(cls, axis: NDArrayFloat) -> SpinOnlyGroup:
        return SpinOnlyGroup(SpinOnlyGroupType.COLLINEAR, axis)

    @classmethod
    def coplanar(cls, axis: NDArrayFloat) -> SpinOnlyGroup:
        return SpinOnlyGroup(SpinOnlyGroupType.COPLANAR, axis)

    @classmethod
    def noncoplanar(cls) -> SpinOnlyGroup:
        return SpinOnlyGroup(SpinOnlyGroupType.NONCOPLANAR, None)


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
    parallel_axis = eigvecs[-1, :] / np.linalg.norm(eigvecs[-1, :])
    residual_collinear = (
        magmoms - (magmoms @ parallel_axis)[:, None] * parallel_axis[None, :]
    )  # (N, 3)
    if np.max(2 * np.linalg.norm(residual_collinear, axis=1)) < mag_symprec:
        return SpinOnlyGroup.collinear(parallel_axis)

    # Coplanar
    vertical_axis = eigvecs[0, :] / np.linalg.norm(eigvecs[0, :])
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
