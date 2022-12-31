from __future__ import annotations

from dataclasses import dataclass

from spinspg.utils import NDArrayFloat, NDArrayInt


@dataclass
class SpinSymmetry:
    """Spin symmetry operation.

    Attributes
    ----------
    rotations: array[int], (3, 3)
    translation: array, (3, )
    spin_rotation: array, (3, 3)
    """

    rotation: NDArrayInt
    translation: NDArrayFloat
    spin_rotation: NDArrayFloat
