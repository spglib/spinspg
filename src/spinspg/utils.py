"""Utility functions."""
from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray
from typing_extensions import TypeAlias  # for Python<3.10

NDArrayInt: TypeAlias = NDArray[np.int_]
NDArrayFloat: TypeAlias = NDArray[np.float_]


def ndarray2d_to_integer_tuple(array: NDArrayFloat) -> tuple[tuple[Any]]:
    """Convert two-dimensional array to tuple of tuple."""
    array_int = np.around(array).astype(int)
    array_t = tuple(map(tuple, array_int.tolist()))
    return array_t  # type: ignore
