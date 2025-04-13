"""Utility functions."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray
from typing_extensions import TypeAlias  # for Python<3.10

NDArrayInt: TypeAlias = NDArray[np.int_]
NDArrayFloat: TypeAlias = NDArray[np.float64]


def ndarray2d_to_integer_tuple(array: NDArrayFloat) -> tuple[tuple[Any]]:
    """Convert two-dimensional array to tuple of tuple."""
    array_int = np.around(array).astype(int)
    array_t = tuple(map(tuple, array_int.tolist()))
    return array_t  # type: ignore


def is_integer_array(array: NDArrayFloat, rtol: float = 1e-5, atol: float = 1e-8) -> bool:
    """Return true if all values of ``array`` are almost integers."""
    array_int = np.around(array).astype(int)
    return np.allclose(array_int, array, rtol=rtol, atol=atol)
