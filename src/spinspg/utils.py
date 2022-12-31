"""Utility functions."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing_extensions import TypeAlias  # for Python<3.10

NDArrayInt: TypeAlias = NDArray[np.int_]
NDArrayFloat: TypeAlias = NDArray[np.float_]
