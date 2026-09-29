"""Import top APIs and version."""

from importlib.metadata import PackageNotFoundError, version

from spinspg.core import (  # noqa: F401
    PreparedSpinSymmetry,
    get_spin_symmetry,
    prepare_spin_symmetry,
)

# https://github.com/pypa/setuptools_scm/#retrieving-package-version-at-runtime
try:
    __version__ = version("spinspg")
except PackageNotFoundError:
    # package is not installed
    pass
