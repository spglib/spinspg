import numpy as np
import pytest

import spinspg.core as core
from spinspg import PreparedSpinSymmetry, get_spin_symmetry, prepare_spin_symmetry


def assert_same_symmetry(actual, expected):
    actual_group, *actual_operations = actual
    expected_group, *expected_operations = expected
    assert actual_group.spin_only_group_type == expected_group.spin_only_group_type
    if expected_group.axis is None:
        assert actual_group.axis is None
    else:
        np.testing.assert_allclose(actual_group.axis, expected_group.axis)
    # Compare complete arrays in order, including input-cell translations.
    for actual_array, expected_array in zip(actual_operations, expected_operations):
        np.testing.assert_allclose(actual_array, expected_array, atol=1e-12)


@pytest.mark.parametrize("backend", ["spglib", "moyopy"])
@pytest.mark.parametrize("fixture", ["rutile", "fcc", "layer_triangular_kagome"])
def test_prepared_fields(request, fixture, backend):
    lattice, positions, numbers, moments = request.getfixturevalue(fixture)
    prepared = prepare_spin_symmetry(lattice, positions, numbers, symprec=1e-3, backend=backend)
    assert isinstance(prepared, PreparedSpinSymmetry)
    fields = [
        moments,
        np.zeros_like(moments),
        moments @ np.array([[0.8, -0.6, 0], [0.6, 0.8, 0], [0, 0, 1]]),
        np.random.default_rng(13).normal(size=moments.shape),
    ]
    for field in fields:
        original = field.copy()
        assert_same_symmetry(
            prepared.get_spin_symmetry(field),
            get_spin_symmetry(lattice, positions, numbers, field, symprec=1e-3, backend=backend),
        )
        np.testing.assert_array_equal(field, original)


@pytest.mark.parametrize("backend", ["spglib", "moyopy"])
def test_prepared_ownership_and_reuse(rutile, backend, monkeypatch):
    lattice, positions, numbers, moments = rutile
    original = [array.copy() for array in rutile]
    expected = get_spin_symmetry(*original, backend=backend)
    expected_changed_moments = get_spin_symmetry(
        *original[:3], np.zeros_like(moments), backend=backend
    )
    prepare = core.get_symmetry_with_cell
    calls = []

    def record_preparation(*args, **kwargs):
        for copied, source in zip(args[:3], rutile[:3]):
            assert not np.shares_memory(copied, source)
        calls.append((args, kwargs))
        return prepare(*args, **kwargs)

    monkeypatch.setattr(core, "get_symmetry_with_cell", record_preparation)
    prepared = prepare_spin_symmetry(lattice, positions, numbers, backend=backend)
    lattice[:] = 0
    positions[:] = 0
    numbers[:] = 99
    result = prepared.get_spin_symmetry(moments)
    assert_same_symmetry(result, expected)
    group, *operations = result
    group.axis[:] = 0
    for array in operations:
        array[:] = 0
    moments[:] = 0
    assert_same_symmetry(prepared.get_spin_symmetry(moments), expected_changed_moments)
    assert_same_symmetry(prepared.get_spin_symmetry(original[3]), expected)
    assert len(calls) == 1


@pytest.mark.parametrize("backend", ["spglib", "moyopy"])
def test_prepared_tolerances(fcc, backend):
    lattice, positions, numbers, _ = fcc
    moments = np.zeros((len(positions), 3))
    moments[:, 2] = 1
    moments[positions[:, 0] == 0, 2] += 1e-3
    prepared = prepare_spin_symmetry(
        lattice, positions, numbers, symprec=1.01e-3, angle_tolerance=0.1, backend=backend
    )
    for mag_symprec, order in [(None, 192), (0.99e-3, 32)]:
        actual = prepared.get_spin_symmetry(moments, mag_symprec)
        expected = get_spin_symmetry(
            lattice,
            positions,
            numbers,
            moments,
            symprec=1.01e-3,
            angle_tolerance=0.1,
            mag_symprec=mag_symprec,
            backend=backend,
        )
        assert len(actual[1]) == order
        assert_same_symmetry(actual, expected)
    assert prepared.symprec == 1.01e-3
    assert prepared.angle_tolerance == 0.1
    assert prepared.backend == backend
    assert prepared.num_sites == len(positions)
    for name in ["symprec", "angle_tolerance", "backend", "num_sites"]:
        with pytest.raises(AttributeError):
            setattr(prepared, name, None)


@pytest.mark.parametrize("shape", [(6,), (3, 6), (5, 3), (7, 3), (6, 3, 1)])
def test_prepared_moment_shape(rutile, shape):
    prepared = prepare_spin_symmetry(*rutile[:3])
    with pytest.raises(ValueError, match=r"magmoms must have shape \(6, 3\)"):
        prepared.get_spin_symmetry(np.zeros(shape))


@pytest.mark.parametrize("backend", ["spglib", "moyopy"])
def test_independent_geometry_tolerances(backend):
    lattice = np.diag([1.0, 1.0, 1.0002])
    positions = np.zeros((1, 3))
    numbers = np.array([0])
    moments = np.zeros((1, 3))
    tetragonal = prepare_spin_symmetry(lattice, positions, numbers, symprec=1e-5, backend=backend)
    cubic = prepare_spin_symmetry(lattice, positions, numbers, symprec=1e-3, backend=backend)
    for prepared, order in [(tetragonal, 16), (cubic, 48), (tetragonal, 16)]:
        _, rotations, _, _ = prepared.get_spin_symmetry(moments)
        assert len(rotations) == order
