from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, Quantity, errors
from symplyphysics.definitions.vector import acceleration_is_velocity_derivative as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "v0 v1 dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    unit = units.meter / units.second
    v0 = QuantityCoordinateVector([1 * unit, 2 * unit, 1 * unit], CARTESIAN)
    v1 = QuantityCoordinateVector([2 * unit, -1 * unit, 3 * unit], CARTESIAN)
    dt = Quantity(0.5 * units.second)

    return Args(v0=v0, v1=v1, dt=dt)


def test_law(test_args: Args) -> None:
    result = law.calculate_acceleration(test_args.v0, test_args.v1, test_args.dt)

    unit = units.meter / units.second**2
    expected = QuantityCoordinateVector([2 * unit, -6 * unit, 4 * unit], CARTESIAN)

    assert_equal_vectors(result, expected)


def test_bad_velocity(test_args: Args) -> None:
    v_bad_vector = QuantityCoordinateVector([1 * units.coulomb, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_acceleration(v_bad_vector, test_args.v1, test_args.dt)
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.v0, v_bad_vector, test_args.dt)

    v_bad_scalar = Quantity(1 * units.meter / units.second)
    with raises(ValueError):
        law.calculate_acceleration(v_bad_scalar, test_args.v1, test_args.dt)
    with raises(ValueError):
        law.calculate_acceleration(test_args.v0, v_bad_scalar, test_args.dt)

    with raises(TypeError):
        law.calculate_acceleration(100, test_args.v1, test_args.dt)
    with raises(TypeError):
        law.calculate_acceleration([100], test_args.v1, test_args.dt)
    with raises(TypeError):
        law.calculate_acceleration(test_args.v0, 100, test_args.dt)
    with raises(TypeError):
        law.calculate_acceleration(test_args.v0, [100], test_args.dt)


def test_bad_time(test_args: Args) -> None:
    t_bad = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_acceleration(test_args.v0, test_args.v1, t_bad)
    with raises(TypeError):
        law.calculate_acceleration(test_args.v0, test_args.v1, 100)
    with raises(TypeError):
        law.calculate_acceleration(test_args.v0, test_args.v1, [100])
