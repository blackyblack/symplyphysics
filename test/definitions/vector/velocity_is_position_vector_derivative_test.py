from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, Quantity, errors
from symplyphysics.definitions.vector import velocity_is_position_vector_derivative as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "dr dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dr = QuantityCoordinateVector([
        3 * units.meter,
        0.5 * units.meter,
        -1 * units.meter,
    ], CARTESIAN)

    dt = Quantity(0.2 * units.second)

    return Args(dr=dr, dt=dt)


def test_velocity_law(test_args: Args) -> None:
    result = law.calculate_velocity(test_args.dr, test_args.dt)

    expected = QuantityCoordinateVector([
        15 * units.meter / units.second,
        2.5 * units.meter / units.second,
        -5 * units.meter / units.second,
    ], CARTESIAN)

    assert_equal_vectors(result, expected)


def test_bad_position(test_args: Args) -> None:
    bad_vector = QuantityCoordinateVector([units.candela, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_velocity(bad_vector, test_args.dt)

    bad_scalar = units.meter
    with raises(ValueError):
        law.calculate_velocity(bad_scalar, test_args.dt)

    with raises(TypeError):
        law.calculate_velocity(100, test_args.dt)
    with raises(TypeError):
        law.calculate_velocity([100], test_args.dt)


def test_bad_time(test_args: Args) -> None:
    bad_scalar = units.candela
    with raises(errors.UnitsError):
        law.calculate_velocity(test_args.dr, bad_scalar)

    with raises(TypeError):
        law.calculate_velocity(test_args.dr, 100)
