from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.kinematics.vector import acceleration_is_normal_plus_tangential_acceleration as acceleration_law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## A body is moving along a curve. At a certain point in time, its radial acceleration
## with respect to the momentary rotational axis is (3, 1, 1) m/s^2, and its tangential
## acceleration is (-1, 2, 1) m/s^2. Its acceleration should amount to (2, 3, 2) m/s^2.

Args = namedtuple("Args", "a_r a_t a_total")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a_r = QuantityCoordinateVector([
        Quantity(3.0 * units.meter / units.second**2),
        Quantity(1.0 * units.meter / units.second**2),
        Quantity(1.0 * units.meter / units.second**2),
    ], CARTESIAN)
    a_t = QuantityCoordinateVector([
        Quantity(-1.0 * units.meter / units.second**2),
        Quantity(2.0 * units.meter / units.second**2),
        Quantity(1.0 * units.meter / units.second**2),
    ], CARTESIAN)
    a_total = QuantityCoordinateVector([
        Quantity(2.0 * units.meter / units.second**2),
        Quantity(3.0 * units.meter / units.second**2),
        Quantity(2.0 * units.meter / units.second**2),
    ], CARTESIAN)
    return Args(a_r=a_r, a_t=a_t, a_total=a_total)


def test_basic_law(test_args: Args) -> None:
    result = acceleration_law.calculate_acceleration(test_args.a_r, test_args.a_t)
    assert_equal_vectors(result, test_args.a_total)


def test_radial_law(test_args: Args) -> None:
    result = acceleration_law.calculate_radial_acceleration(test_args.a_total, test_args.a_t)
    assert_equal_vectors(result, test_args.a_r)


def test_tangential_law(test_args: Args) -> None:
    result = acceleration_law.calculate_tangential_acceleration(test_args.a_total, test_args.a_r)
    assert_equal_vectors(result, test_args.a_t)


def test_bad_acceleration(test_args: Args) -> None:
    a_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        acceleration_law.calculate_acceleration(a_bad_vector, test_args.a_t)
    with raises(errors.UnitsError):
        acceleration_law.calculate_acceleration(test_args.a_r, a_bad_vector)

    a_bad_scalar = Quantity(1.0 * units.meter / units.second**2)
    with raises(ValueError):
        acceleration_law.calculate_acceleration(a_bad_scalar, test_args.a_t)
    with raises(ValueError):
        acceleration_law.calculate_acceleration(test_args.a_r, a_bad_scalar)

    with raises(TypeError):
        acceleration_law.calculate_acceleration(100, test_args.a_t)
    with raises(TypeError):
        acceleration_law.calculate_acceleration(test_args.a_r, 100)

    a_non_orthogonal = QuantityCoordinateVector([
        Quantity(3.0 * units.meter / units.second**2),
        Quantity(4.0 * units.meter / units.second**2),
        Quantity(1.0 * units.meter / units.second**2),
    ], CARTESIAN)
    with raises(ValueError):
        acceleration_law.calculate_acceleration(a_non_orthogonal, test_args.a_t)
    with raises(ValueError):
        acceleration_law.calculate_acceleration(test_args.a_r, a_non_orthogonal)


def test_bad_radial_acceleration(test_args: Args) -> None:
    a_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        acceleration_law.calculate_radial_acceleration(a_bad_vector, test_args.a_t)
    with raises(errors.UnitsError):
        acceleration_law.calculate_radial_acceleration(test_args.a_total, a_bad_vector)

    a_bad_scalar = Quantity(1.0 * units.meter / units.second**2)
    with raises(ValueError):
        acceleration_law.calculate_radial_acceleration(a_bad_scalar, test_args.a_t)
    with raises(ValueError):
        acceleration_law.calculate_radial_acceleration(test_args.a_total, a_bad_scalar)

    with raises(TypeError):
        acceleration_law.calculate_radial_acceleration(100, test_args.a_t)
    with raises(TypeError):
        acceleration_law.calculate_radial_acceleration(test_args.a_total, 100)


def test_bad_tangential_acceleration(test_args: Args) -> None:
    a_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        acceleration_law.calculate_tangential_acceleration(a_bad_vector, test_args.a_r)
    with raises(errors.UnitsError):
        acceleration_law.calculate_tangential_acceleration(test_args.a_total, a_bad_vector)

    a_bad_scalar = Quantity(1.0 * units.meter / units.second**2)
    with raises(ValueError):
        acceleration_law.calculate_tangential_acceleration(a_bad_scalar, test_args.a_r)
    with raises(ValueError):
        acceleration_law.calculate_tangential_acceleration(test_args.a_total, a_bad_scalar)

    with raises(TypeError):
        acceleration_law.calculate_tangential_acceleration(100, test_args.a_r)
    with raises(TypeError):
        acceleration_law.calculate_tangential_acceleration(test_args.a_total, 100)
