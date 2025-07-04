from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, errors, units, Quantity
from symplyphysics.laws.dynamics.vector import instantaneous_power_is_force_dot_velocity as power_law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector

# Description
## A force is acting on an object, and at some time, the force vector is (1, 1, -1) N, and the object
## moves with the velocity vector equal to (2, 0, -1) m/s at that time. The power of the force exerted
## on the object is 3 W.

Args = namedtuple("Args", "f v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = QuantityCoordinateVector([
        Quantity(1.0 * units.newton),
        Quantity(1.0 * units.newton),
        Quantity(-1.0 * units.newton),
    ], CARTESIAN)
    v = QuantityCoordinateVector([
        Quantity(2.0 * units.meter / units.second),
        Quantity(0.0 * units.meter / units.second),
        Quantity(-1.0 * units.meter / units.second),
    ], CARTESIAN)
    return Args(f=f, v=v)


def test_basic_law(test_args: Args) -> None:
    result = power_law.calculate_power(test_args.f, test_args.v)
    assert_equal(result, 3 * units.watt)


def test_bad_force(test_args: Args) -> None:
    f_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        power_law.calculate_power(f_bad_vector, test_args.v)

    f_scalar = Quantity(1.0 * units.newton)
    with raises(ValueError):
        power_law.calculate_power(f_scalar, test_args.v)

    with raises(TypeError):
        power_law.calculate_power(100, test_args.v)
    with raises(TypeError):
        power_law.calculate_power([100], test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    v_bad_vector = QuantityCoordinateVector([
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
        Quantity(1.0 * units.second),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        power_law.calculate_power(test_args.f, v_bad_vector)

    v_scalar = Quantity(1.0 * units.meter / units.second)
    with raises(ValueError):
        power_law.calculate_power(test_args.f, v_scalar)

    with raises(TypeError):
        power_law.calculate_power(test_args.f, 100)
    with raises(TypeError):
        power_law.calculate_power(test_args.f, [100])
