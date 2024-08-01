from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.definitions import angular_acceleration_is_angular_speed_derivative as angular_acceleration_def

# Description
## Suppose that the object starts rotating with some angular acceleration. After 5 seconds, it rotates at 180 radians per second.
## Acceleration should be pi/5 radians/sec^2.

Args = namedtuple("Args", ["w0", "w1", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w0 = Quantity(0 * units.radian / units.second)
    w1 = Quantity(pi * units.radian / units.second)
    t = Quantity(5 * units.second)
    return Args(w0=w0, w1=w1, t=t)


def test_basic_acceleration(test_args: Args) -> None:
    result = angular_acceleration_def.calculate_angular_acceleration(test_args.w0, test_args.w1,
        test_args.t)
    assert_equal(result, (pi / 5) * units.radian / units.second**2)


def test_velocity_with_bad_angle(test_args: Args) -> None:
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_acceleration_def.calculate_angular_acceleration(wb, test_args.w1, test_args.t)
    with raises(errors.UnitsError):
        angular_acceleration_def.calculate_angular_acceleration(test_args.w0, wb, test_args.t)


def test_velocity_with_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_acceleration_def.calculate_angular_acceleration(test_args.w0, test_args.w1, tb)
    with raises(TypeError):
        angular_acceleration_def.calculate_angular_acceleration(test_args.w0, test_args.w1, 100)
