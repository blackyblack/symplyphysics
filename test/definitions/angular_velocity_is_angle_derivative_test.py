from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import angular_velocity_is_angle_derivative as angular_velocity_def

# Description
## Assume object starts rotating with some angular velocity. After 5 seconds it rotates to 180 degrees.
## Velocity should be pi/5 radian/sec.

Args = namedtuple("Args", ["a0", "a1", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a0 = Quantity(0 * units.radian)
    a1 = Quantity(pi * units.radian)
    t = Quantity(5 * units.second)
    return Args(a0=a0, a1=a1, t=t)


def test_basic_velocity(test_args: Args) -> None:
    result = angular_velocity_def.calculate_angular_velocity(test_args.a0, test_args.a1,
        test_args.t)
    assert_equal(result, (pi / 5) * units.radian / units.second)


def test_velocity_with_number(test_args: Args) -> None:
    angular_velocity_def.calculate_angular_velocity(100, test_args.a1, test_args.t)
    angular_velocity_def.calculate_angular_velocity(test_args.a0, 100, test_args.t)


def test_velocity_with_bad_angle(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_def.calculate_angular_velocity(ab, test_args.a1, test_args.t)
    with raises(errors.UnitsError):
        angular_velocity_def.calculate_angular_velocity(test_args.a0, ab, test_args.t)


def test_velocity_with_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_def.calculate_angular_velocity(test_args.a0, test_args.a1, tb)
    with raises(TypeError):
        angular_velocity_def.calculate_angular_velocity(test_args.a0, test_args.a1, 100)
