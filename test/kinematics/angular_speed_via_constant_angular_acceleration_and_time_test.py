from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)
from symplyphysics.laws.kinematics import angular_speed_via_constant_angular_acceleration_and_time as angular_velocity_law

# Description
## A body is rotating about a fixed axis with angular acceleration of 3 rad/s^2. Initially its
## angular velocity is -3 rad/s. In 3 seconds, its angular velocity amounts to 6 rad/s.

Args = namedtuple("Args", "w0 alpha t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w0 = Quantity(-3.0 * units.radian / units.second)
    alpha = Quantity(3.0 * units.radian / units.second**2)
    t = Quantity(3.0 * units.second)
    return Args(w0=w0, alpha=alpha, t=t)


def test_basic_law(test_args: Args) -> None:
    result = angular_velocity_law.calculate_angular_velocity(test_args.w0, test_args.alpha,
        test_args.t)
    assert_equal(result, 6 * units.radian / units.second)


def test_bad_angular_velocity(test_args: Args) -> None:
    w0b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_law.calculate_angular_velocity(w0b, test_args.alpha, test_args.t)
    with raises(TypeError):
        angular_velocity_law.calculate_angular_velocity(100, test_args.alpha, test_args.t)


def test_bad_angular_acceleration(test_args: Args) -> None:
    alpha_b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_law.calculate_angular_velocity(test_args.w0, alpha_b, test_args.t)
    with raises(TypeError):
        angular_velocity_law.calculate_angular_velocity(test_args.w0, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_law.calculate_angular_velocity(test_args.w0, test_args.alpha, tb)
    with raises(TypeError):
        angular_velocity_law.calculate_angular_velocity(test_args.w0, test_args.alpha, 100)
