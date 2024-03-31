from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic import (
    constant_angular_acceleration_rotation_is_parabolic as angular_displacement_law,)

# Description
## A body is rotating about a fixed axis with constant angular acceleration of -2 rad/s**2.
## Its initial angular velocity was 1 rad/s. In 5 s its angular displacement amounts to
## -20 rad.

Args = namedtuple("Args", "w0 alpha t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w0 = Quantity(1.0 * units.radian / units.second)
    alpha = Quantity(-2.0 * units.radian / units.second**2)
    t = Quantity(5.0 * units.second)
    return Args(w0=w0, alpha=alpha, t=t)


def test_basic_law(test_args: Args) -> None:
    result = angular_displacement_law.calculate_angular_displacement(test_args.w0, test_args.alpha,
        test_args.t)
    assert_equal(result, -20 * units.radian)


def test_bad_angular_velocity(test_args: Args) -> None:
    w0b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_displacement_law.calculate_angular_displacement(w0b, test_args.alpha, test_args.t)
    with raises(TypeError):
        angular_displacement_law.calculate_angular_displacement(100, test_args.alpha, test_args.t)


def test_bad_angular_acceleration(test_args: Args) -> None:
    alpha_b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_displacement_law.calculate_angular_displacement(test_args.w0, alpha_b, test_args.t)
    with raises(TypeError):
        angular_displacement_law.calculate_angular_displacement(test_args.w0, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_displacement_law.calculate_angular_displacement(test_args.w0, test_args.alpha, tb)
    with raises(TypeError):
        angular_displacement_law.calculate_angular_displacement(test_args.w0, test_args.alpha, 100)
