from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
    dimensionless,
)
from symplyphysics.laws.kinematic import (
    angular_position_from_constant_angular_velocity as angular_position_law,)

# Description
## A body is rotating about a fixed axis with a constant angular velocity of 0.5 rad/s.
## Initially at an angle of 3 rad, it rotates to an angular position of 15 rad in 24 s.


@fixture(name="test_args")
def test_args_fixture():
    theta0 = 3.0
    w = Quantity(0.5 / units.second)
    t = Quantity(24.0 * units.second)
    Args = namedtuple("Args", "theta0 w t")
    return Args(theta0=theta0, w=w, t=t)


def test_basic_law(test_args):
    result = angular_position_law.calculate_angular_position(test_args.theta0, test_args.w,
        test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, dimensionless)
    result_value = convert_to(result, units.radian).evalf(3)
    assert_approx(result_value, 15)


def test_bad_angle(test_args):
    theta0b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_position_law.calculate_angular_position(theta0b, test_args.w, test_args.t)


def test_bad_angular_velocity(test_args):
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_position_law.calculate_angular_position(test_args.theta0, wb, test_args.t)
    with raises(TypeError):
        angular_position_law.calculate_angular_position(test_args.theta0, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_position_law.calculate_angular_position(test_args.theta0, test_args.w, tb)
    with raises(TypeError):
        angular_position_law.calculate_angular_position(test_args.theta0, test_args.w, 100)
