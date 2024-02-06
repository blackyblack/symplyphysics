from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_approx, errors, units, convert_to, Quantity, SI)
from symplyphysics.definitions import angular_acceleration_is_angular_velocity_derivative as angular_acceleration_def

# Description
## Suppose that the object starts rotating with some angular acceleration. After 5 seconds, it rotates at 180 radians per second.
## Acceleration should be pi/5 radians/sec^2.


@fixture(name="test_args")
def test_args_fixture():
    w0 = Quantity(0 * units.radian / units.second)
    w1 = Quantity(pi * units.radian / units.second)
    t = Quantity(5 * units.second)
    Args = namedtuple("Args", ["w0", "w1", "t"])
    return Args(w0=w0, w1=w1, t=t)


def test_basic_acceleration(test_args):
    result = angular_acceleration_def.calculate_angular_acceleration(test_args.w0, test_args.w1,
        test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / (units.time**2))
    result_velocity = convert_to(result, angular_acceleration_def.definition_units_SI).evalf(4)
    assert_approx(result_velocity, (pi / 5).evalf(4))


def test_velocity_with_bad_angle(test_args):
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_acceleration_def.calculate_angular_acceleration(wb, test_args.w1, test_args.t)
    with raises(errors.UnitsError):
        angular_acceleration_def.calculate_angular_acceleration(test_args.w0, wb, test_args.t)


def test_velocity_with_bad_time(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_acceleration_def.calculate_angular_acceleration(test_args.w0, test_args.w1, tb)
    with raises(TypeError):
        angular_acceleration_def.calculate_angular_acceleration(test_args.w0, test_args.w1, 100)
