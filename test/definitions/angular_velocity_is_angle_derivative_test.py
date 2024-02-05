from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import angular_velocity_is_angle_derivative as angular_velocity_def

# Description
## Assume object starts rotating with some angular velocity. After 5 seconds it rotates to 180 degrees.
## Velocity should be pi/5 radian/sec.


@fixture(name="test_args")
def test_args_fixture():
    a0 = Quantity(0 * units.radian)
    a1 = Quantity(pi * units.radian)
    t = Quantity(5 * units.second)
    Args = namedtuple("Args", ["a0", "a1", "t"])
    return Args(a0=a0, a1=a1, t=t)


def test_basic_velocity(test_args):
    result = angular_velocity_def.calculate_angular_velocity(test_args.a0, test_args.a1,
        test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.time)
    result_velocity = convert_to(result, angular_velocity_def.definition_units_SI).evalf(2)
    assert_approx(result_velocity, (pi / 5).evalf(4))


def test_velocity_with_number(test_args):
    angular_velocity_def.calculate_angular_velocity(100, test_args.a1, test_args.t)
    angular_velocity_def.calculate_angular_velocity(test_args.a0, 100, test_args.t)


def test_velocity_with_bad_angle(test_args):
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_def.calculate_angular_velocity(ab, test_args.a1, test_args.t)
    with raises(errors.UnitsError):
        angular_velocity_def.calculate_angular_velocity(test_args.a0, ab, test_args.t)


def test_velocity_with_bad_time(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_def.calculate_angular_velocity(test_args.a0, test_args.a1, tb)
    with raises(TypeError):
        angular_velocity_def.calculate_angular_velocity(test_args.a0, test_args.a1, 100)
