from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
    angle_type,
)
from symplyphysics.laws.kinematic import accelerated_angular_velocity_from_time as angular_velocity_law

# Description
## A body is rotating about a fixed axis with angular acceleration of 3 rad/s^2. Initially its
## angular velocity is -3 rad/s. In 3 seconds, its angular velocity amounts to 6 rad/s^2.


@fixture(name="test_args")
def test_args_fixture():
    w0 = Quantity(-3.0 * units.radian / units.second)
    alpha = Quantity(3.0 * units.radian / units.second**2)
    t = Quantity(3.0 * units.second)
    Args = namedtuple("Args", "w0 alpha t")
    return Args(w0=w0, alpha=alpha, t=t)


def test_basic_law(test_args):
    result = angular_velocity_law.calculate_angular_velocity(test_args.w0, test_args.alpha, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.time)
    result_value = convert_to(result, units.radian / units.second).evalf(3)
    assert result_value == approx(6.0, 1e-3)


def test_bad_angular_velocity(test_args):
    w0b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_law.calculate_angular_velocity(w0b, test_args.alpha, test_args.t)
    with raises(TypeError):
        angular_velocity_law.calculate_angular_velocity(100, test_args.alpha, test_args.t)


def test_bad_angular_acceleration(test_args):
    alpha_b = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_law.calculate_angular_velocity(test_args.w0, alpha_b, test_args.t)
    with raises(TypeError):
        angular_velocity_law.calculate_angular_velocity(test_args.w0, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        angular_velocity_law.calculate_angular_velocity(test_args.w0, test_args.alpha, tb)
    with raises(TypeError):
        angular_velocity_law.calculate_angular_velocity(test_args.w0, test_args.alpha, 100)
