from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import tangential_acceleration_of_rotating_body as tangential_acceleration_law

# Description
## A body is moving along a curve. Its angular acceleration is 3.5 rad/s^2 and its rotation radius
## around the momentary rotation axis is 0.1 m. The magnitude of its tangential acceleration should
## amount to 0.35 m/s^2.


@fixture(name="test_args")
def test_args_fixture():
    alpha = Quantity(3.5 * units.radian / units.second**2)
    r = Quantity(0.1 * units.meter)
    Args = namedtuple("Args", "alpha r")
    return Args(alpha=alpha, r=r)


def test_basic_law(test_args):
    result = tangential_acceleration_law.calculate_tangential_acceleration(
        test_args.alpha, test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.acceleration)
    result_value = convert_to(result, units.meter / units.second**2).evalf(3)
    assert_approx(result_value, 0.35)


def test_bad_angular_acceleration(test_args):
    alpha_b = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        tangential_acceleration_law.calculate_tangential_acceleration(alpha_b, test_args.r)
    with raises(TypeError):
        tangential_acceleration_law.calculate_tangential_acceleration(100, test_args.r)


def test_bad_rotation_radius(test_args):
    rb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        tangential_acceleration_law.calculate_tangential_acceleration(test_args.alpha, rb)
    with raises(TypeError):
        tangential_acceleration_law.calculate_tangential_acceleration(test_args.alpha, 100)
