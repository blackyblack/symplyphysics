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
from symplyphysics.laws.kinematic import linear_velocity_from_angular_velocity_and_radius as linear_velocity_law

# Description
## A rotating cake stand has a circular birthday cake with a candle on it. When the stand is rotating, the candle has an angular velocity of
## 0.6 rad/s. If the candle is located 5 cm from the center of the cake, what is the linear velocity of the candle?


@fixture(name="test_args")
def test_args_fixture():
    V1 = Quantity(0.6 * units.radian / units.second)
    R1 = Quantity(5 * units.centimeter)
    Args = namedtuple("Args", ["V1", "R1"])
    return Args(V1=V1, R1=R1)


def test_basic_velocity(test_args):
    result = linear_velocity_law.calculate_linear_velocity(test_args.V1, test_args.R1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    result_velocity = convert_to(result, units.centimeter / units.second).evalf(2)
    assert_approx(result_velocity, 3)


def test_bad_angular_velocity(test_args):
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        linear_velocity_law.calculate_linear_velocity(Vb, test_args.R1)
    with raises(TypeError):
        linear_velocity_law.calculate_linear_velocity(100, test_args.R1)


def test_bad_radius(test_args):
    Rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        linear_velocity_law.calculate_linear_velocity(test_args.V1, Rb)
    with raises(TypeError):
        linear_velocity_law.calculate_linear_velocity(test_args.V1, 100)
