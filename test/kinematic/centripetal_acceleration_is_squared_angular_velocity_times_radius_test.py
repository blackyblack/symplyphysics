from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic import (
    centripetal_acceleration_is_squared_angular_velocity_times_radius as centripetal_law
)

# Description
## A body is rotating about an axis with an angular velocity of 5 rad/s. Its perpendicular
## distance to the axis is 2 m. Then its centripetal acceleration is 50 m/s**2.


@fixture(name="test_args")
def test_args_fixture():
    w = Quantity(5.0 * units.radian / units.second)
    r = Quantity(2.0 * units.meter)
    Args = namedtuple("Args", "w r")
    return Args(w=w, r=r)


def test_law(test_args):
    result = centripetal_law.calculate_centripetal_acceleration(test_args.w, test_args.r)
    assert_equal(result, 50.0 * units.meter / units.second**2)


def test_bad_angular_velocity(test_args):
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        centripetal_law.calculate_centripetal_acceleration(wb, test_args.r)
    with raises(TypeError):
        centripetal_law.calculate_centripetal_acceleration(100, test_args.r)


def test_bad_curve_radius(test_args):
    rb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        centripetal_law.calculate_centripetal_acceleration(test_args.w, rb)
    with raises(TypeError):
        centripetal_law.calculate_centripetal_acceleration(test_args.w, 100)
