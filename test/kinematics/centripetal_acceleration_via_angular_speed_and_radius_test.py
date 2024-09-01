from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematics import (
    centripetal_acceleration_via_angular_speed_and_radius as centripetal_law)

# Description
## A body is rotating about an axis with an angular velocity of 5 rad/s. Its perpendicular
## distance to the axis is 2 m. Then its centripetal acceleration is 50 m/s**2.

Args = namedtuple("Args", "w r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w = Quantity(5.0 * units.radian / units.second)
    r = Quantity(2.0 * units.meter)
    return Args(w=w, r=r)


def test_law(test_args: Args) -> None:
    result = centripetal_law.calculate_centripetal_acceleration(test_args.w, test_args.r)
    assert_equal(result, 50.0 * units.meter / units.second**2)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        centripetal_law.calculate_centripetal_acceleration(wb, test_args.r)
    with raises(TypeError):
        centripetal_law.calculate_centripetal_acceleration(100, test_args.r)


def test_bad_curve_radius(test_args: Args) -> None:
    rb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        centripetal_law.calculate_centripetal_acceleration(test_args.w, rb)
    with raises(TypeError):
        centripetal_law.calculate_centripetal_acceleration(test_args.w, 100)
