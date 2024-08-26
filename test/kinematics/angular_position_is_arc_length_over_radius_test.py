from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.laws.kinematics import angular_position_is_arc_length_over_radius as angular_position_def

# Description
## A rigid body is rotating about an axis. At the radius of 1 m, the arc length of the its path is measured
## to be 2 m. Thus, the angular position of the body relative to its initial position is 2 rad.

Args = namedtuple("Args", "s r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s = Quantity(2 * units.meter)
    r = Quantity(1 * units.meter)
    return Args(s=s, r=r)


def test_basic_law(test_args: Args) -> None:
    result = angular_position_def.calculate_angular_position(test_args.s, test_args.r)
    assert_equal(result, 2 * units.radian)


def test_bad_arc_length(test_args: Args) -> None:
    sb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        angular_position_def.calculate_angular_position(sb, test_args.r)
    with raises(TypeError):
        angular_position_def.calculate_angular_position(100, test_args.r)


def test_bad_radius(test_args: Args) -> None:
    rb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        angular_position_def.calculate_angular_position(test_args.s, rb)
    with raises(TypeError):
        angular_position_def.calculate_angular_position(test_args.s, 100)
