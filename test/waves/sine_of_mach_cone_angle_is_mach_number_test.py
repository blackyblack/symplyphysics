from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import sine_of_mach_cone_angle_is_mach_number as mach_angle_law

# Description
## An object is moving with ultrasonic speed and is characterized by Mach number M = 2.5.
## The Mach cone angle is 0.411 rad.

Args = namedtuple("Args", "m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = 2.5
    return Args(m=m)


def test_law(test_args: Args) -> None:
    result = mach_angle_law.calculate_mach_cone_angle(test_args.m)
    assert_equal(result, 0.411 * units.radian, tolerance=2e-3)


def test_bad_mach_number() -> None:
    mb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        mach_angle_law.calculate_mach_cone_angle(mb)

    mb = 0.5
    with raises(ValueError):
        mach_angle_law.calculate_mach_cone_angle(mb)
