from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors,)
from symplyphysics.laws.astronomy import maximum_angle_of_rotation_during_gravitational_maneuver as angle_law

# Description
## The Earth's first cosmic velocity is 7.91 kilometers per second. Let the speed of the rocket be equal to the Earth's
## second cosmic velocity of 11.2 kilometers per second. Then the maximum angle of rotation during a gravitational maneuver
## is 0.463 radians.

Args = namedtuple("Args", ["first_cosmic_velocity_planet", "rocket_speed"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    first_cosmic_velocity_planet = Quantity(7.91 * units.kilometer / units.second)
    rocket_speed = Quantity(11.2 * units.kilometer / units.second)

    return Args(first_cosmic_velocity_planet=first_cosmic_velocity_planet, rocket_speed=rocket_speed)


def test_basic_maximum_angle(test_args: Args) -> None:
    result = angle_law.calculate_maximum_angle(test_args.first_cosmic_velocity_planet, test_args.rocket_speed)
    assert_equal(result, 0.463 * units.radian)


def test_bad_first_cosmic_velocity_planet(test_args: Args) -> None:
    first_cosmic_velocity_planet = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angle_law.calculate_maximum_angle(first_cosmic_velocity_planet, test_args.rocket_speed)
    with raises(TypeError):
        angle_law.calculate_maximum_angle(100, test_args.rocket_speed)


def test_bad_rocket_speed(test_args: Args) -> None:
    rocket_speed = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angle_law.calculate_maximum_angle(test_args.first_cosmic_velocity_planet, rocket_speed)
    with raises(TypeError):
        angle_law.calculate_maximum_angle(test_args.first_cosmic_velocity_planet, 100)
