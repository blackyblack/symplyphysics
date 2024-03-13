from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal)

from symplyphysics.laws.astronomy import angle_of_rotation_during_e_of_rotation_during_gravitational_maneuver as angle_law

## The mass of the Earth is 5.97e24 kilograms, the aiming range is equal to the geodetic orbit of the Earth 43,000 kilometers.
## Then, at a rocket speed equal to the second cosmic velocity of the Earth 11.2 kilometers per second, the angle of rotation will be 0.1475 radians.

Args = namedtuple("Args", ["planet_mass", "aiming_range", "rocket_speed"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    planet_mass = Quantity(5.97e24 * units.kilograms)
    aiming_range = Quantity(43000 * units.kilometers)
    rocket_speed = Quantity(11.2 * units.kilometers / units.second)
    return Args(planet_mass=planet_mass,
        aiming_range=aiming_range,
        rocket_speed=rocket_speed)


def test_basic_angle(test_args: Args) -> None:
    result = angle_law.calculate_angle(test_args.planet_mass, test_args.aiming_range,
        test_args.rocket_speed)
    assert_equal(result, 0.1475 * units.radians)


def test_bad_planet_mass(test_args: Args) -> None:
    bad_planet_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angle_law.calculate_angle(bad_planet_mass, test_args.aiming_range,
            test_args.rocket_speed)
    with raises(TypeError):
        angle_law.calculate_angle(100, test_args.aiming_range,
            test_args.rocket_speed)


def test_bad_aiming_range(test_args: Args) -> None:
    bad_aiming_range = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angle_law.calculate_angle(test_args.planet_mass, bad_aiming_range,
            test_args.rocket_speed)
    with raises(TypeError):
        angle_law.calculate_angle(test_args.planet_mass, 100,
            test_args.rocket_speed)


def test_bad_rocket_speed(test_args: Args) -> None:
    bad_rocket_speed = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        angle_law.calculate_angle(test_args.planet_mass, test_args.aiming_range,
            bad_rocket_speed)
    with raises(TypeError):
        angle_law.calculate_angle(test_args.planet_mass, test_args.aiming_range,
            100)

