from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, Quantity, errors
from symplyphysics.laws.gravity import radius_of_geostationary_orbit as radius_law

# Description
## The mass of the Earth is 5.9742e24 kilograms. At a rotational speed of 7.29e-5 radians per second,
## the radius of the orbit is 42164 kilometers.
## https://ru.wikipedia.org/wiki/Геостационарная_орбита

Args = namedtuple("Args", ["mass_of_planet", "speed_rotation_satellite"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass_of_planet = Quantity(5.9742e24 * units.kilogram)
    speed_rotation_satellite = Quantity(7.29e-5 * units.radian / units.second)

    return Args(mass_of_planet=mass_of_planet, speed_rotation_satellite=speed_rotation_satellite)


def test_basic_radius_of_orbit(test_args: Args) -> None:
    result = radius_law.calculate_radius_of_orbit(test_args.mass_of_planet,
        test_args.speed_rotation_satellite)
    assert_equal(result, 42164 * units.kilometer)


def test_bad_mass_of_planet(test_args: Args) -> None:
    mass_of_planet = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius_of_orbit(mass_of_planet, test_args.speed_rotation_satellite)
    with raises(TypeError):
        radius_law.calculate_radius_of_orbit(100, test_args.speed_rotation_satellite)


def test_bad_speed_rotation_satellite(test_args: Args) -> None:
    speed_rotation_satellite = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius_of_orbit(test_args.mass_of_planet, speed_rotation_satellite)
    with raises(TypeError):
        radius_law.calculate_radius_of_orbit(test_args.mass_of_planet, 100)
