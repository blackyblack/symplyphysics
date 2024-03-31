from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal)

from symplyphysics.laws.astronomy import orbital_velocity_from_length_of_large_half_axis_and_planet_mass as orbital_velocity

## Source of numbers: https://www.astronet.ru/db/msg/1175352/node11.html

Args = namedtuple("Args", ["planet_mass", "distance", "large_half_axis_length"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    planet_mass = Quantity(5.97e24 * units.kilograms)
    distance = Quantity(7302 * units.kilometers)
    large_half_axis_length = Quantity(13900 * units.kilometers)
    return Args(planet_mass=planet_mass,
        distance=distance,
        large_half_axis_length=large_half_axis_length)


def test_basic_law(test_args: Args) -> None:
    result = orbital_velocity.calculate_orbital_velocity(test_args.planet_mass, test_args.distance,
        test_args.large_half_axis_length)
    assert_equal(result, 8970 * units.meter / units.second)


def test_bad_planet_mass(test_args: Args) -> None:
    bad_planet_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        orbital_velocity.calculate_orbital_velocity(bad_planet_mass, test_args.distance,
            test_args.large_half_axis_length)
    with raises(TypeError):
        orbital_velocity.calculate_orbital_velocity(100, test_args.distance,
            test_args.large_half_axis_length)


def test_bad_length(test_args: Args) -> None:
    bad_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        orbital_velocity.calculate_orbital_velocity(test_args.planet_mass, bad_length,
            test_args.large_half_axis_length)
    with raises(errors.UnitsError):
        orbital_velocity.calculate_orbital_velocity(test_args.planet_mass, test_args.distance,
            bad_length)
    with raises(TypeError):
        orbital_velocity.calculate_orbital_velocity(test_args.planet_mass, 100,
            test_args.large_half_axis_length)
    with raises(TypeError):
        orbital_velocity.calculate_orbital_velocity(test_args.planet_mass, test_args.distance, 100)
    with raises(ValueError):
        orbital_velocity.calculate_orbital_velocity(test_args.planet_mass,
            test_args.large_half_axis_length, test_args.distance)
