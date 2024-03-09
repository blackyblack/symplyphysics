from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors,)
from symplyphysics.laws.astronomy import large_semi_axis_of_orbit as axis_law

# Description
## The mass of the sun is 1.989e30 kilograms, and the orbital speed of the Earth is 29.8 kilometers per second.
## Then the semi-major axis of the orbit is equal to 149.5e9 meters.
## https://uchet-jkh.ru/i/kak-naiti-bolsuyu-poluos-orbity/

Args = namedtuple("Args", ["orbital_velocity", "mass"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    orbital_velocity = Quantity(29.8 * (units.kilometer / units.second))
    mass = Quantity(1.989e30 * units.kilogram)

    return Args(orbital_velocity=orbital_velocity, mass=mass)


def test_basic_large_semi_axis(test_args: Args) -> None:
    result = axis_law.calculate_large_semi_axis(test_args.orbital_velocity, test_args.mass)
    assert_equal(result, 149.5e9 * units.meter)


def test_bad_orbital_velocity(test_args: Args) -> None:
    orbital_velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        axis_law.calculate_large_semi_axis(orbital_velocity, test_args.mass)
    with raises(TypeError):
        axis_law.calculate_large_semi_axis(100, test_args.mass)


def test_bad_mass(test_args: Args) -> None:
    mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        axis_law.calculate_large_semi_axis(test_args.orbital_velocity, mass)
    with raises(TypeError):
        axis_law.calculate_large_semi_axis(test_args.orbital_velocity, 100)
