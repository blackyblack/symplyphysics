from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors,)
from symplyphysics.laws.astronomy import approximate_time_of_stars_location_on_the_main_sequence as time_law

# Description
## The mass of the star is 3.978e30 kilograms, and its luminosity is 11.481e26 watts. Then the time of its stay on
## the main sequence will be equal to 6.66e9 years.

Args = namedtuple("Args", ["mass_of_star", "luminosity_of_star"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass_of_star = Quantity(3.978e30 * units.kilogram)
    luminosity_of_star = Quantity(11.481e26 * units.watt)

    return Args(mass_of_star=mass_of_star, luminosity_of_star=luminosity_of_star)


def test_basic_lifetime(test_args: Args) -> None:
    result = time_law.calculate_lifetime(test_args.mass_of_star, test_args.luminosity_of_star)
    assert_equal(result, 6.66e9 * units.common_year)


def test_bad_mass_of_star(test_args: Args) -> None:
    mass_of_star = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_lifetime(mass_of_star, test_args.luminosity_of_star)
    with raises(TypeError):
        time_law.calculate_lifetime(100, test_args.luminosity_of_star)


def test_bad_luminosity_of_star(test_args: Args) -> None:
    luminosity_of_star = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        time_law.calculate_lifetime(test_args.mass_of_star, luminosity_of_star)
    with raises(TypeError):
        time_law.calculate_lifetime(test_args.mass_of_star, 100)
