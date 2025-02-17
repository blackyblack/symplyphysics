from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.astronomy import lifetime_of_star_on_main_sequence as lifetime_law

# Description
## The lifetime of the Sun on the main sequence is approximately 10 billion years. The indicator for the Sun is 4.75.
## The mass of the Sun in the masses of the Sun is 1.
## http://www.spacephys.ru/astronomicheskie-formuly#:~:text=n%20%3D%204%2C75%20%D0%B4%D0%BB%D1%8F%20%D0%B7%D0%B2%D1%91%D0%B7%D0%B4%20%D0%BC%D0%B0%D1%81%D1%81%D0%BE%D0%B9%200%2C7%20%E2%80%93%202%20%D0%9C%D1%81%D0%BE%D0%BB
## https://ru.wikipedia.org/wiki/%D0%A1%D0%BE%D0%BB%D0%BD%D1%86%D0%B5#:~:text=%D0%97%D0%B2%D0%B5%D0%B7%D0%B4%D0%B0%20%D1%82%D0%B0%D0%BA%D0%BE%D0%B9%20%D0%BC%D0%B0%D1%81%D1%81%D1%8B%2C%20%D0%BA%D0%B0%D0%BA%20%D0%A1%D0%BE%D0%BB%D0%BD%D1%86%D0%B5%2C%20%D0%B4%D0%BE%D0%BB%D0%B6%D0%BD%D0%B0%20%D1%81%D1%83%D1%89%D0%B5%D1%81%D1%82%D0%B2%D0%BE%D0%B2%D0%B0%D1%82%D1%8C%20%D0%BD%D0%B0%20%D0%B3%D0%BB%D0%B0%D0%B2%D0%BD%D0%BE%D0%B9%20%D0%BF%D0%BE%D1%81%D0%BB%D0%B5%D0%B4%D0%BE%D0%B2%D0%B0%D1%82%D0%B5%D0%BB%D1%8C%D0%BD%D0%BE%D1%81%D1%82%D0%B8%20%D0%B2%20%D0%BE%D0%B1%D1%89%D0%B5%D0%B9%20%D1%81%D0%BB%D0%BE%D0%B6%D0%BD%D0%BE%D1%81%D1%82%D0%B8%20%D0%BF%D1%80%D0%B8%D0%BC%D0%B5%D1%80%D0%BD%D0%BE%2010%C2%A0%D0%BC%D0%BB%D1%80%D0%B4%20%D0%BB%D0%B5%D1%82.

Args = namedtuple("Args", ["mass_of_star", "indicator"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass_of_star = Quantity(1.989e30 * units.kilogram)
    indicator = 4.75

    return Args(mass_of_star=mass_of_star, indicator=indicator)


def test_basic_lifetime(test_args: Args) -> None:
    result = lifetime_law.calculate_lifetime(test_args.mass_of_star, test_args.indicator)
    assert_equal(result, 10e9 * units.common_year, relative_tolerance=5e-3)


def test_bad_mass_of_star(test_args: Args) -> None:
    mass_of_star = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        lifetime_law.calculate_lifetime(mass_of_star, test_args.indicator)


def test_bad_indicator(test_args: Args) -> None:
    indicator = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        lifetime_law.calculate_lifetime(test_args.mass_of_star, indicator)
