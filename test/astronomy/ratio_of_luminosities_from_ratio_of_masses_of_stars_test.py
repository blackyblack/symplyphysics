from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import ratio_of_luminosities_from_ratio_of_masses_of_stars as illuminance_law

# Description
## The mass of the first star is 1.98847e30 kilograms, and the luminosity is 2 [joule / meter^2]. The mass of the second star is 3.97694e30 kilograms.
## Then the luminosity of the second star will be 32 [joule / meter^2].
## https://college.ru/astronomy/course/content/chapter6/section1/paragraph5/theory.html

Args = namedtuple("Args", ["mass_first", "illuminance_first", "mass_second"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass_first = Quantity(1.98847e30 * units.kilogram)
    mass_second = Quantity(3.97694e30 * units.kilogram)
    illuminance_first = Quantity(2 * units.joule / units.second)

    return Args(
        mass_first=mass_first,
        mass_second=mass_second,
        illuminance_first=illuminance_first,
    )


def test_basic_illuminance_second(test_args: Args) -> None:
    result = illuminance_law.calculate_illuminance_second(test_args.mass_first,
        test_args.mass_second, test_args.illuminance_first)
    assert_equal(result, 32 * units.joule / units.second)


def test_bad_mass(test_args: Args) -> None:
    bad_mass = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        illuminance_law.calculate_illuminance_second(bad_mass, test_args.mass_second,
            test_args.illuminance_first)
    with raises(TypeError):
        illuminance_law.calculate_illuminance_second(100, test_args.mass_second,
            test_args.illuminance_first)
    with raises(errors.UnitsError):
        illuminance_law.calculate_illuminance_second(test_args.mass_first, bad_mass,
            test_args.illuminance_first)
    with raises(TypeError):
        illuminance_law.calculate_illuminance_second(test_args.mass_first, 100,
            test_args.illuminance_first)


def test_bad_illuminance_first(test_args: Args) -> None:
    illuminance_first = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        illuminance_law.calculate_illuminance_second(test_args.mass_first, test_args.mass_second,
            illuminance_first)
    with raises(TypeError):
        illuminance_law.calculate_illuminance_second(test_args.mass_first, test_args.mass_second,
            100)
