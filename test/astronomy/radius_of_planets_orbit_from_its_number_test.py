from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import radius_of_planets_orbit_from_its_number as radius_law

# Description
## The first number of the planet corresponds to the Earth. Then the radius of the orbit will be equal to 1 astronomical unit.

Args = namedtuple("Args", ["number_of_planet"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    number_of_planet = 1

    return Args(number_of_planet=number_of_planet)


def test_basic_radius_of_orbit(test_args: Args) -> None:
    result = radius_law.calculate_radius_of_orbit(test_args.number_of_planet)
    assert_equal(result, 1 * units.astronomical_unit)

def test_bad_number_of_planet() -> None:
    number_of_planet = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius_of_orbit(number_of_planet)
