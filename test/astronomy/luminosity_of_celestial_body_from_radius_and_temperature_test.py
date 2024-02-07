from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal
)

from symplyphysics.laws.astronomy import luminosity_of_celestial_body_from_radius_and_temperature as luminosity_of_celestial_body

## Source of numbers: https://astrogalaxy.ru/065.html

@fixture(name="test_args")
def test_args_fixture():
    radius = Quantity(695700000 * units.meter)
    temperature = Quantity(5772 * units.kelvin)
    Args = namedtuple("Args", ["radius", "temperature"])
    return Args(
        radius=radius,
        temperature=temperature
    )


def test_basic_law(test_args):
    result = luminosity_of_celestial_body.calculate_luminosity(test_args.radius, test_args.temperature)
    assert_equal(result, 3.827e26 * units.watt)


def test_bad_radius(test_args):
    bad_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        luminosity_of_celestial_body.calculate_luminosity(bad_radius, test_args.temperature)
    with raises(TypeError):
        luminosity_of_celestial_body.calculate_luminosity(100, test_args.temperature)


def test_bad_temperature(test_args):
    bad_temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        luminosity_of_celestial_body.calculate_luminosity(test_args.radius, bad_temperature)
    with raises(TypeError):
        luminosity_of_celestial_body.calculate_luminosity(test_args.radius, 100)
