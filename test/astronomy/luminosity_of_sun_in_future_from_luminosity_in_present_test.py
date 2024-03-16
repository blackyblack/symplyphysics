from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors,)
from symplyphysics.laws.astronomy import luminosity_of_sun_in_future_from_luminosity_in_present as luminosity_law

# Description
## Let the time be equal to 15 billion years. Then the luminosity of the Sun in the future will be equal to 2.863 units of the luminosity of the Sun.

Args = namedtuple("Args", ["luminosity_present", "time"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    luminosity_present = 1
    time = Quantity(15e9 * units.common_year)

    return Args(luminosity_present=luminosity_present, time=time)


def test_basic_luminosity_future(test_args: Args) -> None:
    result = luminosity_law.calculate_luminosity_future(test_args.luminosity_present, test_args.time)
    assert_equal(result, 2.863)


def test_bad_luminosity_present(test_args: Args) -> None:
    luminosity_present = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        luminosity_law.calculate_luminosity_future(luminosity_present, test_args.time)


def test_bad_time(test_args: Args) -> None:
    time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        luminosity_law.calculate_luminosity_future(test_args.luminosity_present, time)
    with raises(TypeError):
        luminosity_law.calculate_luminosity_future(test_args.luminosity_present, 100)
