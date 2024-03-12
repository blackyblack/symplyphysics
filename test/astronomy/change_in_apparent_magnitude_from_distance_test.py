from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import change_in_apparent_magnitude_from_distance as magnitude_law

# Description
## Let the illumination of the first object be 1 joule per square meter, the illumination of the second object is 3.512 joules per square meter.
## With the apparent magnitude of the first object equal to 3, the apparent magnitude of the second object is 1.636.
## http://www.astro.spbu.ru/staff/viva/Book/ch4L/node8.html

Args = namedtuple("Args", ["apparent_magnitude_first", "illuminance_first", "illuminance_second"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    apparent_magnitude_first = 3
    illuminance_first = Quantity(1 * units.joule / units.meter**2)
    illuminance_second = Quantity(3.512 * units.joule / units.meter**2)

    return Args(apparent_magnitude_first=apparent_magnitude_first,
        illuminance_first=illuminance_first,
        illuminance_second=illuminance_second)


def test_basic_apparent_magnitude_second(test_args: Args) -> None:
    result = magnitude_law.calculate_apparent_magnitude_second(test_args.apparent_magnitude_first,
        test_args.illuminance_first, test_args.illuminance_second)
    assert_equal(result, 1.636)


def test_bad_apparent_magnitude_first(test_args: Args) -> None:
    apparent_magnitude_first = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        magnitude_law.calculate_apparent_magnitude_second(apparent_magnitude_first, test_args.illuminance_first,
            test_args.illuminance_second)


def test_bad_illuminance(test_args: Args) -> None:
    bad_illuminance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        magnitude_law.calculate_apparent_magnitude_second(test_args.apparent_magnitude_first, bad_illuminance,
            test_args.illuminance_second)
    with raises(TypeError):
        magnitude_law.calculate_apparent_magnitude_second(test_args.apparent_magnitude_first, 100,
            test_args.illuminance_second)
    with raises(errors.UnitsError):
        magnitude_law.calculate_apparent_magnitude_second(test_args.apparent_magnitude_first,
            test_args.illuminance_first, bad_illuminance)
    with raises(TypeError):
        magnitude_law.calculate_apparent_magnitude_second(test_args.apparent_magnitude_first,
            test_args.illuminance_first, 100)
