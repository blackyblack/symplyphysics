from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import absolute_magnitude_of_apparent_magnitude_and_distance as magnitude_law

# Description
## With an apparent magnitude of 29 and a distance of 1.443854e11 astronomical units, the absolute magnitude will be 4.775.
## http://www.astro.spbu.ru/staff/viva/Book/ch4L/node8.html

Args = namedtuple("Args", ["apparent_magnitude", "distance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    apparent_magnitude = 29
    distance = Quantity(1.443854e11 * units.astronomical_unit)

    return Args(apparent_magnitude=apparent_magnitude, distance=distance)


def test_basic_absolute_magnitude(test_args: Args) -> None:
    result = magnitude_law.calculate_absolute_magnitude(test_args.apparent_magnitude,
        test_args.distance)
    assert_equal(result, 4.775)


def test_bad_apparent_magnitude(test_args: Args) -> None:
    apparent_magnitude = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        magnitude_law.calculate_absolute_magnitude(apparent_magnitude, test_args.distance)


def test_bad_distance(test_args: Args) -> None:
    distance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        magnitude_law.calculate_absolute_magnitude(test_args.apparent_magnitude, distance)
    with raises(TypeError):
        magnitude_law.calculate_absolute_magnitude(test_args.apparent_magnitude, 100)
