from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import speed_of_galaxy_from_distance_to_galaxy as speed_law

# Description
## Let the distance to the galaxy be 4.6285e+22 kilometers. Then the speed of the galaxy is 101827 kilometer per second.
## http://www.astro.spbu.ru/staff/viva/Book/ch4L/node17.html

Args = namedtuple("Args", ["distance_to_galaxy"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    distance_to_galaxy = Quantity(4.6285e+22 * units.kilometer)

    return Args(distance_to_galaxy=distance_to_galaxy)


def test_basic_speed(test_args: Args) -> None:
    result = speed_law.calculate_speed(test_args.distance_to_galaxy)
    assert_equal(result, 101827 * units.kilometer / units.second)


def test_bad_distance_to_galaxy() -> None:
    distance_to_galaxy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_speed(distance_to_galaxy)
    with raises(TypeError):
        speed_law.calculate_speed(100)
