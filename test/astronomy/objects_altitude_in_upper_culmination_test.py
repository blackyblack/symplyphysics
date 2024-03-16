from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors,)
from symplyphysics.laws.astronomy import objects_altitude_in_upper_culmination as altitude_law

# Description
## The declination of the object is -26 degrees. The latitude is 50 degrees. Then the height of the luminary
## at the upper culmination is 14 degrees.
## https://infourok.ru/kak-reshat-zadachi-po-astronomii-na-primere-resheniya-zadachi-opredelenie-visoti-svetila-v-verhney-kulminacii-3977106.html

Args = namedtuple("Args", ["latitude", "declination"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    latitude = Quantity(50 * units.deg)
    declination = Quantity(-26 * units.deg)

    return Args(latitude=latitude, declination=declination)


def test_basic_altitude(test_args: Args) -> None:
    result = altitude_law.calculate_altitude(test_args.latitude, test_args.declination)
    assert_equal(result, 14 * units.deg)


def test_bad_latitude(test_args: Args) -> None:
    latitude = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        altitude_law.calculate_altitude(latitude, test_args.declination)


def test_bad_declination(test_args: Args) -> None:
    declination = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        altitude_law.calculate_altitude(test_args.latitude, declination)
