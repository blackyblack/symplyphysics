from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import latitude_from_zenith_distances_and_declinations as latitude_law

# Description
## If the zenith distance and declination of the northern luminary are equal to 30 and 25 degrees, respectively,
## and the zenith distance and declination of the southern luminary are equal to 150 and -15 degrees, respectively,
## then the latitude is 65 degrees.

Args = namedtuple("Args", ["zenith_distance_north", "zenith_distance_south", "declination_north", "declination_south",])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    zenith_distance_north = 30
    zenith_distance_south = 150
    declination_north = 25
    declination_south = -15

    return Args(zenith_distance_north=zenith_distance_north,
        zenith_distance_south=zenith_distance_south,
        declination_north=declination_north,
        declination_south=declination_south,)


def test_basic_latitude(test_args: Args) -> None:
    result = latitude_law.calculate_latitude(test_args.zenith_distance_north,
        test_args.zenith_distance_south, test_args.declination_north, test_args.declination_south)
    assert_equal(result, 65)


def test_bad_zenith_distance(test_args: Args) -> None:
    bad_zenith_distance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        latitude_law.calculate_latitude(bad_zenith_distance, test_args.zenith_distance_south,
            test_args.declination_north, test_args.declination_south)
    with raises(errors.UnitsError):
        latitude_law.calculate_latitude(test_args.zenith_distance_north, bad_zenith_distance,
            test_args.declination_north, test_args.declination_south)


def test_bad_declination(test_args: Args) -> None:
    bad_declination = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        latitude_law.calculate_latitude(test_args.zenith_distance_north, test_args.zenith_distance_south,
            bad_declination, test_args.declination_south)
    with raises(errors.UnitsError):
        latitude_law.calculate_latitude(test_args.zenith_distance_north,
            test_args.zenith_distance_south, test_args.declination_north, bad_declination)
