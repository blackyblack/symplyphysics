from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.optics import interference_minimum as minimum_law

Args = namedtuple("Args", ["wave_length", "minimum_number"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    minimum_number = 2
    wave_length = Quantity(550 * prefixes.nano * units.meters)
    return Args(wave_length=wave_length, minimum_number=minimum_number)


def test_basic_travel_difference(test_args: Args) -> None:
    result = minimum_law.calculate_travel_difference(test_args.wave_length,
        test_args.minimum_number)
    assert_equal(result, 1375 * prefixes.nano * units.meters)


def test_bad_wave_length(test_args: Args) -> None:
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        minimum_law.calculate_travel_difference(lb, test_args.minimum_number)
    with raises(TypeError):
        minimum_law.calculate_travel_difference(100, test_args.minimum_number)


def test_bad_number_minimum(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        minimum_law.calculate_travel_difference(test_args.wave_length, nb)
