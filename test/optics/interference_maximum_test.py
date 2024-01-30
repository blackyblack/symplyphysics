from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to, prefixes,
)
from symplyphysics.laws.optics import interference_maximum as maximum_law


@fixture(name="test_args")
def test_args_fixture():
    maximum_number = 2
    wave_length = Quantity(550 * prefixes.nano * units.meters)
    Args = namedtuple("Args", ["wave_length", "maximum_number"])
    return Args(wave_length=wave_length, maximum_number=maximum_number)


def test_basic_travel_difference(test_args):
    result = maximum_law.calculate_travel_difference(test_args.wave_length, test_args.maximum_number)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_travel_difference = convert_to(result, prefixes.nano * units.meters).evalf(5)
    assert result_travel_difference == approx(1100, 0.1)


def test_bad_wave_length(test_args):
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        maximum_law.calculate_travel_difference(lb, test_args.maximum_number)
    with raises(TypeError):
        maximum_law.calculate_travel_difference(100, test_args.maximum_number)


def test_bad_number_minimum(test_args):
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        maximum_law.calculate_travel_difference(test_args.wave_length, nb)
