from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity import magnetic_field_intensity_of_infinite_wire as intensity_law

# Description
## With a current value of 0.5 ampere and a distance of 2 meters from the wire, the magnetic
## field intensity will be equal to 39.79 milliamperes per meter.
## https://www.fxyz.ru/формулы_по_физике/электричество/магнитное_поле/правило_буравчика_электромагнетизм/напряженность_магнитного_поля/

Args = namedtuple("Args", ["current", "distance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    current = Quantity(0.5 * units.ampere)
    distance = Quantity(2 * units.meter)
    return Args(current=current, distance=distance)


def test_basic_magnetic_intensity(test_args: Args) -> None:
    result = intensity_law.calculate_magnetic_intensity(test_args.current, test_args.distance)
    assert_equal(result, 39.79 * prefixes.milli * units.ampere / units.meter)


def test_bad_current(test_args: Args) -> None:
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_magnetic_intensity(current, test_args.distance)
    with raises(TypeError):
        intensity_law.calculate_magnetic_intensity(100, test_args.distance)


def test_bad_distance(test_args: Args) -> None:
    distance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_magnetic_intensity(test_args.current, distance)
    with raises(TypeError):
        intensity_law.calculate_magnetic_intensity(test_args.current, 100)
