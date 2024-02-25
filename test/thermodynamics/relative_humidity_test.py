from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.thermodynamics import relative_humidity as humidity_law

# Description
## With a water vapor pressure of 1340 pascal and a saturated vapor pressure of 2984 pascal,
## the relative humidity will be 44.9 %.
## http://profil.adu.by/mod/book/view.php?id=3196&chapterid=9248

Args = namedtuple("Args", ["water_vapor_pressure", "saturated_vapor_pressure"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    water_vapor_pressure = Quantity(1340 * units.pascal)
    saturated_vapor_pressure = Quantity(2984 * units.pascal)

    return Args(water_vapor_pressure=water_vapor_pressure,
        saturated_vapor_pressure=saturated_vapor_pressure)


def test_basic_relative_humidity(test_args: Args) -> None:
    result = humidity_law.calculate_relative_humidity(test_args.water_vapor_pressure, test_args.saturated_vapor_pressure)
    assert_equal(result, 0.449)


def test_bad_pressure(test_args: Args) -> None:
    bad_pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        humidity_law.calculate_relative_humidity(bad_pressure, test_args.saturated_vapor_pressure)
    with raises(TypeError):
        humidity_law.calculate_relative_humidity(100, test_args.saturated_vapor_pressure)
    with raises(errors.UnitsError):
        humidity_law.calculate_relative_humidity(test_args.water_vapor_pressure, bad_pressure)
    with raises(TypeError):
        humidity_law.calculate_relative_humidity(test_args.water_vapor_pressure, 100)
