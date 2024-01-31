from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.waves import wavelength_from_temperature as wins_law

# Description
## Assert we are having preheated black object with temperature of 4000K.
## With online calculator
## https://www.calculatoratoz.com/ru/wavelength-corresponding-to-maximum-radiation-emission-calculator/Calc-40123
## we obtain the most intensive radiation at 7.2444e-7 meters wavelength.


@fixture(name="test_args")
def test_args_fixture():
    T = Quantity(4000 * units.kelvin)
    Args = namedtuple("Args", ["T"])
    return Args(T=T)


def test_basic_wavelength(test_args):
    result = wins_law.calculate_intensive_wavelength(test_args.T)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_current = convert_to(result, units.meter).evalf(6)
    assert result_current == approx(7.2444e-7, 0.00001)


def test_bad_temperature():
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wins_law.calculate_intensive_wavelength(tb)
    with raises(TypeError):
        wins_law.calculate_intensive_wavelength(100)

