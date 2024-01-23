from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.condensed_matter import resistance_from_temperature as resistance_law

# Description
## With a temperature coefficient of 15 [1 / kelvin], a temperature of 573.15 [kelvin]
## and an initial resistance of [25 ohm], the resistance will be 112525 [ohm].
## https://www.calculatoratoz.com/ru/temperature-dependence-of-resistance-calculator/Calc-2232


@fixture(name="test_args")
def test_args_fixture():
    resistance_initial = Quantity(25 * units.ohm)
    temperature_coefficient = Quantity(15 * (1 / units.kelvin))
    temperature = Quantity(573.15 * units.kelvin)

    Args = namedtuple("Args",
        ["resistance_initial", "temperature_coefficient", "temperature"])
    return Args(resistance_initial=resistance_initial,
        temperature_coefficient=temperature_coefficient,
        temperature=temperature)


def test_basic_resistance(test_args):
    result = resistance_law.calculate_resistance(test_args.resistance_initial,
        test_args.temperature_coefficient, test_args.temperature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result = convert_to(result, units.ohm).evalf(5)
    assert result == approx(112525, rel=0.01)


def test_bad_resistance_initial(test_args):
    resistance_initial = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(resistance_initial, test_args.temperature_coefficient,
            test_args.temperature)
    with raises(TypeError):
        resistance_law.calculate_resistance(100, test_args.temperature_coefficient, test_args.temperature)


def test_bad_temperature_coefficient(test_args):
    temperature_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.resistance_initial, temperature_coefficient,
            test_args.temperature)
    with raises(TypeError):
        resistance_law.calculate_resistance(test_args.resistance_initial, 100,
            test_args.temperature)


def test_bad_temperature(test_args):
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.resistance_initial,
            test_args.temperature_coefficient, temperature)
    with raises(TypeError):
        resistance_law.calculate_resistance(test_args.resistance_initial,
            test_args.temperature_coefficient, 100)
