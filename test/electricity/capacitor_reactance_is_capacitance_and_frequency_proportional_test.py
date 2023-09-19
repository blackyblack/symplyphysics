from collections import namedtuple
from sympy import re, im
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, convert_to, Quantity, SI)
from symplyphysics.laws.electricity import capacitor_reactance_is_capacitance_and_frequency_proportional as capacitor_reactance_law

# Description
## Assert we have a capacitor with 5F capacitance. In 0.05Hz (0.314159 rad/s) circuit reactance of this element should be 636.61mOhm.
## (https://www.translatorscafe.com/unit-converter/ru-RU/calculator/capacitor-impedance/)

@fixture(name="test_args")
def test_args_fixture():
    capacitance = Quantity(5 * units.farad)
    frequency = Quantity(0.314159 * units.radian / units.second)
    Args = namedtuple("Args", ["capacitance", "frequency"])
    return Args(capacitance=capacitance, frequency=frequency)


def test_basic_impedance(test_args):
    result = capacitor_reactance_law.calculate_reactance(test_args.capacitance, test_args.frequency)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_impedance = convert_to(result, units.ohm).evalf(5)
    result_re = re(result_impedance)
    result_im = abs(im(result_impedance))
    assert result_re == 0
    assert result_im == approx(0.63662, 0.001)

def test_bad_capacitance(test_args):
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitor_reactance_law.calculate_reactance(cb, test_args.frequency)
    with raises(TypeError):
        capacitor_reactance_law.calculate_reactance(100, test_args.frequency)


def test_bad_frequency(test_args):
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitor_reactance_law.calculate_reactance(test_args.capacitance, wb)
    with raises(TypeError):
        capacitor_reactance_law.calculate_reactance(test_args.capacitance, 100)
