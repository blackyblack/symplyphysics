from collections import namedtuple
from sympy import re, im
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, convert_to, Quantity, SI)
from symplyphysics.laws.electricity import coil_impedance_is_inductivity_and_frequency_proportional as coil_impedance_law

# Description
## Assert we have a coil with 5H inductivity. In 50Hz (314.159 rad/s) circuit impedance of this coil should be 1575.8 Ohm.
## (https://www.translatorscafe.com/unit-converter/ru-RU/calculator/inductor-impedance/)

@fixture(name="test_args")
def test_args_fixture():
    inductivity = Quantity(5 * units.henry)
    frequency = Quantity(314.159 * units.radian / units.second)
    Args = namedtuple("Args", ["inductivity", "frequency"])
    return Args(inductivity=inductivity, frequency=frequency)


def test_basic_impedance(test_args):
    result = coil_impedance_law.calculate_impedance(test_args.inductivity, test_args.frequency)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_impedance = convert_to(result, units.ohm).evalf(5)
    result_re = re(result_impedance)
    result_im = im(result_impedance)
    assert result_re == 0
    assert result_im == approx(1575.8, 0.1)

def test_bad_inductivity(test_args):
    ib = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coil_impedance_law.calculate_impedance(ib, test_args.frequency)
    with raises(TypeError):
        coil_impedance_law.calculate_impedance(100, test_args.frequency)


def test_bad_frequency(test_args):
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coil_impedance_law.calculate_impedance(test_args.inductivity, wb)
    with raises(TypeError):
        coil_impedance_law.calculate_impedance(test_args.inductivity, 100)
