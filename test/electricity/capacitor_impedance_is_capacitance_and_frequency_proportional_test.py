from collections import namedtuple
from sympy import re, im
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, convert_to, Quantity, SI)
from symplyphysics.laws.electricity import capacitor_impedance_is_capacitance_and_frequency_proportional as capacitor_impedance_law

# Description
## Assert we have a capacitor with 5F capacitance. In 0.05Hz (0.314159 rad/s) circuit impedance of this element should be 636.61mOhm.
## (https://www.translatorscafe.com/unit-converter/ru-RU/calculator/capacitor-impedance/)

@fixture(name="test_args")
def test_args_fixture():
    c = Quantity(5 * units.farad)
    w = Quantity(0.314159 * units.radian / units.second)
    Args = namedtuple("Args", ["c", "w"])
    return Args(c=c, w=w)


def test_basic_impedance(test_args):
    result = capacitor_impedance_law.calculate_impedance(test_args.c, test_args.w)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_impedance = convert_to(result, units.ohm).evalf(5)
    result_re = re(result_impedance)
    result_im = im(result_impedance)
    assert result_re == 0
    assert result_im == approx(-0.63662, 0.001)

def test_bad_capacitance(test_args):
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitor_impedance_law.calculate_impedance(cb, test_args.w)
    with raises(TypeError):
        capacitor_impedance_law.calculate_impedance(100, test_args.w)


def test_bad_frequency(test_args):
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitor_impedance_law.calculate_impedance(test_args.c, wb)
    with raises(TypeError):
        capacitor_impedance_law.calculate_impedance(test_args.c, 100)
