from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, convert_to, Quantity, SI)
from symplyphysics.laws.electricity import coil_impedance_is_inductivity_and_frequency_proportional as coil_impedance_law

# Description
## Assert we have a coil with 5H inductivity. In 50Hz (314.159 s**-1) circuit impedance of this coil should be 1570.8 Ohm.
## (https://www.translatorscafe.com/unit-converter/ru-RU/calculator/inductor-impedance/)

@fixture(name="test_args")
def test_args_fixture():
    inductivity = Quantity(5 * units.henry)
    frequency = Quantity(315.159 * units.hertz)
    Args = namedtuple("Args", ["inductivity", "frequency"])
    return Args(inductivity=inductivity, frequency=frequency)


def test_basic_impedance(test_args):
    result = coil_impedance_law.calculate_impedance(test_args.inductivity, test_args.frequency)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_current = convert_to(result, units.ohm).evalf(6)
    assert result_current == approx(1570.8, 0.01)

'''
def test_bad_resistivity(test_args):
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coil_impedance_law.calculate_resistance(rb, test_args.wire_length, test_args.cross_section)
    with raises(TypeError):
        coil_impedance_law.calculate_resistance(100, test_args.wire_length, test_args.cross_section)


def test_bad_length(test_args):
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        wire_law.calculate_resistance(test_args.resistivity, lb, test_args.cross_section)
    with raises(TypeError):
        wire_law.calculate_resistance(test_args.resistivity, 100, test_args.cross_section)

'''
