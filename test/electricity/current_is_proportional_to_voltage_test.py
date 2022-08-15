# Description
## Assert we have 3 Volts applied to 2-Ohm resistor.
## According to Ohm's Law we should have 3/2 = 1.5Amps current flowing through this resistor.

from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.electricity import current_through_element_is_proporional_to_voltage_und_reversly_proportional_to_resistance as ohms_law

@fixture
def test_args():
    Voltage = units.Quantity('Voltage')
    Resistance = units.Quantity('Resistance')
    SI.set_quantity_dimension(Voltage, units.voltage)
    SI.set_quantity_scale_factor(Voltage, 3 * units.volt)
    SI.set_quantity_dimension(Resistance, units.impedance)
    SI.set_quantity_scale_factor(Resistance, 2 * units.ohm)

    Args = namedtuple('Args', ['Voltage', 'Resistance'])
    return Args(Voltage = Voltage, Resistance = Resistance)

def test_basic_current(test_args):
    result = ohms_law.calculate_current(test_args.Voltage, test_args.Resistance)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current)
    result_current = convert_to(result, units.ampere).subs(units.ampere, 1).evalf(2)
    assert result_current == approx(1.5, 0.01)

def test_bad_voltage(test_args):
    Vb = units.Quantity('Vb')
    SI.set_quantity_dimension(Vb, units.length)
    SI.set_quantity_scale_factor(Vb, 1 * units.meter)

    with raises(errors.UnitsError):
        ohms_law.calculate_current(Vb, test_args.Resistance)

    with raises(TypeError):
        ohms_law.calculate_current(100, test_args.Resistance)

def test_bad_resistance(test_args):
    Rb = units.Quantity('Rb')
    SI.set_quantity_dimension(Rb, units.length)
    SI.set_quantity_scale_factor(Rb, 1 * units.meter)

    with raises(errors.UnitsError):
        ohms_law.calculate_current(test_args.Voltage, Rb)

    with raises(TypeError):
        ohms_law.calculate_current(test_args.Voltage, 100)