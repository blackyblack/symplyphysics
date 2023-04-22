from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors, units, convert_to, Quantity, SI,
)
from symplyphysics.laws.electricity import current_is_proportional_to_voltage as ohms_law

# Description
## Assert we have 3 Volts applied to 2-Ohm resistor.
## According to Ohm's Law we should have 3/2 = 1.5Amps current flowing through this resistor.

@fixture
def test_args():
    Voltage = Quantity(3 * units.volt)
    Resistance = Quantity(2 * units.ohm)
    Args = namedtuple("Args", ["Voltage", "Resistance"])
    return Args(Voltage=Voltage, Resistance=Resistance)

def test_basic_current(test_args):
    result = ohms_law.calculate_current(test_args.Voltage, test_args.Resistance)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current)
    result_current = convert_to(result, units.ampere).subs(units.ampere, 1).evalf(2)
    assert result_current == approx(1.5, 0.01)

def test_bad_voltage(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        ohms_law.calculate_current(Vb, test_args.Resistance)
    with raises(TypeError):
        ohms_law.calculate_current(100, test_args.Resistance)

def test_bad_resistance(test_args):
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        ohms_law.calculate_current(test_args.Voltage, Rb)
    with raises(TypeError):
        ohms_law.calculate_current(test_args.Voltage, 100)
