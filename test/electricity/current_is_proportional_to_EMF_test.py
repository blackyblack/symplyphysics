from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import current_is_proportional_to_EMF as ohms_law_full_circuit

# Description
## Assert we have 4 Volts applied to 3-Ohm resistor and 2-Ohm resistor in power supply.
## According to Ohm's Law for full circuit we should have 4/(3+2) = 0.8Amps current flowing through this circuit.


@fixture(name="test_args")
def test_args_fixture():
    V = Quantity(4 * units.volt)
    R = Quantity(3 * units.ohm)
    r = Quantity(2 * units.ohm)
    Args = namedtuple("Args", ["V", "R", "r"])
    return Args(V=V, R=R, r=r)


def test_basic_current(test_args):
    result = ohms_law_full_circuit.calculate_current(test_args.V, test_args.R, test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current)
    result_current = convert_to(result, units.ampere).evalf(2)
    assert result_current == approx(0.8, 0.01)

def test_bad_voltage(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        ohms_law_full_circuit.calculate_current(Vb, test_args.R, test_args.r)
    with raises(TypeError):
        ohms_law_full_circuit.calculate_current(100, test_args.R, test_args.r)

def test_bad_resistance_R(test_args):
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        ohms_law_full_circuit.calculate_current(test_args.V, Rb, test_args.r)
    with raises(TypeError):
        ohms_law_full_circuit.calculate_current(test_args.V, 100, test_args.r)

def test_bad_resistance_r(test_args):
    rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        ohms_law_full_circuit.calculate_current(test_args.V, test_args.R, rb)
    with raises(TypeError):
        ohms_law_full_circuit.calculate_current(test_args.V, test_args.R, 100)