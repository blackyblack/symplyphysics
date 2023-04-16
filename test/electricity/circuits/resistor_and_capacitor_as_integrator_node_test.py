# Description
## Assert we have 3 Volts applied to 2-Ohm resistor in series with 2 Farads capacitor.
## After 1 Tau seconds capacitor voltage should be 63% of initial voltage. Tau = R * C = 4.

from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.electricity.circuits import resistor_and_capacitor_as_integrator_node as rc_node

@fixture
def test_args():
    V0 = Quantity(units.voltage, 3 * units.volt)
    R = Quantity(units.impedance, 2 * units.ohm)
    C = Quantity(units.capacitance, 2 * units.farad)
    Time = Quantity(units.time, 4 * units.second)
    Args = namedtuple("Args", ["V0", "R", "C", "Time"])
    return Args(V0=V0, R=R, C=C, Time=Time)

def test_basic_voltage(test_args):
    result = rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, test_args.R, test_args.Time)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result_voltage = convert_to(result, units.volt).subs(units.volt, 1).evalf(2)
    assert result_voltage == approx(1.89, 0.01)

def test_bad_voltage(test_args):
    Vb = Quantity(units.length)
    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(Vb, test_args.C, test_args.R, test_args.Time)
    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(100, test_args.C, test_args.R, test_args.Time)

def test_bad_capacity(test_args):
    Cb = Quantity(units.length)
    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(test_args.V0, Cb, test_args.R, test_args.Time)
    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(test_args.V0, 100, test_args.R, test_args.Time)

def test_bad_resistance(test_args):
    Rb = Quantity(units.length)
    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, Rb, test_args.Time)
    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, 100, test_args.Time)

def test_bad_time(test_args):
    Tb = Quantity(units.length)
    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, test_args.R, Tb)
    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, test_args.R, 100)
