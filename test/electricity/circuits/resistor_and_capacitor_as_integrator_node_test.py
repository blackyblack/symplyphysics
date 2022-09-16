# Description
## Assert we have 3 Volts applied to 2-Ohm resistor in series with 2 Farads capacitor.
## After 1 Tau seconds capacitor voltage should be 63% of initial voltage. Tau = R * C = 4.

from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.electricity.circuits import resistor_and_capacitor_as_integrator_node as rc_node

@fixture
def test_args():
    V0 = units.Quantity('V0')
    R = units.Quantity('R')
    C = units.Quantity('C')
    Time = units.Quantity('Time')
    SI.set_quantity_dimension(V0, units.voltage)
    SI.set_quantity_scale_factor(V0, 3 * units.volt)
    SI.set_quantity_dimension(R, units.impedance)
    SI.set_quantity_scale_factor(R, 2 * units.ohm)
    SI.set_quantity_dimension(C, units.capacitance)
    SI.set_quantity_scale_factor(C, 2 * units.farad)
    SI.set_quantity_dimension(Time, units.time)
    SI.set_quantity_scale_factor(Time, 4 * units.second)

    Args = namedtuple('Args', ['V0', 'R', 'C', 'Time'])
    return Args(V0 = V0, R = R, C = C, Time = Time)

def test_basic_voltage(test_args):
    result = rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, test_args.R, test_args.Time)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result_voltage = convert_to(result, units.volt).subs(units.volt, 1).evalf(2)
    assert result_voltage == approx(1.89, 0.01)

def test_bad_voltage(test_args):
    Vb = units.Quantity('Vb')
    SI.set_quantity_dimension(Vb, units.length)
    SI.set_quantity_scale_factor(Vb, 1 * units.meter)

    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(Vb, test_args.C, test_args.R, test_args.Time)

    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(100, test_args.C, test_args.R, test_args.Time)

def test_bad_capacity(test_args):
    Cb = units.Quantity('Cb')
    SI.set_quantity_dimension(Cb, units.length)
    SI.set_quantity_scale_factor(Cb, 1 * units.meter)

    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(test_args.V0, Cb, test_args.R, test_args.Time)

    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(test_args.V0, 100, test_args.R, test_args.Time)

def test_bad_resistance(test_args):
    Rb = units.Quantity('Rb')
    SI.set_quantity_dimension(Rb, units.length)
    SI.set_quantity_scale_factor(Rb, 1 * units.meter)

    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, Rb, test_args.Time)

    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, 100, test_args.Time)

def test_bad_time(test_args):
    Tb = units.Quantity('Tb')
    SI.set_quantity_dimension(Tb, units.length)
    SI.set_quantity_scale_factor(Tb, 1 * units.meter)

    with raises(errors.UnitsError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, test_args.R, Tb)

    with raises(TypeError):
        rc_node.calculate_capacitor_voltage(test_args.V0, test_args.C, test_args.R, 100)