from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.definitions import self_induction_voltage_is_current_derivative as self_induction_def

# Description
## Current through 0.0025 henry inductor increases from 0 to 0.5A in 5 seconds. 
## Self- induction voltage should be 0.00025 volts

@fixture
def test_args():
    L = units.Quantity('L')
    SI.set_quantity_dimension(L, units.inductance)
    SI.set_quantity_scale_factor(L, 0.0025 * units.henry)

    I0 = units.Quantity('I0')
    SI.set_quantity_dimension(I0, units.current)
    SI.set_quantity_scale_factor(I0, 0 * units.ampere)

    I1 = units.Quantity('I1')
    SI.set_quantity_dimension(I1, units.current)
    SI.set_quantity_scale_factor(I1, 0.5 * units.ampere)

    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.time)
    SI.set_quantity_scale_factor(t, 5 * units.second)

    Args = namedtuple('Args', ['L', 'I0', 'I1', 't'])
    return Args(L=L, I0=I0, I1=I1, t=t)


def test_basic_voltage(test_args):
    result = self_induction_def.calculate_voltage(
        test_args.L, test_args.I0, test_args.I1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(
        result.dimension, units.voltage)

    result_current = convert_to(result, self_induction_def.definition_dimension_SI).subs({
        units.voltage: 1}).evalf(6)
    assert result_current == approx(0.00025, 0.000001)

'''
def test_current_with_bad_charge(test_args):
    Q0b = units.Quantity('Q0b')
    SI.set_quantity_dimension(Q0b, units.length)
    SI.set_quantity_scale_factor(Q0b, 1 * units.meter)
    with raises(errors.UnitsError):
        current_def.calculate_current(
            Q0b, test_args.Q1, test_args.t)

    # Let Q1 be invalid    
    with raises(errors.UnitsError):
        current_def.calculate_current(
            test_args.Q0, Q0b, test_args.t)

    with raises(TypeError):
        current_def.calculate_current(
            100, test_args.Q1, test_args.t)

    with raises(TypeError):
        current_def.calculate_current(
            test_args.Q0, 100, test_args.t)


def test_current_with_bad_time(test_args):
    tb = units.Quantity('tb')
    SI.set_quantity_dimension(tb, units.length)
    SI.set_quantity_scale_factor(tb, 1 * units.meter)

    with raises(errors.UnitsError):
        current_def.calculate_current(
            test_args.Q0, test_args.Q1, tb)

    with raises(TypeError):
        current_def.calculate_current(
            test_args.Q0, test_args.Q1, 100)

'''
