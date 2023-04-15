from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.definitions import self_induction_voltage_is_current_derivative as self_induction_def

# Description
## Current through 0.0025 henry inductor increases from 0 to 0.5A in 5 seconds. 
## Self- induction voltage should be -0.00025 volts

@fixture
def test_args():
    L = Quantity(units.inductance, 0.0025 * units.henry)
    I0 = Quantity(units.current, 0 * units.ampere)
    I1 = Quantity(units.current, 0.5 * units.ampere)
    t = Quantity(units.time, 5 * units.second)
    Args = namedtuple("Args", ["L", "I0", "I1", "t"])
    return Args(L=L, I0=I0, I1=I1, t=t)


def test_basic_voltage(test_args):
    result = self_induction_def.calculate_voltage(
        test_args.L, test_args.I0, test_args.I1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(
        result.dimension, units.voltage)
    result_current = convert_to(result, self_induction_def.definition_units_SI).subs(units.volt, 1).evalf(6)
    assert result_current == approx(-0.00025, 0.000001)

def test_voltage_with_bad_induction(test_args):
    Lb = Quantity(units.length)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(Lb, test_args.I0, test_args.I1, test_args.t)
    with raises(TypeError):
        self_induction_def.calculate_voltage(100, test_args.I0, test_args.I1, test_args.t)

def test_voltage_with_bad_current(test_args):
    Ib = Quantity(units.length)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(test_args.L, Ib, test_args.I1, test_args.t)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, Ib, test_args.t)
    with raises(TypeError):
        self_induction_def.calculate_voltage(test_args.L, 100, test_args.I1, test_args.t)
    with raises(TypeError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, 100, test_args.t)

def test_voltage_with_bad_time(test_args):
    tb = Quantity(units.length)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, test_args.I1, tb)
    with raises(TypeError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, test_args.I1, 100)
