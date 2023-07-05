from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, convert_to, Quantity, SI, prefixes)
from symplyphysics.definitions import self_induction_voltage_is_current_derivative as self_induction_def

# Description
## Current through 2.5 millihenry inductor increases from 0 to 0.5A in 5 seconds.
## Self-induction voltage should be -0.00025 volts


@fixture(name="test_args")
def test_args_fixture():
    L = Quantity(2.5 * prefixes.milli * units.henry)
    I0 = Quantity(0 * units.ampere)
    I1 = Quantity(0.5 * units.ampere)
    t = Quantity(5 * units.second)
    Args = namedtuple("Args", ["L", "I0", "I1", "t"])
    return Args(L=L, I0=I0, I1=I1, t=t)


def test_basic_voltage(test_args):
    result = self_induction_def.calculate_voltage(test_args.L, test_args.I0, test_args.I1,
        test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result_current = convert_to(result, self_induction_def.definition_units_SI).evalf(6)
    assert result_current == approx(-0.00025, 0.000001)


def test_voltage_with_bad_induction(test_args):
    Lb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(Lb, test_args.I0, test_args.I1, test_args.t)
    with raises(TypeError):
        self_induction_def.calculate_voltage(100, test_args.I0, test_args.I1, test_args.t)


def test_voltage_with_bad_current(test_args):
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(test_args.L, Ib, test_args.I1, test_args.t)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, Ib, test_args.t)
    with raises(TypeError):
        self_induction_def.calculate_voltage(test_args.L, 100, test_args.I1, test_args.t)
    with raises(TypeError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, 100, test_args.t)


def test_voltage_with_bad_time(test_args):
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, test_args.I1, tb)
    with raises(TypeError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, test_args.I1, 100)
