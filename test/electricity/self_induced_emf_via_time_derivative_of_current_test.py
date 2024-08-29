from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity, prefixes)
from symplyphysics.laws.electricity import self_induced_emf_via_time_derivative_of_current as self_induction_def

# Description
## Current through 2.5 millihenry inductor increases from 0 to 0.5A in 5 seconds.
## Self-induction voltage should be -0.00025 volts

Args = namedtuple("Args", ["L", "I0", "I1", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    L = Quantity(2.5 * prefixes.milli * units.henry)
    I0 = Quantity(0 * units.ampere)
    I1 = Quantity(0.5 * units.ampere)
    t = Quantity(5 * units.second)
    return Args(L=L, I0=I0, I1=I1, t=t)


def test_basic_voltage(test_args: Args) -> None:
    result = self_induction_def.calculate_voltage(test_args.L, test_args.I0, test_args.I1,
        test_args.t)
    assert_equal(result, -0.00025 * units.volt)


def test_voltage_with_bad_induction(test_args: Args) -> None:
    Lb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(Lb, test_args.I0, test_args.I1, test_args.t)
    with raises(TypeError):
        self_induction_def.calculate_voltage(100, test_args.I0, test_args.I1, test_args.t)


def test_voltage_with_bad_current(test_args: Args) -> None:
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(test_args.L, Ib, test_args.I1, test_args.t)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, Ib, test_args.t)
    with raises(TypeError):
        self_induction_def.calculate_voltage(test_args.L, 100, test_args.I1, test_args.t)
    with raises(TypeError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, 100, test_args.t)


def test_voltage_with_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, test_args.I1, tb)
    with raises(TypeError):
        self_induction_def.calculate_voltage(test_args.L, test_args.I0, test_args.I1, 100)
