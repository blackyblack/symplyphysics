from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import power_via_voltage_and_current as power_law

# What is the power can be released by the battery at 9 V and a current of 0.5 A?

Args = namedtuple("Args", ["U", "I"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    I = Quantity(0.5 * units.ampere)
    U = Quantity(9 * units.volt)
    return Args(I=I, U=U)


def test_basic_power(test_args: Args) -> None:
    result = power_law.calculate_power(test_args.I, test_args.U)
    assert_equal(result, 4.5 * units.watt)


def test_bad_current(test_args: Args) -> None:
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_law.calculate_power(Ib, test_args.U)
    with raises(TypeError):
        power_law.calculate_power(100, test_args.U)


def test_bad_voltage(test_args: Args) -> None:
    Ub = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        power_law.calculate_power(test_args.I, Ub)
    with raises(TypeError):
        power_law.calculate_power(test_args.I, 100)
