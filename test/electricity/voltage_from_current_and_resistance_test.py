from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import voltage_from_current_and_resistance as voltage_law

# Description
## With an internal resistance of 20 ohms and an external resistance of 100 ohms,
## the voltage will be 1.2 volts with a current of 0.01 amperes.
## https://www.chipdip.ru/calc/ohm-law?i=0.01&r=120

Args = namedtuple("Args", ["current", "inner_resistance", "outer_resistance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    current = Quantity(0.01 * units.ampere)
    inner_resistance = Quantity(20 * units.ohm)
    outer_resistance = Quantity(100 * units.ohm)
    return Args(current=current,
        inner_resistance=inner_resistance,
        outer_resistance=outer_resistance)


def test_basic_voltage(test_args: Args) -> None:
    result = voltage_law.calculate_voltage(test_args.current, test_args.inner_resistance,
        test_args.outer_resistance)
    assert_equal(result, 1.2 * units.volt)


def test_bad_current(test_args: Args) -> None:
    current = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage(current, test_args.inner_resistance,
            test_args.outer_resistance)
    with raises(TypeError):
        voltage_law.calculate_voltage(100, test_args.inner_resistance, test_args.outer_resistance)


def test_bad_resistance(test_args: Args) -> None:
    resistance = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage(test_args.current, resistance, test_args.outer_resistance)
    with raises(TypeError):
        voltage_law.calculate_voltage(test_args.current, 100, test_args.outer_resistance)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage(test_args.current, test_args.inner_resistance, resistance)
    with raises(TypeError):
        voltage_law.calculate_voltage(test_args.current, test_args.inner_resistance, 100)
