from collections import namedtuple
from pytest import fixture, raises
from sympy import I
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import serial_impedance as serial_resistor

# Description
## Assert we have two elements of circuit with (1 + 1 * I) Ohm and (2 + 2 * I) Ohm impedances.
## Resulting impedance should be 3 Ohm.

Args = namedtuple("Args", ["R1", "R2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    R1 = Quantity((1 + 1 * I) * units.ohm)
    R2 = Quantity((2 + 2 * I) * units.ohm)
    return Args(R1=R1, R2=R2)


def test_basic_impedance(test_args: Args) -> None:
    result = serial_resistor.calculate_serial_impedance([test_args.R1, test_args.R2])
    assert_equal(result, (3 + 3 * I) * units.ohm)


def test_three_resistors_array(test_args: Args) -> None:
    R3 = Quantity(3 * units.ohm)
    result = serial_resistor.calculate_serial_impedance([test_args.R1, test_args.R2, R3])
    assert_equal(result, (6 + 3 * I) * units.ohm)


def test_bad_impedance(test_args: Args) -> None:
    Zb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        serial_resistor.calculate_serial_impedance([Zb, test_args.R2])
    with raises(TypeError):
        serial_resistor.calculate_serial_impedance([100, test_args.R2])
    with raises(errors.UnitsError):
        serial_resistor.calculate_serial_impedance([test_args.R1, Zb])
    with raises(TypeError):
        serial_resistor.calculate_serial_impedance([test_args.R1, 100])
    with raises(errors.UnitsError):
        serial_resistor.calculate_serial_impedance([Zb, Zb])
    with raises(TypeError):
        serial_resistor.calculate_serial_impedance([100, 100])
    with raises(TypeError):
        serial_resistor.calculate_serial_impedance(test_args.R1)
