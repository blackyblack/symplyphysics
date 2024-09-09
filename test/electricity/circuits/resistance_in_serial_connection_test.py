from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import resistance_in_serial_connection as serial_resistor

# Description
## Assert we have two resistors with 1 Ohm and 2 Ohm resistances.
## Resulting resistance should be 3 Ohm.

Args = namedtuple("Args", ["R1", "R2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    R1 = Quantity(1 * units.ohm)
    R2 = Quantity(2 * units.ohm)
    return Args(R1=R1, R2=R2)


def test_basic_resistance(test_args: Args) -> None:
    result = serial_resistor.calculate_serial_resistance([test_args.R1, test_args.R2])
    assert_equal(result, 3 * units.ohm)


def test_three_resistors_array(test_args: Args) -> None:
    R3 = Quantity(3 * units.ohm)
    result = serial_resistor.calculate_serial_resistance([test_args.R1, test_args.R2, R3])
    assert_equal(result, 6 * units.ohm)


def test_bad_resistance(test_args: Args) -> None:
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        serial_resistor.calculate_serial_resistance([Rb, test_args.R2])
    with raises(TypeError):
        serial_resistor.calculate_serial_resistance([100, test_args.R2])
    with raises(errors.UnitsError):
        serial_resistor.calculate_serial_resistance([test_args.R1, Rb])
    with raises(TypeError):
        serial_resistor.calculate_serial_resistance([test_args.R1, 100])
    with raises(errors.UnitsError):
        serial_resistor.calculate_serial_resistance([Rb, Rb])
    with raises(TypeError):
        serial_resistor.calculate_serial_resistance([100, 100])
    with raises(TypeError):
        serial_resistor.calculate_serial_resistance(test_args.R1)
