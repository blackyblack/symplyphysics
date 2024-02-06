from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import resistivity_of_serial_resistors as serial_resistor

# Description
## Assert we have two resistors with 1 Ohm and 2 Ohm resistances.
## Resulting resistance should be 3 Ohm.


@fixture(name="test_args")
def test_args_fixture():
    R1 = Quantity(1 * units.ohm)
    R2 = Quantity(2 * units.ohm)
    Args = namedtuple("Args", ["R1", "R2"])
    return Args(R1=R1, R2=R2)


def test_basic_resistance(test_args):
    result = serial_resistor.calculate_serial_resistance([test_args.R1, test_args.R2])
    assert_equal(result, 3 * units.ohm)


def test_three_resistors_array(test_args):
    R3 = Quantity(3 * units.ohm)
    result = serial_resistor.calculate_serial_resistance([test_args.R1, test_args.R2, R3])
    assert_equal(result, 6 * units.ohm)


def test_bad_resistance(test_args):
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
