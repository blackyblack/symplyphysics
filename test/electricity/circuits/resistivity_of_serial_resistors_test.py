from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity.circuits import resistivity_of_serial_resistors as serial_resistor

# Description
## Assert we have two resistors with 1 Ohm and 2 Ohm resistances.
## Resulting conductance should be 3 Ohm.

@fixture(name="test_args")
def test_args_fixture():
    R1 = Quantity(1 * units.ohm)
    R2 = Quantity(2 * units.ohm)
    Args = namedtuple("Args", ["R1", "R2"])
    return Args(R1=R1, R2=R2)

def test_basic_conductivity(test_args):
    result = serial_resistor.calculate_serial_resistance([test_args.R1, test_args.R2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_resistance = convert_to(result, units.ohm).evalf(3)
    assert result_resistance == approx(3, 0.001)

def test_three_resistors_array(test_args):
    S3 = Quantity(3 * units.ohm)
    result = serial_resistor.calculate_serial_resistance([test_args.R1, test_args.R2, S3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_resistance = convert_to(result, units.ohm).evalf(3)
    assert result_resistance == approx(6, 0.01)


def test_bad_conductivity(test_args):
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

