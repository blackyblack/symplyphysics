from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits import total_resistance_matrix as resistance_matrix_law

## The input voltage is 100 volts, the output voltage is 25 volts. Z-parameters Z11, Z12, Z21, Z22 are equal
## to 100, 50, 50, 1000 ohms, respectively.
## Then the values of the input and output currents will be equal to 1.013 ampere and -25.64 milliampere, respectively.

Args = namedtuple("Args", ["input_voltage", "output_voltage", "impedances"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    input_voltage = Quantity(100 * units.volt)
    output_voltage = Quantity(25 * units.volt)
    impedances = ((Quantity(100 * units.ohm), Quantity(50 * units.ohm)), (Quantity(50 * units.ohm), Quantity(1000 * units.ohm)))
    return Args(input_voltage=input_voltage,
        output_voltage=output_voltage,
        impedances=impedances
        )


def test_basic_currents(test_args: Args) -> None:
    result = resistance_matrix_law.calculate_currents(test_args.input_voltage, test_args.output_voltage, test_args.impedances)
    assert_equal(result[0], 1.013 * units.ampere)
    assert_equal(result[1], -25.64 * prefixes.milli * units.ampere)


def test_bad_voltages(test_args: Args) -> None:
    bad_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_matrix_law.calculate_currents(bad_voltage, test_args.output_voltage, test_args.impedances)
    with raises(TypeError):
        resistance_matrix_law.calculate_currents(100, test_args.output_voltage, test_args.impedances)
    with raises(errors.UnitsError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, bad_voltage, test_args.impedances)
    with raises(TypeError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, 100, test_args.impedances)


def test_bad_impedances(test_args: Args) -> None:
    bad_impedance = Quantity(1 * units.coulomb)
    with raises(AssertionError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, test_args.output_voltage, ((bad_impedance, test_args.impedances[0][1]), (test_args.impedances[1][0], test_args.impedances[1][1])))
    with raises(AttributeError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, test_args.output_voltage, ((100, test_args.impedances[0][1]), (test_args.impedances[1][0], test_args.impedances[1][1])))
    with raises(AssertionError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, test_args.output_voltage, ((test_args.impedances[0][0], bad_impedance), (test_args.impedances[1][0], test_args.impedances[1][1])))
    with raises(AttributeError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, test_args.output_voltage, ((test_args.impedances[0][0], 100), (test_args.impedances[1][0], test_args.impedances[1][1])))
    with raises(AssertionError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, test_args.output_voltage, ((test_args.impedances[0][0], test_args.impedances[0][1]), (bad_impedance, test_args.impedances[1][1])))
    with raises(AttributeError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, test_args.output_voltage, ((test_args.impedances[0][0], test_args.impedances[0][1]), (100, test_args.impedances[1][1])))
    with raises(AssertionError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, test_args.output_voltage, ((test_args.impedances[0][0], test_args.impedances[0][1]), (test_args.impedances[1][0], bad_impedance)))
    with raises(AttributeError):
        resistance_matrix_law.calculate_currents(test_args.input_voltage, test_args.output_voltage, ((test_args.impedances[0][0], test_args.impedances[0][1]), (test_args.impedances[1][0], 100)))
