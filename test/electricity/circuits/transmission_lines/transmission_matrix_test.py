from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits.transmission_lines import transmission_matrix as parameter_matrix_law

## The input voltage is 100 volt, the input current is 25 ampere. ABCD-parameters: A is equal to 50, B is equal to 100 ohm,
## C is equal to 1000 siemens, D is equal to 50.
## Then the values of the output voltage and output current will be equal to -25.64 millivolt and 1.013 ampere, respectively.

Args = namedtuple("Args", ["input_voltage", "input_current", "parameters"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    input_voltage = Quantity(100 * units.volt)
    input_current = Quantity(25 * units.ampere)
    parameters = ((50, Quantity(100 * units.ohm)), (Quantity(1000 * units.siemens), 50))
    return Args(input_voltage=input_voltage, input_current=input_current, parameters=parameters)


def test_basic_current_and_voltage(test_args: Args) -> None:
    result = parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage,
        test_args.input_current, test_args.parameters)
    assert_equal(result[0], -25.64 * prefixes.milli * units.volt)
    assert_equal(result[1], 1.013 * units.ampere)


def test_bad_voltage(test_args: Args) -> None:
    bad_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        parameter_matrix_law.calculate_current_and_voltage(bad_voltage, test_args.input_current,
            test_args.parameters)
    with raises(TypeError):
        parameter_matrix_law.calculate_current_and_voltage(100, test_args.input_current,
            test_args.parameters)


def test_bad_current(test_args: Args) -> None:
    bad_current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, bad_current,
            test_args.parameters)
    with raises(TypeError):
        parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, 100,
            test_args.parameters)


def test_bad_parameters(test_args: Args) -> None:
    bad_parameter = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage,
            test_args.input_current, ((test_args.parameters[0][0], bad_parameter),
            (test_args.parameters[1][0], test_args.parameters[1][1])))
    with raises(TypeError):
        parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage,
            test_args.input_current, ((test_args.parameters[0][0], 100),
            (test_args.parameters[1][0], test_args.parameters[1][1])))
    with raises(errors.UnitsError):
        parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage,
            test_args.input_current, ((test_args.parameters[0][0], test_args.parameters[0][1]),
            (bad_parameter, test_args.parameters[1][1])))
    with raises(TypeError):
        parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage,
            test_args.input_current, ((test_args.parameters[0][0], test_args.parameters[0][1]),
            (100, test_args.parameters[1][1])))
