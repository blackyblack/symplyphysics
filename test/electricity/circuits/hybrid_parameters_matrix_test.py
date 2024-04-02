from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits import hybrid_parameters_matrix as parameter_matrix_law

## The input voltage is 100 volt, the output current is 25 ampere. H-parameters: H11 is equal to 100 ohm, H12 is equal to 50,
## H21 is equal to 50, H22 is equal to 1000 siemens.
## Then the values of the input current and output voltage will be equal to 1.013 ampere and -25.64 millivolt, respectively.

Args = namedtuple("Args", ["input_voltage", "output_current", "parameter_11", "parameter_12", "parameter_21", "parameter_22"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    input_voltage = Quantity(100 * units.volt)
    output_current = Quantity(25 * units.ampere)
    parameter_11 = Quantity(100 * units.ohm)
    parameter_12 = 50
    parameter_21 = 50
    parameter_22 = Quantity(1000 * units.siemens)
    return Args(input_voltage=input_voltage,
        output_current=output_current,
        parameter_11=parameter_11,
        parameter_12=parameter_12,
        parameter_21=parameter_21,
        parameter_22=parameter_22)


def test_basic_current_and_voltage(test_args: Args) -> None:
    result = parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, test_args.output_current, ((test_args.parameter_11, test_args.parameter_12), (test_args.parameter_21, test_args.parameter_22)))
    assert_equal(result[0], 1.013 * units.ampere)
    assert_equal(result[1], -25.64 * prefixes.milli * units.volt)


def test_bad_voltage(test_args: Args) -> None:
    bad_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        parameter_matrix_law.calculate_current_and_voltage(bad_voltage, test_args.output_current, ((test_args.parameter_11, test_args.parameter_12), (test_args.parameter_21, test_args.parameter_22)))
    with raises(TypeError):
        parameter_matrix_law.calculate_current_and_voltage(100, test_args.output_current, ((test_args.parameter_11, test_args.parameter_12), (test_args.parameter_21, test_args.parameter_22)))


def test_bad_current(test_args: Args) -> None:
    bad_current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, bad_current, ((test_args.parameter_11, test_args.parameter_12), (test_args.parameter_21, test_args.parameter_22)))
    with raises(TypeError):
        parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, 100, ((test_args.parameter_11, test_args.parameter_12), (test_args.parameter_21, test_args.parameter_22)))


# def test_bad_parametrs(test_args: Args) -> None:
#     bad_parameter = Quantity(1 * units.coulomb)
#     with raises(errors.UnitsError):
#         parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, test_args.output_current, ((bad_parameter, test_args.parameter_12), (test_args.parameter_21, test_args.parameter_22)))
#     with raises(TypeError):
#         parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, test_args.output_current, ((100, test_args.parameter_12), (test_args.parameter_21, test_args.parameter_22)))
#     with raises(errors.UnitsError):
#         parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, test_args.output_current, ((test_args.parameter_11, bad_parameter), (test_args.parameter_21, test_args.parameter_22)))
#     with raises(errors.UnitsError):
#         parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, test_args.output_current, ((test_args.parameter_11, test_args.parameter_12), (bad_parameter, test_args.parameter_22)))
#     with raises(errors.UnitsError):
#         parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, test_args.output_current, ((test_args.parameter_11, test_args.parameter_12), (test_args.parameter_21, bad_parameter)))
#     with raises(TypeError):
#         parameter_matrix_law.calculate_current_and_voltage(test_args.input_voltage, test_args.output_current, ((test_args.parameter_11, test_args.parameter_12), (test_args.parameter_21, 100)))
