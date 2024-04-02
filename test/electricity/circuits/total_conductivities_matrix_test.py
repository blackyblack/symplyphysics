from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits import total_conductivities_matrix as conductivity_matrix_law

## The input current is 100 ampere, the output current is 25 ampere. Y-parameters Y11, Y12, Y21, Y22 are equal
## to 100, 50, 50, 1000 siemens, respectively.
## Then the values of the input and output voltages will be equal to 1.013 volt and -25.64 millivolt, respectively.

Args = namedtuple("Args", ["input_current", "output_current", "conductivity_11", "conductivity_12", "conductivity_21", "conductivity_22"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    input_current = Quantity(100 * units.ampere)
    output_current = Quantity(25 * units.ampere)
    conductivity_11 = Quantity(100 * units.siemens)
    conductivity_12 = Quantity(50 * units.siemens)
    conductivity_21 = Quantity(50 * units.siemens)
    conductivity_22 = Quantity(1000 * units.siemens)
    return Args(input_current=input_current,
        output_current=output_current,
        conductivity_11=conductivity_11,
        conductivity_12=conductivity_12,
        conductivity_21=conductivity_21,
        conductivity_22=conductivity_22)


def test_basic_voltages(test_args: Args) -> None:
    result = conductivity_matrix_law.calculate_voltages(test_args.input_current, test_args.output_current, ((test_args.conductivity_11, test_args.conductivity_12), (test_args.conductivity_21, test_args.conductivity_22)))
    assert_equal(result[0], 1.013 * units.volt)
    assert_equal(result[1], -25.64 * prefixes.milli * units.volt)


def test_bad_voltages(test_args: Args) -> None:
    bad_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        conductivity_matrix_law.calculate_voltages(bad_voltage, test_args.output_current, ((test_args.conductivity_11, test_args.conductivity_12), (test_args.conductivity_21, test_args.conductivity_22)))
    with raises(TypeError):
        conductivity_matrix_law.calculate_voltages(100, test_args.output_current, ((test_args.conductivity_11, test_args.conductivity_12), (test_args.conductivity_21, test_args.conductivity_22)))
    with raises(errors.UnitsError):
        conductivity_matrix_law.calculate_voltages(test_args.input_current, bad_voltage, ((test_args.conductivity_11, test_args.conductivity_12), (test_args.conductivity_21, test_args.conductivity_22)))
    with raises(TypeError):
        conductivity_matrix_law.calculate_voltages(test_args.input_current, 100, ((test_args.conductivity_11, test_args.conductivity_12), (test_args.conductivity_21, test_args.conductivity_22)))


# def test_bad_conductivities(test_args: Args) -> None:
#     bad_conductivity = Quantity(1 * units.coulomb)
#     with raises(errors.UnitsError):
#         conductivity_matrix_law.calculate_voltages(test_args.input_current, test_args.output_current, ((bad_conductivity, test_args.conductivity_12), (test_args.conductivity_21, test_args.conductivity_22)))
#     with raises(TypeError):
#         conductivity_matrix_law.calculate_voltages(test_args.input_current, test_args.output_current, ((100, test_args.conductivity_12), (test_args.conductivity_21, test_args.conductivity_22)))
#     with raises(errors.UnitsError):
#         conductivity_matrix_law.calculate_voltages(test_args.input_current, test_args.output_current, ((test_args.conductivity_11, bad_conductivity), (test_args.conductivity_21, test_args.conductivity_22)))
#     with raises(TypeError):
#         conductivity_matrix_law.calculate_voltages(test_args.input_current, test_args.output_current, ((test_args.conductivity_11, 100), (test_args.conductivity_21, test_args.conductivity_22)))
#     with raises(errors.UnitsError):
#         conductivity_matrix_law.calculate_voltages(test_args.input_current, test_args.output_current, ((test_args.conductivity_11, test_args.conductivity_12), (bad_conductivity, test_args.conductivity_22)))
#     with raises(TypeError):
#         conductivity_matrix_law.calculate_voltages(test_args.input_current, test_args.output_current, ((test_args.conductivity_11, test_args.conductivity_12), (100, test_args.conductivity_22)))
#     with raises(errors.UnitsError):
#         conductivity_matrix_law.calculate_voltages(test_args.input_current, test_args.output_current, ((test_args.conductivity_11, test_args.conductivity_12), (test_args.conductivity_21, bad_conductivity)))
#     with raises(TypeError):
#         conductivity_matrix_law.calculate_voltages(test_args.input_current, test_args.output_current, ((test_args.conductivity_11, test_args.conductivity_12), (test_args.conductivity_21, 100)))
