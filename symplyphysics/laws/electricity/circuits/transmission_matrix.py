from sympy import Eq, solve, Matrix
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    dimensionless,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

## Description
## The transmission parameters matrix is one of the ways to describe a microwave device. The ABCD-parameters of the device act as elements
## of this matrix. The matrix equation relates the input voltage and input current to the output voltage and output current.

## Law is: [U1, I1] = Matrix([[A, B], [C, D]]) * [U2, I2], where
## U1 - input voltage,
## U2 - output voltage,
## A - ratio of input voltage to output voltage at idle at the output,
## B - ratio of input voltage to output current in case of a short circuit at the output,
## C - ratio of input current to output voltage at idle at the output,
## D - ratio of input current to output current in case of a short circuit at the output,
## I1 - input current,
## I2 - output current.

input_voltage = Symbol("input_voltage", units.voltage)
output_voltage = Symbol("output_voltage", units.voltage)
input_current = Symbol("input_current", units.current)
output_current = Symbol("output_current", units.current)
parameter_voltage_to_voltage = Symbol("parameter_voltage_to_voltage", dimensionless)
parameter_impedance = Symbol("parameter_impedance", units.impedance)
parameter_conductance = Symbol("parameter_conductance", units.conductance)
parameter_current_to_current = Symbol("parameter_current_to_current", dimensionless)

law = Eq(Matrix([input_voltage, input_current]), Matrix([[parameter_voltage_to_voltage, parameter_impedance], [parameter_conductance, parameter_current_to_current]]) * Matrix([output_voltage, output_current]))


def print_law() -> str:
    return print_expression(law)


@validate_input(input_voltage_=input_voltage, input_current_=input_current,)
def calculate_current_and_voltage(input_voltage_: Quantity, input_current_: Quantity, parameters_: tuple[tuple[float, Quantity], tuple[Quantity, float]]) -> tuple[Quantity, Quantity]:
    assert_equivalent_dimension(parameters_[0][1], f"parameters_[{0}][{1}]", "calculate_current_and_voltage", units.impedance)
    assert_equivalent_dimension(parameters_[1][0], f"parameters_[{1}][{0}]", "calculate_current_and_voltage", units.conductance)
    result = solve(law, [output_voltage, output_current], dict=True)[0]
    result_output_current = result[output_current]
    result_output_voltage = result[output_voltage]
    substitutions = {
        input_voltage: input_voltage_,
        input_current: input_current_,
        parameter_voltage_to_voltage: parameters_[0][0],
        parameter_impedance: parameters_[0][1],
        parameter_conductance: parameters_[1][0],
        parameter_current_to_current: parameters_[1][1],
    }
    result_output_current = Quantity(result_output_current.subs(substitutions))
    result_output_voltage = Quantity(result_output_voltage.subs(substitutions))
    assert_equivalent_dimension(result_output_current, 'result_output_current', "calculate_current_and_voltage", units.current)
    assert_equivalent_dimension(result_output_voltage, 'result_output_voltage', "calculate_current_and_voltage", units.voltage)
    return (result_output_voltage, result_output_current)
