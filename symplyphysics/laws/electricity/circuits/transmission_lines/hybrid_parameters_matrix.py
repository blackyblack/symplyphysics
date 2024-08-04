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
## The hybrid parameters matrix is one of the ways to describe a microwave device. The H-parameters of the device act as elements
## of this matrix. The matrix equation relates the input voltage and output current to the input current and output voltage.

## Law is: [U1, I2] = Matrix([[H11, H12], [H21, H22]]) * [I1, U2], where
## U1 - input voltage,
## U2 - output voltage,
## H11 - ratio of input voltage to input current in case of a short circuit at the output,
## H12 - ratio of input voltage to output voltage at idle at the input,
## H21 - ratio of output current to input current in case of a short circuit at the output,
## H22 - ratio of output current to output voltage at idle at the input,
## I1 - input current,
## I2 - output current.

input_voltage = Symbol("input_voltage", units.voltage)
output_voltage = Symbol("output_voltage", units.voltage)
input_current = Symbol("input_current", units.current)
output_current = Symbol("output_current", units.current)
parameter_input_input = Symbol("parameter_input_input", units.impedance)
parameter_input_output = Symbol("parameter_input_output", dimensionless)
parameter_output_input = Symbol("parameter_output_input", dimensionless)
parameter_output_output = Symbol("parameter_output_output", units.conductance)

law = Eq(
    Matrix([input_voltage, output_current]),
    Matrix([[parameter_input_input, parameter_input_output],
    [parameter_output_input, parameter_output_output]]) * Matrix([input_current, output_voltage]))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    input_voltage_=input_voltage,
    output_current_=output_current,
)
def calculate_current_and_voltage(
    input_voltage_: Quantity, output_current_: Quantity, parameters_: tuple[tuple[Quantity, float],
    tuple[float, Quantity]]) -> tuple[Quantity, Quantity]:
    assert_equivalent_dimension(parameters_[0][0], "parameters_[0][0]",
        "calculate_current_and_voltage", units.impedance)
    assert_equivalent_dimension(parameters_[1][1], "parameters_[1][1]",
        "calculate_current_and_voltage", units.conductance)
    result = solve(law, [input_current, output_voltage], dict=True)[0]
    result_input_current = result[input_current]
    result_output_voltage = result[output_voltage]
    substitutions = {
        input_voltage: input_voltage_,
        output_current: output_current_,
        parameter_input_input: parameters_[0][0],
        parameter_input_output: parameters_[0][1],
        parameter_output_input: parameters_[1][0],
        parameter_output_output: parameters_[1][1],
    }
    result_input_current = Quantity(result_input_current.subs(substitutions))
    result_output_voltage = Quantity(result_output_voltage.subs(substitutions))
    assert_equivalent_dimension(result_input_current, 'result_input_current',
        "calculate_current_and_voltage", units.current)
    assert_equivalent_dimension(result_output_voltage, 'result_output_voltage',
        "calculate_current_and_voltage", units.voltage)
    return (result_input_current, result_output_voltage)
