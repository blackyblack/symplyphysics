from sympy import Eq, solve, Matrix
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)

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
parameter_11 = Symbol("parameter_11", units.impedance)
parameter_12 = Symbol("parameter_12", dimensionless)
parameter_21 = Symbol("parameter_21", dimensionless)
parameter_22 = Symbol("parameter_22", units.conductance)

law = Eq(Matrix([input_voltage, output_current]), Matrix([[parameter_11, parameter_12], [parameter_21, parameter_22]]) * Matrix([input_current, output_voltage]))


def print_law() -> str:
    return print_expression(law)


@validate_input(input_voltage_=input_voltage, output_current_=output_current,)# parameters_=parameter_11)
@validate_output((units.current, units.voltage))
def calculate_current_and_voltage(input_voltage_: Quantity, output_current_: Quantity, parameters_: tuple[tuple[Quantity, float], tuple[float, Quantity]]) -> tuple[Quantity, Quantity]:
    result = solve(law, [input_current, output_voltage], dict=True)[0]
    result_input_current = result[input_current]
    result_output_voltage = result[output_voltage]
    substitutions = {
        input_voltage: input_voltage_,
        output_current: output_current_,
        parameter_11: parameters_[0][0],
        parameter_12: parameters_[0][1],
        parameter_21: parameters_[1][0],
        parameter_22: parameters_[1][1],
    }
    result_input_current = result_input_current.subs(substitutions)
    result_output_voltage = result_output_voltage.subs(substitutions)
    return (Quantity(result_input_current), Quantity(result_output_voltage))
