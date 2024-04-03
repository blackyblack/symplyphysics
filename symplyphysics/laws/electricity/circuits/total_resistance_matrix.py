from sympy import Eq, solve, Matrix
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    SI,
)

## Description
## The impedance matrix is one of the ways to describe a microwave device. The Z-parameters of the device act as elements
## of this matrix. The matrix equation relates the input and output voltages to the input and output currents.

## Law is: [U1, U2] = Matrix([[Z11, Z12], [Z21, Z22]]) * [I1, I2], where
## U1 - input voltage,
## U2 - output voltage,
## Z11 - ratio of input voltage to input current at idle at the output,
## Z12 - ratio of input voltage to output current at idle at the input,
## Z21 - ratio of output voltage to input current at idle at the output,
## Z22 - ratio of output voltage to output current at idle at the input,
## I1 - input current,
## I2 - output current.

input_voltage = Symbol("input_voltage", units.voltage)
output_voltage = Symbol("output_voltage", units.voltage)
input_current = Symbol("input_current", units.current)
output_current = Symbol("output_current", units.current)
impedance_input_input = Symbol("impedance_input_input", units.impedance)
impedance_input_output = Symbol("impedance_input_output", units.impedance)
impedance_output_input = Symbol("impedance_output_input", units.impedance)
impedance_output_output = Symbol("impedance_output_output", units.impedance)

law = Eq(Matrix([input_voltage, output_voltage]), Matrix([[impedance_input_input, impedance_input_output], [impedance_output_input, impedance_output_output]]) * Matrix([input_current, output_current]))


def print_law() -> str:
    return print_expression(law)


@validate_input(input_voltage_=input_voltage, output_voltage_=output_voltage,)
@validate_output(units.current)
def calculate_currents(input_voltage_: Quantity, output_voltage_: Quantity, impedances_: tuple[tuple[Quantity, Quantity], tuple[Quantity, Quantity]]) -> tuple[Quantity, Quantity]:
    for impedance_1, impedance_2 in impedances_:
        assert SI.get_dimension_system().equivalent_dims(impedance_1.dimension,
            units.impedance)
        assert SI.get_dimension_system().equivalent_dims(impedance_2.dimension,
            units.impedance)
    result_currents = solve(law, [input_current, output_current], dict=True)[0]
    result_input_current = result_currents[input_current]
    result_output_current = result_currents[output_current]
    substitutions = {
        input_voltage: input_voltage_,
        output_voltage: output_voltage_,
        impedance_input_input: impedances_[0][0],
        impedance_input_output: impedances_[0][1],
        impedance_output_input: impedances_[1][0],
        impedance_output_output: impedances_[1][1],
    }
    result_input_current = result_input_current.subs(substitutions)
    result_output_current = result_output_current.subs(substitutions)
    return (Quantity(result_input_current), Quantity(result_output_current))
