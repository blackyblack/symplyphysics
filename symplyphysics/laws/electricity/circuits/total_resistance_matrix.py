from sympy import Eq, solve, Matrix
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

## Description
## The impedance matrix is one of the ways to describe a microwave device. The Z-parameters of the device act as elements
## of this matrix. The matrix equation relates the input and output voltages to the input and output currents.

## Law is: [U1, U2] = Matrix([[Z11, Z12], [Z21, Z22]]), where
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
impedance_11 = Symbol("impedance_11", units.impedance)
impedance_12 = Symbol("impedance_12", units.impedance)
impedance_21 = Symbol("impedance_21", units.impedance)
impedance_22 = Symbol("impedance_22", units.impedance)

law = Eq(Matrix([input_voltage, output_voltage]), Matrix([[impedance_11, impedance_12], [impedance_21, impedance_22]]) * Matrix([input_current, output_current]))


def print_law() -> str:
    return print_expression(law)


@validate_input(input_voltage_=input_voltage, output_voltage_=output_voltage, impedance_11_=impedance_11, impedance_12_=impedance_12, impedance_21_=impedance_21, impedance_22_=impedance_22)
@validate_output(units.current)
def calculate_currents(input_voltage_: Quantity, output_voltage_: Quantity, impedance_11_: Quantity,
    impedance_12_: Quantity, impedance_21_: Quantity, impedance_22_: Quantity) -> list[Quantity, Quantity]:
    result_input_current = solve(law, [input_current, output_current], dict=True)[0][input_current]
    result_output_current = solve(law, [input_current, output_current], dict=True)[0][output_current]
    result_input_current = result_input_current.subs({
        input_voltage: input_voltage_,
        output_voltage: output_voltage_,
        impedance_11: impedance_11_,
        impedance_12: impedance_12_,
        impedance_21: impedance_21_,
        impedance_22: impedance_22_,
    })
    result_output_current = result_output_current.subs({
        input_voltage: input_voltage_,
        output_voltage: output_voltage_,
        impedance_11: impedance_11_,
        impedance_12: impedance_12_,
        impedance_21: impedance_21_,
        impedance_22: impedance_22_,
    })
    return [Quantity(result_input_current), Quantity(result_output_current)]
