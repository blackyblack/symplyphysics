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
## The conductivity matrix is one of the ways to describe a microwave device. The Y-parameters of the device act as elements
## of this matrix. The matrix equation relates the input and output currents to the input and output voltages.

## Law is: [I1, I2] = Matrix([[Y11, Y12], [Y21, Y22]]) * [U1, U2], where
## U1 - input voltage,
## U2 - output voltage,
## Y11 - ratio of input current to input voltage in case of a short circuit at the output,
## Y12 - ratio of input current to output voltage in case of a short circuit at the output,
## Y21 - ratio of output current to input voltage in case of a short circuit at the output,
## Y22 - ratio of output current to output voltage in case of a short circuit at the output,
## I1 - input current,
## I2 - output current.

input_voltage = Symbol("input_voltage", units.voltage)
output_voltage = Symbol("output_voltage", units.voltage)
input_current = Symbol("input_current", units.current)
output_current = Symbol("output_current", units.current)
conductivity_11 = Symbol("conductivity_11", units.conductance)
conductivity_12 = Symbol("conductivity_12", units.conductance)
conductivity_21 = Symbol("conductivity_21", units.conductance)
conductivity_22 = Symbol("conductivity_22", units.conductance)

law = Eq(Matrix([input_current, output_current]), Matrix([[conductivity_11, conductivity_12], [conductivity_21, conductivity_22]]) * Matrix([input_voltage, output_voltage]))


def print_law() -> str:
    return print_expression(law)


@validate_input(input_current_=input_current, output_current_=output_current,) #conductivities_=conductivity_11)
@validate_output(units.voltage)
def calculate_voltages(input_current_: Quantity, output_current_: Quantity, conductivities_: tuple[tuple[Quantity, Quantity], tuple[Quantity, Quantity]]) -> tuple[Quantity, Quantity]:
    result_voltages = solve(law, [input_voltage, output_voltage], dict=True)[0]
    result_input_voltage = result_voltages[input_voltage]
    result_output_voltage = result_voltages[output_voltage]
    substitutions = {
        input_current: input_current_,
        output_current: output_current_,
        conductivity_11: conductivities_[0][0],
        conductivity_12: conductivities_[0][1],
        conductivity_21: conductivities_[1][0],
        conductivity_22: conductivities_[1][1],
    }
    result_input_voltage = result_input_voltage.subs(substitutions)
    result_output_voltage = result_output_voltage.subs(substitutions)
    return (Quantity(result_input_voltage), Quantity(result_output_voltage))
